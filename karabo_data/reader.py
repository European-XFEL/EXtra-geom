# coding: utf-8
"""
Collection of classes and functions to help reading HDF5 file generated at
The European XFEL.

Copyright (c) 2017, European X-Ray Free-Electron Laser Facility GmbH
All rights reserved.

You should have received a copy of the 3-Clause BSD License along with this
program. If not, see <https://opensource.org/licenses/BSD-3-Clause>
"""

from collections import defaultdict
import datetime
import fnmatch
from glob import glob
import h5py
import numpy as np
import os.path as osp
import pandas as pd
import re
import xarray


__all__ = ['H5File', 'RunDirectory', 'RunHandler', 'stack_data',
           'stack_detector_data', 'by_id', 'by_index', 'SourceNameError',
           'PropertyNameError',
          ]


RUN_DATA = 'RUN'
INDEX_DATA = 'INDEX'
METADATA = 'METADATA'

DETECTOR_NAMES = {'AGIPD', 'LPD'}
DETECTOR_SOURCE_RE = re.compile(r'(.+)/DET/(\d+)CH')


class FilenameInfo:
    is_detector = False
    detector_name = None
    detector_moduleno = -1

    _rawcorr_descr = {'RAW': 'Raw', 'CORR': 'Corrected'}

    def __init__(self, path):
        self.basename = osp.basename(path)
        nameparts = self.basename[:-3].split('-')
        assert len(nameparts) == 4, self.basename
        rawcorr, runno, datasrc, segment = nameparts
        m = re.match(r'([A-Z]+)(\d+)', datasrc)

        if m and m.group(1) == 'DA':
            self.description = "Aggregated data"
        elif m and m.group(1) in DETECTOR_NAMES:
            self.is_detector = True
            name, moduleno = m.groups()
            self.detector_name = name
            self.detector_moduleno = moduleno
            self.description = "{} detector data from {} module {}".format(
                self._rawcorr_descr.get(rawcorr, '?'), name, moduleno
            )
        else:
            self.description = "Unknown data source ({})", datasrc

class _SliceConstructor(type):
    """Allows instantiation like subclass[1:5]
    """
    def __getitem__(self, item):
        return self(item)

class by_id(metaclass=_SliceConstructor):
    def __init__(self, value):
        self.value = value

class by_index(metaclass=_SliceConstructor):
    def __init__(self, value):
        self.value = value

def _tid_to_slice_ix(tid, train_ids, stop=False):
    """Convert a train ID to an integer index for slicing the dataset

    Throws ValueError if the slice won't overlap the trains in the data.
    The *stop* parameter tells it which end of the slice it is making.
    """
    if tid is None:
        return None

    try:
        return train_ids.index(tid)
    except ValueError:
        pass

    if tid < train_ids[0]:
        if stop:
            raise ValueError("Train ID {} is before this run (starts at {})"
                             .format(tid, train_ids[0]))
        else:
            return None
    elif tid > train_ids[-1]:
        if stop:
            return None
        else:
            raise ValueError("Train ID {} is after this run (ends at {})"
                             .format(tid, train_ids[-1]))
    else:
        # This train ID is within the run, but doesn't have an entry.
        # Find the first ID in the run greater than the one given.
        return (train_ids > tid).nonzero()[0][0]

class SourceNameError(KeyError):
    def __init__(self, source, run=True):
        self.source = source
        self.run = run

    def __str__(self):
        run_file = 'run' if self.run else 'file'
        return "This {0} has no source named {1!r}.\n" \
               "See {0}.all_sources for available sources.".format(
            run_file, self.source
        )

class PropertyNameError(KeyError):
    def __init__(self, prop, source):
        self.prop = prop
        self.source = source

    def __str__(self):
        return "No property {!r} for source {!r}".format(self.prop, self.source)

class FileAccess:
    """Accesses a single Karabo HDF5 file

    Parameters
    ----------
    file: h5py.File
        Open h5py file object
    """

    def __init__(self, file):
        self.file = file
        self.filename = file.filename
        tid_data = file['INDEX/trainId'].value
        self.train_ids = tid_data[tid_data != 0]

        self.control_sources = set()
        self.instrument_sources = set()

        for source in file['METADATA/dataSourceId'].value:
            if not source:
                continue
            source = source.decode()
            category, _, h5_source = source.partition('/')
            if category == 'INSTRUMENT':
                device, _, chan_grp = h5_source.partition(':')
                chan, _, group = chan_grp.partition('/')
                source = device + ':' + chan
                self.instrument_sources.add(source)
                # TODO: Do something with groups?
            elif category == 'CONTROL':
                self.control_sources.add(h5_source)
            else:
                raise ValueError("Unknown data category %r" % category)

        # {(file, source, group): (firsts, counts)}
        self._index_cache = {}
        # {source: set(keys)}
        self._keys_cache = {}

    def __hash__(self):
        return hash(self.filename)

    def __eq__(self, other):
        return isinstance(other, FileAccess) and (other.filename == self.filename)

    def get_index(self, source, group):
        """Get first index & count for a source and for a specific train ID.

        Indices are cached; this appears to provide some performance benefit.
        """
        try:
            return self._index_cache[(source, group)]
        except KeyError:
            ix = self._read_index(source, group)
            self._index_cache[(source, group)] = ix
            return ix

    def _read_index(self, source, group):
        """Get first index & count for a source.

        This is 'real' reading when the requested index is not in the cache.
        """
        ix_group = self.file['/INDEX/{}/{}'.format(source, group)]
        firsts = ix_group['first'][:]
        if 'count' in ix_group:
            counts = ix_group['count'][:]
        else:
            status = ix_group['status'][:]
            counts = np.uint64((ix_group['last'][:] - firsts + 1) * status)
        return firsts, counts

    def get_keys(self, source):
        """Get keys for a given source name

        Keys are found by walking the HDF5 file, and cached for reuse.
        """
        try:
            return self._keys_cache[source]
        except KeyError:
            pass

        if source in self.control_sources:
            group = '/CONTROL/' + source
        elif source in self.instrument_sources:
            group = '/INSTRUMENT/' + source
        else:
            raise SourceNameError(source)

        res = set()

        def add_key(key, value):
            if isinstance(value, h5py.Dataset):
                res.add(key.replace('/', '.'))

        self.file[group].visititems(add_key)
        self._keys_cache[source] = res
        return res

class DataCollection:
    """An assemblage of data generated at European XFEL

    Data consists of *sources* which each have *keys*. It is further
    organised by *trains*, which are identified by train IDs.

    You normally get an instance of this class by calling :func:`H5File`
    for a single file or :func:`RunDirectory` for a directory.
    """
    def __init__(self, files, selection=None, train_ids=None):
        self.files = list(files)

        # selection: {source: set(keys)}
        # None as value -> all keys for this source
        if selection is None:
            selection = {}
            for f in self.files:
                selection.update(dict.fromkeys(f.control_sources))
                selection.update(dict.fromkeys(f.instrument_sources))
        self.selection = selection

        self.control_sources = set()
        self.instrument_sources = set()
        self._source_index = defaultdict(list)
        for f in self.files:
            self.control_sources.update(f.control_sources.intersection(selection))
            self.instrument_sources.update(f.instrument_sources.intersection(selection))
            for source in (f.control_sources | f.instrument_sources):
                self._source_index[source].append(f)

        if train_ids is None:
            train_ids = sorted(set().union(*(f.train_ids for f in files)))
        self.train_ids = train_ids

    @classmethod
    def from_paths(cls, paths):
        files = [FileAccess(h5py.File(path)) for path in paths]
        return cls(files)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass  # TODO: Should this close files?

    @property
    def all_sources(self):
        return self.control_sources | self.instrument_sources

    @property
    def detector_sources(self):
        return set(filter(DETECTOR_SOURCE_RE.match, self.instrument_sources))

    def _check_field(self, source, key):
        if source not in self.all_sources:
            raise SourceNameError(source)
        if key not in self._keys_for_source(source):
            raise PropertyNameError(key, source)

    def _keys_for_source(self, source):
        selected_keys = self.selection[source]
        if selected_keys is not None:
            return selected_keys

        # The same source may be in multiple files, but this assumes it has
        # the same keys in all files that it appears in.
        for f in self._source_index[source]:
            return f.get_keys(source)

    def _check_data_missing(self, tid) -> bool:
        """Return True if a train does not have data for all sources"""
        for source in self.control_sources:
            file, _ = self._find_data(source, tid)
            if file is None:
                return True

        for source in self.instrument_sources:
            file, pos = self._find_data(source, tid)
            if file is None:
                return True

            groups = {k.partition('.')[0] for k in self._keys_for_source(source)}
            for group in groups:
                _, counts = file.get_index(source, group)
                if counts[pos] == 0:
                    return True

        return False

    def trains(self, devices=None, train_range=None, *, require_all=False):
        """Iterate over all trains in the run and gather all sources.

        ::

            run = Run('/path/to/my/run/r0123')
            for train_id, data in run.trains():
                value = data['device']['parameter']

        Parameters
        ----------

        devices: dict or list, optional
            Filter data by devices and by parameters.

            Refer to :meth:`H5File.trains` for how to use this.

        train_range: by_id or by_index object, optional
            Iterate over only selected trains, by train ID or by index::

                f.trains(train_range=by_index[20:])

        require_all: bool
            False (default) returns any data available for the requested trains.
            True skips trains which don't have all the requested data;
            this requires that you specify required data using *devices*.

        Yields
        ------

        tid : int
            The train ID of the returned train
        data : dict
            The data for this train, keyed by device name
        """
        dc = self
        if devices is not None:
            dc = dc.select(devices)
        if train_range is not None:
            dc = dc.select_trains(train_range)
        return iter(TrainIterator(dc, require_all=require_all))

    def train_from_id(self, train_id, devices=None):
        """Get Train data for specified train ID.

        Parameters
        ----------

        train_id: int
            The train ID
        devices: dict or list, optional
            Filter data by devices and by parameters.

            Refer to :meth:`H5File.trains` for how to use this.

        Returns
        -------

        tid : int
            The train ID of the returned train
        data : dict
            The data for this train, keyed by device name

        Raises
        ------
        KeyError
            if `train_id` is not found in the run.
        """
        if devices is not None:
            return self.select(devices).train_from_id(train_id)

        res = {}
        for source in self.control_sources:
            source_data = res[source] = {}
            file, pos = self._find_data(source, train_id)
            if file is None:
                continue

            for key in self._keys_for_source(source):
                path = '/CONTROL/{}/{}'.format(source, key.replace('.', '/'))
                source_data[key] = file.file[path][pos]

        for source in self.instrument_sources:
            source_data = res[source] = {}
            file, pos = self._find_data(source, train_id)
            if file is None:
                continue

            for key in self._keys_for_source(source):
                group = key.partition('.')[0]
                firsts, counts = file.get_index(source, group)
                first, count = firsts[pos], counts[pos]
                if not count:
                    continue

                path = '/INSTRUMENT/{}/{}'.format(source, key.replace('.', '/'))
                if count == 1:
                    source_data[key] = file.file[path][first]
                else:
                    source_data[key] = file.file[path][first:first + count]

        return train_id, res

    def train_from_index(self, train_index, devices=None):
        """Get train data of the nth train in file.

        Parameters
        ----------
        index: int
            Index of the train in the file.
        devices: dict or list, optional
            Filter data by devices and by parameters.

            Refer to :meth:`~.H5File.trains` for how to use this.

        Returns
        -------

        tid : int
            The train ID of the returned train
        data : dict
            The data for this train, keyed by device name
        """
        train_id = self.train_ids[train_index]
        return self.train_from_id(train_id, devices=devices)

    def get_series(self, source, key):
        """Return a pandas Series for a particular data field.

        Parameters
        ----------

        source: str
            Device name with optional output channel, e.g.
            "SA1_XTD2_XGM/DOOCS/MAIN" or "SPB_DET_AGIPD1M-1/DET/7CH0:xtdf"
        key: str
            Key of parameter within that device, e.g. "beamPosition.iyPos.value"
            or "header.linkId". The data must be 1D in the file.
        """
        self._check_field(source, key)
        name = source + '/' + key
        if name.endswith('.value'):
            name = name[:-6]

        seq_series = []

        if source in self.control_sources:
            data_path = "/CONTROL/{}/{}".format(source, key.replace('.', '/'))
            for f in self._source_index[source]:
                data = f.file[data_path][:len(f.train_ids), ...]
                index = pd.Index(f.train_ids, name='trainId')

                seq_series.append(pd.Series(data, name=name, index=index))

        elif source in self.instrument_sources:
            data_path = "/INSTRUMENT/{}/{}".format(source, key.replace('.', '/'))
            for f in self._source_index[source]:
                group = key.partition('.')[0]
                firsts, counts = f.get_index(source, group)
                trainids = self._expand_trainids(counts, f.train_ids)

                index = pd.Index(trainids, name='trainId')
                data = f.file[data_path][:]
                if not index.is_unique:
                    pulse_id = f.file['/INSTRUMENT/{}/{}/pulseId'
                        .format(source, group)]
                    pulse_id = pulse_id[:len(index), 0]
                    index = pd.MultiIndex.from_arrays([trainids, pulse_id],
                                                      names=['trainId', 'pulseId'])
                    # Does pulse-oriented data always have an extra dimension?
                    assert data.shape[1] == 1
                    data = data[:, 0]
                data = data[:len(index)]

                seq_series.append(pd.Series(data, name=name, index=index))
        else:
            raise Exception("Unknown source category")

        return pd.concat(sorted(seq_series, key=lambda s: s.index[0]))

    def get_dataframe(self, fields=None, *, timestamps=False):
        """Return a pandas dataframe for given data fields.

        Parameters
        ----------
        fields : list of 2-tuples
            Glob patterns to match device and field names, e.g.
            ``("*_XGM/*", "*.i[xy]Pos")`` matches ixPos and iyPos from any XGM devices.
            By default, all fields from all control devices are matched.
        timestamps : bool
            If false (the default), exclude the timestamps associated with each
            control data field.
        """
        if fields is not None:
            return self.select(fields).get_dataframe(timestamps=timestamps)

        series = []
        for source in self.all_sources:
            for key in self._keys_for_source(source):
                if (not timestamps) and key.endswith('.timestamp'):
                    continue
                series.append(self.get_series(source, key))

        return pd.concat(series, axis=1)

    def get_array(self, source, key, extra_dims=None):
        """Return a labelled array for a particular data field.

        The first axis of the returned data will be the train IDs.
        Datasets which are per-pulse in the first dimension are not supported.

        Parameters
        ----------

        source: str
            Device name with optional output channel, e.g.
            "SA1_XTD2_XGM/DOOCS/MAIN" or "SPB_DET_AGIPD1M-1/DET/7CH0:xtdf"
        key: str
            Key of parameter within that device, e.g. "beamPosition.iyPos.value"
            or "header.linkId". The data must be 1D in the file.
        extra_dims: list of str
            Name extra dimensions in the array. The first dimension is
            automatically called 'train'. The default for extra dimensions
            is dim_0, dim_1, ...
        """
        self._check_field(source, key)
        seq_arrays = []

        if source in self.control_sources:
            data_path = "/CONTROL/{}/{}".format(source, key.replace('.', '/'))
            for f in self._source_index[source]:
                data = f.file[data_path][:len(f.train_ids), ...]
                if extra_dims is None:
                    extra_dims = ['dim_%d' % i for i in range(data.ndim - 1)]
                dims = ['trainId'] + extra_dims

                seq_arrays.append(xarray.DataArray(data, dims=dims,
                                                   coords={'trainId': f.train_ids}))

        elif source in self.instrument_sources:
            data_path = "/INSTRUMENT/{}/{}".format(source, key.replace('.', '/'))
            for f in self._source_index[source]:
                group = key.partition('.')[0]
                firsts, counts = f.get_index(source, group)
                if (counts > 1).any():
                    raise ValueError("{}/{} data has more than one data point per train"
                                     .format(source, group))
                trainids = self._expand_trainids(counts, f.train_ids)

                data = f.file[data_path][:len(trainids), ...]

                if extra_dims is None:
                    extra_dims = ['dim_%d' % i for i in range(data.ndim - 1)]
                dims = ['trainId'] + extra_dims

                seq_arrays.append(
                    xarray.DataArray(data, dims=dims, coords={'trainId': trainids}))
        else:
            raise SourceNameError(source)

        non_empty = [a for a in seq_arrays if (a.size > 0)]
        if not non_empty:
            if seq_arrays:
                # All per-file arrays are empty, so just return the first one.
                return seq_arrays[0]

            raise Exception(("Unable to get data for source {!r}, key {!r}. "
                             "Please report an issue so we can investigate")
                            .format(source, key))

        return xarray.concat(sorted(non_empty,
                                    key=lambda a: a.coords['trainId'][0]),
                             dim='trainId')

    def union(self, *others):
        files = set(self.files)
        train_ids = set(self.train_ids)

        for other in others:
            files.update(other.files)
            train_ids.update(other.train_ids)

        train_ids = sorted(train_ids)
        selection = union_selections([self.selection] +
                                     [o.selection for o in others])

        return DataCollection(files, selection=selection, train_ids=train_ids)

    def _expand_selection(self, selection):
        res = defaultdict(set)
        if isinstance(selection, dict):
            # {source: {key1, key2}}
            # {source: {}} -> all keys for this source
            for source, keys in selection.items():  #
                if source not in self.all_sources:
                    raise SourceNameError(source)

                res[source].update(keys or None)

        elif isinstance(selection, list):
            # selection = [('src_glob', 'key_glob'), ...]
            res = union_selections(self._select_glob(src_glob, key_glob)
                                   for (src_glob, key_glob) in selection)
        else:
            TypeError("Unknown selection type: {}".format(type(selection)))

        return dict(res)

    def _select_glob(self, source_glob, key_glob):
        source_re = re.compile(fnmatch.translate(source_glob))
        key_re = re.compile(fnmatch.translate(key_glob))
        if key_glob.endswith(('.value', '*')):
            ctrl_key_re = key_re
        else:
            # The translated pattern ends with "\Z" - insert before this
            p = key_re.pattern
            end_ix = p.rindex('\Z')
            ctrl_key_re = re.compile(p[:end_ix] + r'(\.value)?' + p[end_ix:])

        matched = {}
        for source in self.all_sources:
            if not source_re.match(source):
                continue

            if key_glob == '*':
                matched[source] = None
            else:
                r = ctrl_key_re if source in self.control_sources else key_re
                keys = set(filter(r.match, self._keys_for_source(source)))
                if keys:
                    matched[source] = keys

        if not matched:
            raise ValueError("No matches for pattern {}"
                             .format((source_glob, key_glob)))
        return matched

    def select(self, seln_or_source_glob, key_glob='*'):
        """Return a new DataCollection with selected sources & keys
        """
        if isinstance(seln_or_source_glob, str):
            seln_or_source_glob = [(seln_or_source_glob, key_glob)]
        selection = self._expand_selection(seln_or_source_glob)

        return DataCollection(self.files, selection=selection, train_ids=self.train_ids)

    def select_trains(self, train_range):
        if isinstance(train_range, by_id):
            start_ix = _tid_to_slice_ix(train_range.value.start, self.train_ids, stop=False)
            stop_ix = _tid_to_slice_ix(train_range.value.stop, self.train_ids, stop=True)
            ix_slice = slice(start_ix, stop_ix, train_range.value.step)
        elif isinstance(train_range, by_index):
            ix_slice = train_range.value
        else:
            raise TypeError(type(train_range))

        new_train_ids = self.train_ids[ix_slice]
        files = [f for f in self.files
                 if np.intersect1d(f.train_ids, new_train_ids).size > 0]

        return DataCollection(files, selection=self.selection, train_ids=new_train_ids)

    def _expand_trainids(self, counts, trainIds):
        n = min(len(counts), len(trainIds))
        return np.repeat(trainIds[:n], counts.astype(np.intp)[:n])

    def _find_data(self, source, train_id) -> (FileAccess, int):
        for f in self._source_index[source]:
            ixs = (f.train_ids == train_id).nonzero()[0]
            if ixs.size > 0:
                return f, ixs[0]

        return None, None

    def info(self):
        """Show information about the run.
        """
        # time info
        first_train = self.train_ids[0]
        last_train = self.train_ids[-1]
        train_count = len(self.train_ids)
        span_sec = (last_train - first_train) / 10
        span_txt = str(datetime.timedelta(seconds=span_sec))

        detector_srcs, non_detector_inst_srcs = [], []
        detector_modules = {}
        for source in self.instrument_sources:
            m = DETECTOR_SOURCE_RE.match(source)
            if m:
                detector_srcs.append(source)
                name, modno = m.groups((1, 2))
                detector_modules[(name, modno)] = source
            else:
                non_detector_inst_srcs.append(source)

        # A run should only have one detector, but if that changes, don't hide it
        detector_name = ','.join(sorted(set(k[0] for k in detector_modules)))

        # disp
        print('# of trains:   ', train_count)
        print('Duration:      ', span_txt)
        print('First train ID:', first_train)
        print('Last train ID: ', last_train)
        print()

        print("{} detector modules ({})".format(
            len(detector_srcs), detector_name
        ))
        if len(detector_modules) > 0:
            # Show detail on the first module (the others should be similar)
            mod_key = sorted(detector_modules)[0]
            mod_source = detector_modules[mod_key]
            dinfo = self.detector_info(mod_source)
            module = ' '.join(mod_key)
            dims = ' x '.join(str(d) for d in dinfo['dims'])
            print("  e.g. module {} : {} pixels".format(module, dims))
            print("  {} frames per train, {} total frames".format(
                dinfo['frames_per_train'], dinfo['total_frames']
            ))
        print()

        print(len(non_detector_inst_srcs), 'instrument sources (excluding detectors):')
        for d in sorted(non_detector_inst_srcs):
            print('  -', d)
        print()
        print(len(self.control_sources), 'control sources:')
        for d in sorted(self.control_sources):
            print('  -', d)
        print()

    def detector_info(self, source):
        """Get statistics about the detector data.

        Returns a dictionary with keys:
        - 'dims' (pixel dimensions)
        - 'frames_per_train'
        - 'total_frames'
        """
        all_counts = []
        for file in self._source_index[source]:
            _, counts = file.get_index(source, 'image')
            all_counts.append(counts)

        all_counts = np.concatenate(all_counts)
        dims = file.file['/INSTRUMENT/{}/image/data'.format(source)].shape[-2:]

        return {
            'dims': dims,
            # Some trains have 0 frames; max is the interesting value
            'frames_per_train': all_counts.max(),
            'total_frames': all_counts.sum(),
        }

    def train_info(self, train_id):
        """Show information about a specific train in the run.

        Parameters
        ----------
        train_id: int
            The specific train ID you get details information.

        Raises
        ------
        ValueError
            if `train_id` is not found in the run.
        """
        tid, files = next((t for t in self.ordered_trains
                          if t[0] == train_id), (None, None))
        if tid is None:
            raise ValueError("train {} not found in run.".format(train_id))
        ctrl, inst = self._get_sources(files)

        # disp
        print('Train [{}] information'.format(train_id))
        print('Devices')
        print('\tInstruments')
        [print('\t-', d) for d in sorted(inst)] or print('\t-')
        print('\tControls')
        [print('\t-', d) for d in sorted(ctrl)] or print('\t-')

class TrainIterator:
    def __init__(self, data, require_all=True):
        self.data = data
        self.require_all = require_all
        # {(source, key): (f, dataset)}
        self._datasets_cache = {}

    def _find_data(self, source, key, tid):
        try:
            file, ds = self._datasets_cache[(source, key)]
        except KeyError:
            pass
        else:
            ixs = (file.train_ids == tid).nonzero()[0]
            if ixs.size > 0:
                return file, ixs[0], ds

        data = self.data
        section = 'CONTROL' if source in data.control_sources else 'INSTRUMENT'
        path = '/{}/{}/{}'.format(section, source, key.replace('.', '/'))
        f, pos = data._find_data(source, tid)
        if f is not None:
            ds = f.file[path]
            self._datasets_cache[(source, key)] = (f, ds)
            return f, pos, ds

        return None, None, None

    def _assemble_data(self, tid):
        res = {}
        for source in self.data.control_sources:
            source_data = res[source] = {}
            for key in self.data._keys_for_source(source):
                _, pos, ds = self._find_data(source, key, tid)
                if ds is None:
                    continue
                source_data[key] = ds[pos]

        for source in self.data.instrument_sources:
            source_data = res[source] = {}
            for key in self.data._keys_for_source(source):
                file, pos, ds = self._find_data(source, key, tid)
                if ds is None:
                    continue
                group = key.partition('.')[0]
                firsts, counts = file.get_index(source, group)
                first, count = firsts[pos], counts[pos]
                if count == 1:
                    source_data[key] = ds[first]
                else:
                    source_data[key] = ds[first:first + count]

        return res

    def __iter__(self):
        for tid in self.data.train_ids:
            if self.require_all and self.data._check_data_missing(tid):
                continue
            yield tid, self._assemble_data(tid)

def H5File(path):
    """Open a single HDF5 file generated at European XFEL.

    Returns a :class:`DataCollection` object.

    Parameters
    ----------
    path: str
        Path to the HDF5 file
    """
    return DataCollection.from_paths([path])

def RunDirectory(path):
    """Open data files from a 'run' at European XFEL.

    A 'run' is a directory containing a number of HDF5 files with data from the
    same time period.

    Returns a :class:`DataCollection` object.

    Parameters
    ----------
    path: str
        Path to the run directory containing HDF5 files.
    """
    files = list(filter(h5py.is_hdf5, glob(osp.join(path, '*.h5'))))
    if not files:
        raise Exception("No HDF5 files found in {}".format(path))
    return DataCollection.from_paths(files)

def union_selections(selections):
    """Merge together different selections

    A selection is a dict of {source: set(keys)}, or {source: None}
    to include all keys for a given source.
    """
    selection_multi = defaultdict(list)

    for seln in selections:
        for source, keys in seln.items():
            selection_multi[source].append(keys)

    # Merge selected keys; None -> all keys selected
    return {source: None if (None in keygroups) else set().union(*keygroups)
            for (source, keygroups) in selection_multi.items()}

# RunDirectory was previously RunHandler; we'll leave it accessible in case
# any code was already using it.
RunHandler = RunDirectory


def stack_data(train, data, axis=-3, xcept=()):
    """Stack data from devices in a train.

    For detector data, use stack_detector_data instead: it can handle missing
    modules, which this function cannot.

    The returned array will have an extra dimension. The data will be ordered
    according to any groups of digits in the source name, interpreted as
    integers. Other characters do not affect sorting. So:

        "B_7_0" < "A_12_0" < "A_12_1"

    Parameters
    ----------
    train: dict
        Train data.
    data: str
        The path to the device parameter of the data you want to stack.
    axis: int, optional
        Array axis on which you wish to stack.
    xcept: list
        List of devices to ignore (useful if you have reccored slow data with
        detector data in the same run).

    Returns
    -------
    combined: numpy.array
        Stacked data for requested data path.
    """
    devices = [dev for dev in train.keys() if dev not in xcept]

    if not devices:
        raise ValueError("No data after filtering by 'xcept' argument.")

    dtypes, shapes = set(), set()
    ordered_arrays = []
    for device in sorted(devices, key=lambda d: list(map(int, re.findall(r'\d+', d)))):
        array = train[device][data]
        dtypes.add(array.dtype)
        ordered_arrays.append(array)

    if len(dtypes) > 1:
        raise ValueError("Arrays have mismatched dtypes: {}".format(dtypes))

    return np.stack(ordered_arrays, axis=axis)


def stack_detector_data(train, data, axis=-3, modules=16, only='', xcept=()):
    """Stack data from detector modules in a train.

    Parameters
    ----------
    train: dict
        Train data.
    data: str
        The path to the device parameter of the data you want to stack, e.g. 'image.data'.
    axis: int
        Array axis on which you wish to stack (default is -3).
    modules: int
        Number of modules composing a detector (default is 16).
    only: str
        Only use devices in train containing this substring.
    xcept: list
        List of devices to ignore (useful if you have reccored slow data with
        detector data in the same run).

    Returns
    -------
    combined: numpy.array
        Stacked data for requested data path.
    """
    devices = [dev for dev in train.keys() if only in dev and dev not in xcept]

    if not devices:
        raise ValueError("No data after filtering by 'only' and 'xcept' arguments.")

    dtypes, shapes = set(), set()
    modno_arrays = {}
    for device in devices:
        modno = int(re.search(r'/DET/(\d+)CH', device).group(1))
        array = train[device][data]
        dtypes.add(array.dtype)
        shapes.add(array.shape)
        modno_arrays[modno] = array

    if len(dtypes) > 1:
        raise ValueError("Arrays have mismatched dtypes: {}".format(dtypes))
    if len(shapes) > 1:
        raise ValueError("Arrays have mismatched shapes: {}".format(shapes))
    if max(modno_arrays) >= modules:
        raise IndexError("Module {} is out of range for a detector with {} modules"
                         .format(max(modno_arrays), modules))

    dtype = dtypes.pop()
    shape = shapes.pop()

    combined = np.full((modules, ) + shape, np.nan, dtype=dtype)
    for modno, array in modno_arrays.items():
        combined[modno] = array

    return np.moveaxis(combined, 0, axis)


if __name__ == '__main__':
    r = RunDirectory('./data/r0185')
    for tid, d in r.trains():
        print(tid)
