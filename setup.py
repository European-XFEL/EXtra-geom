#!/usr/bin/env python
import os.path as osp
import re
import sys

from setuptools import find_packages, setup


def get_script_path():
    return osp.dirname(osp.realpath(sys.argv[0]))


def read(*parts):
    return open(osp.join(get_script_path(), *parts)).read()


def find_version(*parts):
    vers_file = read(*parts)
    match = re.search(r'^__version__ = "(\d+\.\d+\.\d+)"', vers_file, re.M)
    if match is not None:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(name="EXtra-geom",
      version=find_version("extra_geom", "__init__.py"),
      author="European XFEL GmbH",
      author_email="da-support@xfel.eu",
      maintainer="Thomas Michelat",
      url="https://github.com/European-XFEL/EXtra-geom",
      description="Tools to work with EuXFEL detector geometry and assemble detector images",
      long_description=read("README.md"),
      long_description_content_type='text/markdown',
      license="BSD-3-Clause",
      packages=find_packages(),
      package_data={
          'extra_geom.tests': ['dssc_geo_june19.h5', 'lpd_mar_18.h5'],
      },
      install_requires=[
          'cfelpyutils>=0.92, <2.0',
          'h5py>=2.7.1',
          'matplotlib',
          'numpy',
      ],
      extras_require={
          'interpolate': ['scipy'],
          'docs': [
              'sphinx',
              'nbsphinx',
              'ipython',  # For nbsphinx syntax highlighting
              'sphinxcontrib_github_alt',
              'sphinx-rtd-theme',
              # Pin docutils until sphinx-rtd-theme issue 1115 is fixed
              'docutils==0.16',
          ],
          'test': [
              'pytest',
              'pytest-cov',
              'coverage',
              'nbval',
              'testpath',
              'xarray',
              'EXtra-data',
          ]
      },
      python_requires='>=3.6',
      classifiers=[
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      project_urls={
          'Documentation': 'https://extra-geom.readthedocs.io/en/latest/',
          'Changelog': 'https://extra-geom.readthedocs.io/en/latest/changelog.html',
      },
)
