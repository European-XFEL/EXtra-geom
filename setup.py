#!/usr/bin/env python
import os.path as osp
import re
from setuptools import setup, find_packages
import sys


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
          'cfelpyutils>=0.92',
          'h5py>=2.7.1',
          'matplotlib',
          'numpy',
          'scipy',
      ],
      extras_require={
          'docs': [
              'sphinx',
              'nbsphinx',
              'ipython',  # For nbsphinx syntax highlighting
              'sphinxcontrib_github_alt',
          ],
          'test': [
              'pytest',
              'pytest-cov',
              'coverage',
              'nbval',
              'testpath',
          ]
      },
      python_requires='>=3.6',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Physics',
      ]
)
