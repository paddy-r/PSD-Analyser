''' HR 01/12/22 To test setup file '''

from setuptools import setup, find_packages
import sys
sys.path.append(".")
setup(name = 'psdanalyser',
      version = '1.3.0', # Random, but current release executable is 1.2
      packages = find_packages(include = ['psdanalyser','psdanalyser.*']),
      install_requires = [])
