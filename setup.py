#!/usr/bin/env python

from setuptools import setup, Command

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py


class DendroTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import sys
        import subprocess
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

setup(name='dendro-core',
      version='0.0.1',
      description='Python package for computation of astronomical dendrograms',
      author='Braden MacDonald and Thomas Robitaille',
      author_email='braden@bradenmacdonald.com',
      packages=['astrodendro', 'astrodendro.io', 'astrodendro.test'],
      provides=['astrodendro'],
      requires=['numpy'],
      cmdclass={'build_py': build_py, 'test': DendroTest},
      keywords=['Scientific/Engineering'],
      classifiers=[
                   "Development Status :: 4 - Beta",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: MIT License",
                  ],
     )
