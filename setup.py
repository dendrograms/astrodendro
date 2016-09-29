#!/usr/bin/env python

from setuptools import setup, Command

from distutils.command.build_py import build_py


class DendroTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):

        import os
        import shutil
        import tempfile

        # First ensure that we build the package so that 2to3 gets executed
        self.reinitialize_command('build', inplace=False)
        self.run_command('build')
        build_cmd = self.get_finalized_command('build')
        new_path = os.path.abspath(build_cmd.build_lib)

        # Copy the build to a temporary directory for the purposes of testing
        # - this avoids creating pyc and __pycache__ directories inside the
        # build directory
        tmp_dir = tempfile.mkdtemp(prefix='astrodendro-test-')
        testing_path = os.path.join(tmp_dir, os.path.basename(new_path))
        shutil.copytree(new_path, testing_path)

        import sys
        import subprocess

        errno = subprocess.call([sys.executable, os.path.abspath('runtests.py')], cwd=testing_path)
        raise SystemExit(errno)


setup(name='astrodendro',
      version='0.2.0',
      url='http://www.dendrograms.org',
      description='Python package for computation of astronomical dendrograms',
      author='Thomas Robitaille, Chris Beaumont, Adam Ginsburg, Braden MacDonald, and Erik Rosolowsky',
      author_email='thomas.robitaille@gmail.com',
      packages=['astrodendro', 'astrodendro.io', 'astrodendro.tests'],
      package_data={'astrodendro.tests':['*.npz', 'benchmark_data/*fits']},
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
