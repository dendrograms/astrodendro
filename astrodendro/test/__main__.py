"""
Find and run all unit tests in this package.
This allows you to run all unit tests with the command
python -m astrodendro.test
"""
from os.path import dirname
import unittest

suite = unittest.TestLoader().discover(dirname(__file__), 'test_*.py')
print("Running entire test suite...")
unittest.TextTestRunner().run(suite)