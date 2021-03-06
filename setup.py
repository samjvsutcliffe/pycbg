from setuptools import setup, find_packages

setup(
    name='pycbg',
    version='1.0.0',
    packages=find_packages(include=['preprocessing.py', 'postprocessing.py'])
)
