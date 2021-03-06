from setuptools import setup, find_packages

setup(
    name='pycbg',
    version='1.0.0',
    description='Python interface to generate CB-Geo mpm input files',
    author='Sacha Duverger',
    packages=find_packages(include=['pycbg/*'])
)
