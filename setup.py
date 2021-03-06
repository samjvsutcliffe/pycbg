from setuptools import setup, find_packages

setup(
    name='pycbg',
    version='1.0.0',
    description='Python interface to generate CB-Geo mpm input files',
    author='Sacha Duverger',
    url='git_repos@axp20009:~/pycbg.git',
    install_requires=[
        'numpy',
        'gmsh',
        'pandas',
        'matplotlib',
        'sphinx>=3.3.1',
        'sphinx_rtd_theme'
    ],
    # extras_require={
    #     'documentation': ['sphinx>=3.3.1', 'sphinx_rtd_theme']
    # },
    packages=find_packages(include=['pycbg/*'])
)
