import setuptools
from version import get_git_version

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pycbg',
    version=get_git_version(),
    description='Python scripts able to generate easily CB-Geo mpm input files',
    author='Sacha Duverger',
    author_email='sacha.duverger@inrae.fr',
    url='git_repos@axp20009:~/pycbg.git',
    long_description=long_description,
    long_description_content_type='text/markdown',
    licence='MIT',
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
    packages=setuptools.find_packages(include=['pycbg/*'])
)