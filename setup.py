import setuptools
import versioneer

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pycbg',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Python scripts able to generate easily CB-Geo mpm input files',
    author='Sacha Duverger, Jérôme Duriez',
    author_email='sacha.duverger@inrae.fr, jerome.duriez@inrae.fr',
    url='https://forgemia.inra.fr/mpm-at-recover/pycbg.git',
    long_description=long_description,
    long_description_content_type='text/markdown',
    licence='MIT',
    install_requires=[
        'numpy',
        'gmsh',
        'pandas',
        'matplotlib',
        'versioneer',
        'sphinx>=3.3.1',
        'sphinx_rtd_theme'
    ],
    # extras_require={
    #     'documentation': ['sphinx>=3.3.1', 'sphinx_rtd_theme']
    # },
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src")
)
