"""The PyRhO package setup script"""

#from __future__ import print_function        # Added for Python 2.x support
from setuptools import setup, find_packages  # Prefer setuptools over distutils
from codecs import open                      # To use a consistent encoding
import os

# Download and install setuptools if not installed
#from ez_setup import use_setuptools
#use_setuptools()
#python -m ensurepip --upgrade

#from setuptools import setup
#from distutils import setup

here = os.path.abspath(os.path.dirname(__file__))
home = os.path.expanduser("~")
print(home)
prwd = os.path.join(home, 'pyrho')  # pyrho working directory

# TODO: Test changes to package_data and include notebooks and license without MANIFEST

# TODO: Fix this to remove redundant long_description text
# Get the long description from the relevant file
#with open(os.path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
#with open('DESCRIPTION.rst', encoding='utf-8') as f:
#    long_description = f.read()

long_description = """
PyRhO - A Virtual Optogenetics Laboratory
=========================================

A Python module to fit and characterise rhodopsin photocurrents.

Background
----------

Optogenetics has become a key tool for understanding the function of neural circuits and controlling their behaviour. An array of directly light driven opsins have been genetically isolated from several families of organisms, with a wide range of temporal and spectral properties. In order to characterize, understand and apply these rhodopsins, we present an integrated suite of open-source, multi-scale computational tools called PyRhO.

PyRhO enables users to:

(i) characterize new (and existing) rhodopsins by automatically fitting a minimal set of experimental data to three, four or six-state kinetic models,
(ii) simulate these models at the channel, neuron & network levels and
(iii) gain functional insights through model selection and virtual experiments *in silico*.

The module is written in Python with an additional IPython/Jupyter notebook based GUI, allowing models to be fit, simulations to be run and results to be shared through simply interacting with a webpage. The seamless integration of model fitting algorithms with simulation environments for these virtual opsins will enable (neuro)scientists to gain a comprehensive understanding of their behaviour and rapidly identify the most suitable variant for application in a particular biological system. This process may thereby guide not only experimental design and opsin choice but also alterations of the rhodopsin genetic code in a neuro-engineering feed-back loop. In this way, we hope PyRhO will help to significantly improve optogenetics as a tool for transforming biological sciences.

Further Information
-------------------

If you use PyRhO please cite our paper:

Evans, B. D., Jarvis, S., Schultz, S. R. & Nikolic K. (2016) "PyRhO: A Multiscale Optogenetics Simulation Platform", *Frontiers in Neuroinformatics, 10* (8). `doi:10.3389/fninf.2016.00008 <https://dx.doi.org/10.3389/fninf.2016.00008>`_

The PyRhO project website with additional documentation may be found here: `www.imperial.ac.uk/bio-modelling/pyrho <http://www.imperial.ac.uk/a-z-research/bio-modelling/pyrho>`_

Finally, don't forget to follow us on twitter for updates: `@ProjectPyRhO <https://twitter.com/ProjectPyRhO>`_!

"""


setup(
    name='PyRhO',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.9.5',

    description='Fit and characterise rhodopsin photocurrents',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/ProjectPyRhO/PyRhO/',
    # download_url='https://github.com/ProjectPyRhO/PyRhO/archive/master.zip',
    # download_url='https://github.com/ProjectPyRhO/PyRhO/tarball/' + version,

    # Author details
    author='Benjamin D. Evans',
    author_email='ben.d.evans@gmail.com',

    license='BSD',
    platforms=['Linux', 'Mac OS X', 'Windows'],
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Life',
        'Topic :: Scientific/Engineering :: Human Machine Interfaces',

        # The license should match "license" above
        'License :: OSI Approved :: BSD License',

        # Supported Python versions
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        # 3.6 EOL: 23/12/21
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Framework :: IPython',
        'Natural Language :: English',
        'Operating System :: OS Independent',
    ],

    #keywords='optogenetics rhodopsin opsin brain neuroscience neuron brian jupyter',
    keywords=['optogenetics', 'rhodopsin', 'opsin', 'brain', 'neuroscience',
              'neuron', 'brian', 'jupyter'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    # package_dir = {'':'.'},
    # package_dir = {'pyrho': 'pyrho'}, # Relative to this script

    # List run-time dependencies here.  These will be installed by pip when your
    # project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # ipython is used for latex repr - remove from requirements and have a fallback repr?
    install_requires=['numpy>=1.8', 'scipy>=0.15', 'matplotlib>=1.3',
                      'lmfit>=0.9.3,<1.0.3', 'brian2>=2.0'],  # 'ipython>=4.1'

    # List additional groups of dependencies here (e.g. development dependencies).
    # You can install these using the following syntax, for example:
    # $ pip install -e .[dev,test]
    extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #    'brian' : ['brian2'],
    #    'docs' : ['sphinx>=1.3'],
        'extras': ['seaborn>=0.7', 'pandas>=0.17'],  # 'cython>=0.23'
    # traitlets is a dependency of ipywidgets and can be removed if 4.1 entails traitlets>=4.1
        'GUI' : ['jupyter>=1.0', 'notebook>=4.1', 'ipywidgets>=4.1,<5',
                 'seaborn>=0.7'],  # , 'traitlets>=4.1,<5'
        'full': ['jupyter>=1.0', 'notebook>=4.1', 'ipywidgets>=4.1,<5',
                 'seaborn>=0.7', 'pandas>=0.17'],  # 'cython>=0.23'
    },

    include_package_data=True,
    # If there are data files included in your packages that need to be
    # installed, specify them here. If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        # TODO: Try this without MANIFEST
        'NEURON'  : ['*.mod', '*.hoc', '*.sh'],
        'gui'     : ['*.png'],
        'datasets': ['*.pkl'],
    },


    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages.
    # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],
    #data_files=[#(prwd, []),
    #            (prwd, [os.path.join(prwd, 'gui/*.png'), ])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)
