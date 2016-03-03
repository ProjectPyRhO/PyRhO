PyRhO - A Virtual Optogenetics Laboratory
=========================================

A Python module to fit and characterise rhodopsin photocurrents

Optogenetics has become a key tool for understanding the function of neural circuits and controlling their behaviour. An array of directly light driven opsins have been genetically isolated from several families of organisms, with a wide range of temporal and spectral properties. In order to characterize, understand and apply these rhodopsins, we present an integrated suite of open-source, multi-scale computational tools called PyRhO. 

The purpose of developing PyRhO is threefold: 

(i) to characterize new (and existing) rhodopsins by automatically fitting a minimal set of experimental data to three, four or six-state kinetic models, 
(ii) to simulate these models at the channel, neuron & network levels and 
(iii) provide functional insights through model selection and virtual experiments *in silico*. 

The module is written in Python with an additional IPython/Jupyter notebook based GUI, allowing models to be fit, simulations to be run and results to be shared through simply interacting with a webpage. The seamless integration of model fitting algorithms with simulation environments for these virtual opsins will enable neuroscientists to gain a comprehensive understanding of their behaviour and rapidly identify the most suitable variant for application in a particular biological system. This process may thereby guide not only experimental design and opsin choice but also alterations of the rhodopsin genetic code in a neuro-engineering feed-back loop. In this way, we expect PyRhO will help to significantly improve optogenetics as a tool for transforming biological sciences. 

If you use PyRhO please cite our paper: 

Evans, B. D., Jarvis, S., Schultz, S. R. & Nikolic K. (2016) "PyRhO: A Multiscale Optogenetics Simulation Platform", *Front. Neuroinform., 10* (8). `doi:10.3389/fninf.2016.00008 <https://dx.doi.org/10.3389/fninf.2016.00008>`_

The PyRhO project website with additional documentation may be found here: `www.imperial.ac.uk/bio-modelling/pyrho <http://www.imperial.ac.uk/a-z-research/bio-modelling/pyrho>`_

Installation
------------

To install PyRhO (including the GUI) from PyPI use the command:
::
    pip install pyrho[full]
    
Alternatively, to install the latest code from Github (including the GUI) use the command:
::
    pip install git+https://github.com/ProjectPyRhO/PyRhO.git#egg=PyRhO[full]

Currently PyRhO only supports Python 3. To use PyRhO with the `NEURON simulator <http://www.neuron.yale.edu/neuron/>`_, NEURON must be compiled from its source code so that it works with Python 3. An installation script is provided for doing this on Mac OS X or Linux.  
The shell script may be called after importing PyRhO with the following function:
::
    from pyrho import *
    setupNEURON('/installation/path/for/NEURON/')
This will attempt to compile NEURON from source, copy the supplied ``mod`` and ``hoc`` files into place (the current working directory by default) finally compiling the ``mod`` files describing the opsin models ready for inclusion in simulations. 

The `Brian simulator <http://briansimulator.org/>`_ is included with the PyRhO installation for modelling networks of optogenetically transfected spiking neurons. 
