.. PyRhO documentation master file, created by
   sphinx-quickstart on Tue Mar 22 16:26:52 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyRhO documentation
===================

.. automodule:: pyrho
.. currentmodule:: pyrho

Introduction
============

.. include:: ../README.rst

Contents
--------

.. toctree::
   :maxdepth: 2
   
   models
   protocols
   simulators
   expdata
   fitting
   parameters
   utilities
   config


Abstraction Layers
==================

PyRhO contains three abstraction layers which can be used independently of one another. 

* :doc:`models`
    .. automodule:: pyrho.models
* :doc:`protocols`
    .. automodule:: pyrho.protocols
* :doc:`simulators`
    .. automodule:: pyrho.simulators
	
.. image:: figs/architecture.png
    
Model Fitting
=============

In addition PyRhO has a fitting subpackage for optionally parameterising the rhodopsin 
models with experimental datasets. These must first be loaded into the appropriate data
structures as detailed below. 

:doc:`expdata`
    .. automodule:: pyrho.expdata
    
:doc:`fitting`
    .. automodule:: pyrho.fitting
    

Alternatively, parameters for the models (and protocols and simulators) can be loaded from pre-defined sets. 

:doc:`parameters`
    .. automodule:: pyrho.parameters
    

To aid in the use of these abilities, several utility functions are also provided:

:doc:`utilities`
    .. automodule:: pyrho.utilities
    
:doc:`config`
    .. automodule:: pyrho.config
    

Finally there is an IPython/Jupyter notebook based graphical user interface (GUI) which can be used to fit and simulate the models. 



Examples
========

Several common usage examples are given below which are available as a Jupyter notebook. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

