.. PyRhO documentation master file, created by
   sphinx-quickstart on Tue Mar 22 16:26:52 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyRhO's documentation!
=================================

Contents:

.. toctree::
   :maxdepth: 2

.. automodule:: pyrho


Abstraction Layers
==================

PyRhO contains three abstraction layers which can be used independently of one another. 

.. automodule:: pyrho.models
	:members:
.. automodule:: pyrho.protocols
	:members:
.. automodule:: pyrho.simulators
	:members:

In addition PyRhO has a fitting subpackage for optionally parameterising the rhodopsin 
models with experimental datasets. These must first be loaded into the appropriate data
structures as detailed below. 

.. automodule:: pyrho.expdata
    :members:
.. automodule:: pyrho.fitting
    :members:

Alternatively, parameters for the models (and protocols and simulators) can be loaded from
pre-defined sets. 

.. automodule:: pyrho.parameters
    :members:

To aid in the use of these abilities, several utility functions are also provided:

.. automodule:: pyrho.utilities
    :members:
.. automodule:: pyrho.config
    :members:

Finally there is an IPython/Jupyter notebook based graphical user interface (GUI) which 
can be used to fit and simulate the models. 


PyRhO API
=========

* fitModels(dataSet, nStates=3, params=modelParams['3'], method=methods[3])
    """Routine to fit as many models as possible and select between them according to some parsimony criterion"""
* plotFit(I, t, onInd, offInd, phi, V, nStates, params, fitRates, index)

* PhotoCurrent
    
* ProtocolData


Examples
========

Several common usage examples are given below which are available as a Jupyter notebook. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

