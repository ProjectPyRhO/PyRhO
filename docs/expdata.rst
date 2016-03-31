Experimental Data Structures
============================

Each voltage-clamp recorded photocurrent should be loaded into a ``PhotoCurrent`` object as follows: 
.. code-block::

    pc = PhotoCurrent(I=i0, t=t, pulses=[[t\_on, t\_off]], phi=2e15, V=-70)
    
Here, ``I`` is the array of photocurrent values in nanoamperes, ``t`` is the corresponding array of recording times (or a scalar representing the time-step) in milliseconds, ``pulses`` is a nested list of :math:`n` lists (or an :math:`n \times 2` array), where :math:`n` corresponds to the number of pulses and each inner list contains the on-time and off-time for each pulse in milliseconds, ``phi`` represents the stimulating flux value in :math:`\mathrm{photons\cdot mm^{-2}\cdot s^{-1}}` and ``V`` is the clamp voltage in millivolts (or ``None`` if the voltage was not clamped). 

The ``PhotoCurrent`` class contains methods which automatically check the data and extract the key features from it, which may then be accessed as properties of the object with the ``.`` operator. Any properties which are data-derived are suffixed with ``_``, for example, the peak and steady-state current are accessed with ``pc.peak_`` and ``pc.ss_`` respectively. These photocurrents may easily be plotted, along with their main features and the light stimulus using the ``plot()`` method. 

The ``PhotoCurrent`` objects are then combined into ``ProtocolData`` objects. These sets of photocurrents also provide several convenient methods for plotting and extracting parameters from the data set as a whole. 

.. code-block::

    stepPD = ProtocolData(protocol="step", nRuns=1, 
                    phis=[1e14,1e15,1e16,1e17,1e18], Vs=[-70])

    for iPhi, phi in enumerate(phis):
        for iV, V in enumerate(Vs):
            pc = PhotoCurrent(Is[iPhi][iV], t, pulses, phi, V)
            stepPD.addTrial(pc)


Finally, the data sets are combined into a dictionary using the protocol names as keys: 

.. code-block::

    ChR2DataSet = {"step" : stepPD, 
                "recovery" : recovPD, 
                "rectifier" : rectiPD, 
                "shortPulse" : shortPD} 


This dictionary contains all the data necessary to parameterise all three models, however, if only the three and four-state models are of interest then the ``"shortPulse"`` protocol may be omitted. 

The data are now in a form suitable to be passed to the fitting algorithms as demonstrated in :doc:`fitting`. 

Package: ``expdata`` API
========================

.. automodule:: pyrho.expdata
    :members: