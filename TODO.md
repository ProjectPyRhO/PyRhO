PyRhO ToDo List
===============

Priority
--------
- [ ] Finish docstrings (rST) and generate more documentation with [Sphinx](http://sphinx-doc.org/tutorial.html)
- [ ] Add more tutorial notebooks
- [ ] Make more functions and variables private and use `__all__` to clear namespace

Mid-Term
--------
- [ ] Change flux units to mW/mm^2
- [ ] Add more asserts, warnings and write unit tests
- [ ] Standardise testing with tox
- [ ] Add Continuous Integration with Travis
- [ ] Move GUI to submodule
- [ ] Remove v_1 as a parameter (and add a dummy scaling factor for FV)
- [ ] Suffix all parameters to avoid clashes with Brian?

Long-Term
---------
- [ ] GUI for Brian
- [ ] Implement automatic unit handling and conversions
- [ ] Backport to Python 2?
- [ ] Incorporate additional independent variables
  - [ ] Temperature (Q10)
  - [ ] Wavelength
  - [ ] pH (intracellular and extracellular)
