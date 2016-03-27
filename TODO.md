PyRhO TODO List
===============

Priority
--------
- [x] Make more functions and variables private and use `__all__` and `del` to tidy namespace
- [ ] Documentation
  - [ ] Finish docstrings (rST)
  - [ ] Generate more documentation with [Sphinx](http://sphinx-doc.org/tutorial.html)
  - [ ] Add to http://pythonhosted.org/PyRhO and/or https://readthedocs.org/
- [ ] Testing
  - [ ] Add more asserts, warnings and write unit tests with [`py.test`](http://pytest.org/latest/) or [`nose`](http://nose.readthedocs.org/en/latest/)
  - [ ] Standardise testing different versions with [`tox`](http://tox.readthedocs.org/en/latest/)
  - [ ] Add Continuous Integration with [Travis](https://travis-ci.org/)

Mid-Term
--------
- [ ] Add more tutorial notebooks
- [ ] Extend database of `modelFits` and add `general` with a wide range
- [ ] Change flux units to mW/mm^2
- [ ] Move GUI to submodule and add `modelFits` dropdowns (fitting and simulating)
- [ ] Remove v_1 as a parameter, add procedure to mod files (and add a dummy scaling factor for FV)
- [ ] Suffix all parameters to avoid clashes with Brian?

Long-Term
---------
- [ ] GUI for Brian
- [ ] Implement automatic unit handling and conversions
- [x] Backport to Python 2: 2.7
- [ ] Incorporate additional independent variables
  - [ ] Temperature (Q10)
  - [ ] Wavelength
  - [ ] pH (intracellular and extracellular)
