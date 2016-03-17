# Collection of experimental rhodopsin data sets

# http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
# foo_config = open(os.path.join(os.path.dirname(__file__),'foo.conf').read() # Old way

# from pkg_resources import resource_string # New way
# foo_config = resource_string(__name__, 'foo.conf')

from pkg_resources import resource_stream #, resource_string
import pickle

#from pyrho.dataSets import load_ChR2
# c.f. https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/datasets/base.py

#def load(resource='test.txt'):
#    return resource_string(__name__, resource)

def loadChR2():
    """Return dictionary of ChR2 protocol data.
        The dataset is a dictionary containing three ProtocolData objects:
        'step'      : Six PhotoCurrents recorded over flux values from
                      2.21e15 to 2.65e17 photons/mm2/s to steady-state
        'delta'     : The PhotoCurrent from the 'step' set with the strongest peak
        'shortPulse': Ten PhotoCurrents recorded with pulse durations of
                      [1, 2, 3, 4, 5, 6, 8, 10, 20, 30] ms
        With thanks to Nir Grossman, Juan Burrone and Matthew Grub.
    """
    return pickle.load(resource_stream(__name__, 'ChR2.pkl'))

'''
def loadpkl(pkl, path=None):
    # with, try, finaly etc...
    #import os
    pklFile = pkl+".pkl"
    if path is None:
        dataSet = resource_string(__name__, pklFile)

        cwd = os.getcwd()
        if pkl in cwd: ### Finish!!!
            pass #pklFile = pklFile
            #fh = open(pklFile, "rb")
        else:
            pklFile = os.path.join(config., pklFile)
            #fh = open(os.path.join(config., pklFile), "rb")
    else:
        pklFile = os.path.join(path, pklFile)
    with open(pklFile, "rb") as fh :
        dataSet = pickle.load(fh)
    #fh.close()
    return dataSet
'''