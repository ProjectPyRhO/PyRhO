
import os.path
import os
import ast
import copy
from collections import OrderedDict

import numpy as np
import ipywidgets as widgets
from IPython.display import display
from IPython.display import clear_output
from lmfit import Parameters

from pyrho.parameters import (modelParams, modelList, modelLabels,
                              protList, protParams,
                              simList, simParams, simUnitLabels,
                              unitLabels, stateLabs,
                              protParamLabels, protUnitLabels, protParamNotes,
                              PyRhOparameters)
from pyrho.models import *
from pyrho.simulators import *
from pyrho.protocols import *
from pyrho.expdata import *
from pyrho.fitting import *
from pyrho.config import simAvailable, GUIdir, setupGUI  # , dDir
from pyrho.utilities import *
from pyrho import config
from pyrho.config import verbose

__all__ = ['loadGUI']

has_gui = {'Python': True, 'NEURON': True, 'Brian': False}
statesDict = OrderedDict([(s, i) for i, s in enumerate(list(modelParams))])  # .keys()
statesArray = modelList  # statesArray = list(statesDict)


boolDict = OrderedDict([('True', True), ('False', False)])
visibilityMap = ['hidden', '']  # [False, True]
displayMap = ['none', '']  # [False, True]

class ParamWidgets(object):
    """Common base class for all sets of parameter widgets"""

    def __init__(self, pSet, params_type=None):
        self.defaults = pSet
        self.type = params_type

        # Create lists and dictionaries of titles
        TabGroups = {'Fit': 0, 'Models': 1, 'Protocols': 2, 'Simulators': 3}
        #self.fit_methods = OrderedDict([(m, i) for i, m in enumerate(methods)])
        self.fit_methods = methods  # Use self.fit_methods.index('powell')

        self.titles_models = ['Three-state model', 'Four-state model', 'Six-state model']
        # Get lists from pyrho.parameters
        self.models = modelList
        # mParamsK2I = OrderedDict([(model,OrderedDict([(p,i) for i,p in enumerate(list(modelParams[model]))])) for model in modelList])
        # mParamsI2K = OrderedDict([(model,list(modelParams[model])) for model in modelList])
        self.params_models = OrderedDict([(model, list(modelParams[model])) for model in modelList])

        self.paramsK2I = OrderedDict([(s, OrderedDict([(p, i) for i, p in enumerate(list(modelParams[s]))]))
                                     for s in statesArray])
        self.paramsI2K = OrderedDict([(s, list(modelParams[s])) for s in statesArray])

        self.protocols = protList
        # pParamsK2I = OrderedDict([(prot,OrderedDict([(p,i) for i,p in enumerate(list(protParams[prot]))])) for prot in protList])
        # pParamsI2K = OrderedDict([(prot,list(protParams[prot])) for prot in protList])
        self.params_protocols = OrderedDict([(prot, list(protParams[prot])) for prot in protList])

        self.sims = [sim for sim in simList if simAvailable(sim) and has_gui(sim)]
        # sParamsK2I = OrderedDict([(sim,OrderedDict([(p,i) for i,p in enumerate(list(simParams[sim]))])) for sim in simList])
        # sParamsI2K = OrderedDict([(sim,list(simParams[sim])) for sim in simList])
        simParamsList = OrderedDict([(sim, list(simParams[sim])) for sim in simList])

        self.widgets_models = OrderedDict([([widgets], list(modelParams[model])) for model in modelList])




    def __str__(self):
        return "Parameter set: " + self.type

    def getParams(self):  # , pSet, valueList, varyList=None, minList=None, maxList=None, exprList=None
        '''Read values of parameters from the widgets'''
        userParams = Parameters()
        pSet = self.defaults
        for i, param in enumerate(pSet):
            if isinstance(pSet[param].value, list):  # TODO: Change to not number!!!
                userParams.add(param, value=ast.literal_eval(widgetList[i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            else:
                userParams.add(param, value=widgetList[i].value)
            if varyList is not None:
                userParams[param].set(vary=varyList[i].value) #userParams[param].vary = varyList[i].value
            if minList is not None:
                userParams[param].set(min=minList[i].value)
            if maxList is not None:
                userParams[param].set(max=maxList[i].value)
            if exprList is not None:
                userParams[param].set(expr=maxList[i].value)

        return userParams

    def setGUIparams(pSet, widgetList, varyList=None, minList=None, maxList=None, exprList=None):
        for i, param in enumerate(pSet):
            print(i, param)
            if isinstance(pSet[param].value, list):  # TODO: Change to not number!!!
                widgetList[i].value = str(pSet[param].value)
            else:
                widgetList[i].value = pSet[param].value
            if varyList is not None:
                varyList[i].value = pSet[param].vary
            if minList is not None:
                minList[i].value = pSet[param].min
            if maxList is not None:
                maxList[i].value = pSet[param].max
            if exprList is not None and pSet[param].expr is not None:
                exprList[i].value = pSet[param].expr  ### Change units handling

    def setParams(self, params):
        for p in params.keys():
            #if p in self.__dict__:
            #self.__dict__[p] = params[p].value
            setattr(self, p, params[p].value)
            #else:
            #    warnings.warn('Warning: "{p}" not found in {self}'.format(p,self))

    def exportParams(self, params):
        """Export parameters to lmfit dictionary"""
        for p in self.__dict__.keys():
            params[p].value = self.__dict__[p]
        return params

    ##### Not used yet... #####
    def buildLayerTab(paramGroup):

        pluginBoxes = [None for plugin in paramGroup] #{plugin: None for plugin in paramGroup}
        pluginParamsBox = [None for plugin in paramGroup] #{plugin: None for plugin in paramGroup} # Left side container
        pluginNotesBox = [None for plugin in paramGroup] #{plugin: None for plugin in paramGroup} # Right side container for figure and equations
        #pluginStimHTML = [None for plugin in paramGroup]
        pluginFigHTML = [None for plugin in paramGroup] #{plugin: None for plugin in paramGroup}
        #eqBox = [None for m in range(len(modelParams))]
        pluginValues = [[None for param in paramGroup[plugin]] for plugin in paramGroup] #{plugin:{param: None for param in paramGroup[plugin]} for plugin in paramGroup} # Array of parameter values
        pluginUnits = [[None for param in paramGroup[plugin]] for plugin in paramGroup] #{plugin:{param: None for param in paramGroup[plugin]} for plugin in paramGroup} # Array of units
        pluginParams = [[None for param in paramGroup[plugin]] for plugin in paramGroup] #{plugin:{param: None for param in paramGroup[plugin]} for plugin in paramGroup} # Array of parameter boxes

        for pluginInd, plugin in enumerate(paramGroup):
            #pluginInd = paramGroup.keys().index(plugin)
            pSet = paramGroup[plugin]
            for pInd, param in enumerate(pSet):
                #paramInd = pSet.keys().index(param)
                if isinstance(pSet[param].value, list): # list ==> Text
                    pluginValues[pluginInd][pInd] = widgets.Text(value=str(pSet[param].value), description=param) # np.asarray
                elif isinstance(pSet[param].value, str): # str ==> Text
                    pluginValues[pluginInd][pInd] = widgets.Text(value=str(pSet[param].value), description=param)
                elif isinstance(pSet[param].value, bool):
                    pluginValues[pluginInd][pInd] = widgets.Dropdown(options=boolDict,value=pSet[param].value,description=param)
                else:  # Numeric
                    if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf): # No limits
                        pluginValues[pluginInd][pInd] = widgets.FloatText(value=pSet[param].value, description=param)
                    else: # Bounded # ==> widgets.FloatSlider() ?
                        pluginValues[pluginInd][pInd] = widgets.BoundedFloatText(value=pSet[param].value, min=pSet[param].min, max=pSet[param].max, description=param)
                    pluginValues[pluginInd][pInd].width = '150px'
                if pSet[param].expr is None: # No units
                    pluginParams[pluginInd][pInd] = widgets.HBox(children=[pluginValues[pluginInd][pInd]])
                else:
                    pluginUnits[pluginInd][pInd] = widgets.Dropdown(options=[pSet[param].expr],value=pSet[param].expr) ### Change units handling
                    pluginParams[pluginInd][pInd] = widgets.HBox(children=[pluginValues[pluginInd][pInd],pluginUnits[pluginInd][pInd]])



            pluginFigHTML[pluginInd] = widgets.HTML()
            #exampleProt = '{}{}6s.{}'.format(fDir,prot,'png')#saveFigFormat)
            #if os.path.isfile(exampleProt):
            #    protFigHTML[pInd].value='<img src="{}" alt=Example {} width=200px>'.format(exampleProt,prot)
            #else:
            #    protFigHTML[pInd].value='Example Figure'
            pluginParamsBox[pluginInd] = widgets.Box(children=pluginParams[pluginInd])
            pluginNotesBox[pluginInd] = widgets.HBox(children=[pluginFigHTML[pluginInd]])# simStimHTML[sInd]  , ])#[figHTML[prot],eqBox[prot]])

            pluginBoxes[pluginInd] = widgets.HBox(children=[pluginParamsBox[pluginInd],pluginNotesBox[pluginInd]])#modelBox
            #display(protBoxes[pInd])
            pluginBoxes[pluginInd].margin = '5px'


        ##### Plugin parameters tab #####
        pluginParamsTabs = widgets.Tab(description='Plugin Settings', children=pluginBoxes)# \
        pluginParamsTabs.margin = '5px'

        return pluginParamsTabs
