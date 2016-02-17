#TODO: Refactor the whole GUI! 

#from IPython.html.widgets import interact
#from IPython.html import widgets # IPython < 4

from __future__ import print_function
import ipywidgets as widgets
#from IPython.utils.traitlets import link # IPython < 4
#from traitlets import link
from traitlets import Unicode
from IPython.display import display
from IPython.display import clear_output
#from IPython import get_ipython

#from IPython.core.pylabtools import print_figure
#import base64

from pyrho.models import *
from pyrho.simulators import *
from pyrho.protocols import *
from pyrho.loadData import *
from pyrho.fitting import *
from pyrho.config import GUIdir, setupGUI, dDir # For dataSet loading
from pyrho.utilities import *
import warnings
import time
from collections import OrderedDict
import os.path
import os
import ast
import copy


# LaTeX in widget descriptions/labels
# FloatRangeSlider and IntRangeSlider
# An Output widget was added, which allows you to print and display within widgets - replace Popup
# A SelectMultiple widget was added
# Suppress widget warning
# Add placeholder attribute to text widgets
# Tooltip on toggle button
# Dropdown options can be a dict, tuple or list


### Enable these for NEURON!!!
#import neuron
#from neuron import h

#%precision %.6g
#import numpy as np
#np.set_printoptions(precision=6)
#pprint()


# Create lists and dictionaries of titles
modelTitles = ['Three-state model', 'Four-state model', 'Six-state model']

fitMethodsDict = OrderedDict([(m,i) for i,m in enumerate(methods)])

# State keys must be padded with a leading ' ' to avoid a widgets bug: https://github.com/ipython/ipython/issues/6469
#statesDict = OrderedDict([(' '+s,i) for i,s in enumerate(list(modelParams.keys()))]) # enumerate(modelList)
#statesDict = OrderedDict([(' 3',0), (' 4',1), (' 6',2)]) 
statesDict = OrderedDict([(s,i) for i,s in enumerate(list(modelParams))]) #.keys()
statesArray = modelList #statesArray = list(statesDict) #.keys() #[' 3', ' 4', ' 6'] # [u' 3',u' 4',u' 6'] ### Redundant!

TabGroups = {'Fit':0, 'Models':1, 'Protocols':2, 'Simulators':3}
#TabGroups = {'Models':0, 'Simulators':1, 'Protocols':2}

#clearDelay = 1.5 # Pause [s] before clearing text entry fields

# Structures for cross-referencing arrays of widgets to their corresponding parameters
# http://stackoverflow.com/questions/18809482/python-nesting-dictionary-ordereddict-from-collections
#mParamsK2I = OrderedDict([ (model,OrderedDict([(p,i) for i,p in enumerate(list(modelParams[model]))])) for model in modelList ])
#mParamsI2K = OrderedDict([ (model,list(modelParams[model])) for model in modelList ]) 
modelParamsList = OrderedDict([ (model, list(modelParams[model])) for model in modelList ])

#sParamsK2I = OrderedDict([ (sim,OrderedDict([(p,i) for i,p in enumerate(list(simParams[sim]))])) for sim in simList ])
#sParamsI2K = OrderedDict([ (sim,list(simParams[sim])) for sim in simList ])
simParamsList = OrderedDict([ (sim, list(simParams[sim])) for sim in simList ]) 
loadedSims = [ sim for sim in simList if simAvailable(sim) ]

#pParamsK2I = OrderedDict([ (prot,OrderedDict([(p,i) for i,p in enumerate(list(protParams[prot]))])) for prot in protList ])
#pParamsI2K = OrderedDict([ (prot,list(protParams[prot])) for prot in protList ]) 
protParamsList = OrderedDict([ (prot, list(protParams[prot])) for prot in protList ])

boolDict = OrderedDict([('True',True), ('False',False)])



### TODO: Replace GUI with object oriented code...
class ParamWidgets(object):
    """Common base class for all sets of parameter widgets"""
    
    def __init__(self, pSet, type=None):
        self.defaults = pSet
        self.type = type
        self.paramsK2I = OrderedDict([ (set, OrderedDict([(p,i) for i,p in enumerate(list(modelParams[set]))])) for set in statesArray ])
        self.paramsI2K = OrderedDict([ (set, list(modelParams[s])) for set in statesArray ]) 
        
    def __str__(self):
        return "Parameter set: "+self.type
    
    def getParams(self): #, pSet, valueList, varyList=None, minList=None, maxList=None, exprList=None
        userParams = Parameters()
        for i, param in enumerate(pSet):
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
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
            print(i,param)
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
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
                exprList[i].value = pSet[param].expr ### Change units handling
    
    def setParams(self, params):
        for p in params.keys():
            #if p in self.__dict__:
            self.__dict__[p] = params[p].value
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
                else: # Numeric
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
    
    #pluginParamsTabs.on_trait_change(onChangeSimTab,'selected_index') # External
    


    
    

def loadGUI(IPythonWorkspace=None):
    
    path = os.path.join(os.getcwd(), GUIdir)
    if not os.path.isdir(path):
        setupGUI(path)
    
    if IPythonWorkspace is None:
        pass
    else:
        globals().update(IPythonWorkspace)
    
    ##### Model fitting bar Functions #####
    '''
    def fitToggle(name, value):
        if value == True:
            fitBar.visible = True
            #time.sleep(clearDelay) # Pause then clear the input field
            #dataVar.value = ''
        else:
            fitBar.visible = False
            #dataVar.value='<variable name>'
        return
    '''
    
    #dataSet = None
    def onDataLoad(name):
        global dataSet
        print('Loading: "', dataVar.value, '"...', end=' ')
        #print(vars())
        #print(globals())
        if dataVar.value in vars(): ### locals()? # http://stackoverflow.com/questions/7969949/whats-the-difference-between-globals-locals-and-vars
            dataSet = vars()[dataVar.value]
            #useFitCheck.value = True
            print('Successfully loaded from vars!')
        elif dataVar.value in globals():
            dataSet = globals()[dataVar.value] # eval(dataVar.value) ##### Is this safe?
            #useFitCheck.value = True
            print('Successfully loaded from globals!')
        else:
            # Try treating it as a file name instead?
            fh = open(path.join(dDir, dataVar.value), "rb") #dDir+'expData'+".pkl"
            dataSet = pickle.load(fh)
            fh.close()
            print('Successfully loaded "{}"!'.format(path.join(dDir, dataVar.value)))
            #dataSet = None
            #useFitCheck.value = False
            #warnings.warn('Warning! Variable: {} not found!'.format(dataVar.value))
        #return dataSet
        
    def onClickFitButton(b): # on_button_clicked
        """Main function to fit a model to a supplied data set"""
        
        #global fitParamsPopup
        #global clearOutput
        if clearOutput.value:
            clear_output()
            #if 'fitParamsPopup' in vars() or 'fitParamsPopup' in globals():
            #    fitParamsPopup.close()
        
        model = str(statesToFitButtons.value)
        mInd = statesDict[model] #statesDict[' '+str(nStates)]
        pSet = modelParams[modelList[mInd]]

        initialParams = getGUIparams(pSet, pfValArr[mInd][:], varyList=fVaryArr[mInd][:], minList=pfMinArr[mInd][:], maxList=pfMaxArr[mInd][:], exprList=fExprArr[mInd][:])
        fittedParams = fitModels(dataSet, nStates=int(statesToFitButtons.value), params=initialParams, postFitOpt=runPostOpt.value, relaxFact=relaxFactWid.value, method=methods[fitMethods.value], postFitOptMethod=methods[postOptFitMethods.value])
        
        #fitParamReport = widgets.TextareaWidget(description='Report:',value=fitRhO.reportParams())
        #fitParamsPopup = widgets.PopupWidget(children=[fitParamReport],button_text='Fitted Parameters',description='Fitted {} state model parameters from: {}'.format(int(statesToFitButtons.value),dataVar.value))
        #display(fitParamsPopup)
        
        #[:]???
        setGUIparams(fittedParams, modelParamsList[model], pfValArr[mInd], varyList=fVaryArr[mInd], minList=pfMinArr[mInd], maxList=pfMaxArr[mInd], exprList=fExprArr[mInd])
        
        # if useFitCheck.value: # Set the run parameters too
            # setGUIparams(fittedParams, modelParamsList[model], pValArr[mInd])
        
        if runSSAcheck.value == True:
            fitRhO = models[modelList[mInd]]()
            fitRhO.updateParams(fittedParams)
            fitRhO.plotRates()
            characterise(fitRhO)
        
        return
    
    def onClickCharacteriseButton(b):
        model = str(statesToFitButtons.value)
        mInd = statesDict[model]
        pSet = modelParams[modelList[mInd]]
        
        fittedParams = getGUIparams(pSet, pfValArr[mInd][:], varyList=fVaryArr[mInd][:], minList=pfMinArr[mInd][:], maxList=pfMaxArr[mInd][:], exprList=fExprArr[mInd][:])
        
        fitRhO = models[modelList[mInd]]()
        fitRhO.updateParams(fittedParams)
        fitRhO.plotRates()
        characterise(fitRhO)
        return
    
    def onClickExportFitButton(b):
        model = str(statesToFitButtons.value)
        mInd = statesDict[model]
        pSet = modelParams[modelList[mInd]]
        
        fittedParams = getGUIparams(pSet, pfValArr[mInd][:], varyList=fVaryArr[mInd][:], minList=pfMinArr[mInd][:], maxList=pfMaxArr[mInd][:], exprList=fExprArr[mInd][:])
        setGUIparams(fittedParams, modelParamsList[model], pValArr[mInd])
        return
    
    

    '''    
    ##### NEURON bar functions #####
    # def onHocLoad(name):
        # print('Loading: "',hocFile.value, '"...', end=' ')
        # try: # Load mechanism and set some appropriate parameters
            # h = hoc.HocObject()
            # h.xopen(hocFile.value)
            # h.load_mechanisms() #h('nrn_load_dll("libnrnmech.so")')
            
            # #from neuron import gui # Load NEURON GUI
            
        # except:
            # print('Error! File: {} not found!'.format(dataVar.value))
        # return
    
        
    ##### Brian bar functions #####
    def onBrianLoad(name):
        print('Loading: "',brianFile.value, '"...', end=' ')
        try: # Load mechanism and set some appropriate parameters
            print('Finish me!!!')
            
        except:
            print('Error! File: {} not found!'.format(dataVar.value))
        return
    '''    

        
    ##### Parameter Bar and Tabs functions #####
    
    def onClickParamsLoad(name):
        global paramSet
        #from IPython import get_ipython # for output widgets
        #with paramOutput:
        print('Loading: "', paramVar.value, '"...', end=' ')
        
        if paramVar.value in vars():
            paramSet = vars()[paramVar.value]
            #with paramOutput:
            print('Successfully loaded from vars!')
            
        elif paramVar.value in globals():
            paramSet = globals()[paramVar.value]
            #with paramOutput:
            print('Successfully loaded from globals!')
            
        else:
            paramSet = None
            #with paramOutput:
            print("Unable to find '{}'!".format(paramVar.value))
            #warnings.warn('Warning! Variable: {} not found!'.format(paramVar.value))
            
        if paramTabs.selected_index == TabGroups['Models']:
            mInd = modelParamsTabs.selected_index
            model = statesArray[mInd]
            #mInd = statesDict[' '+str(nStates)]
            #pSet = modelParams[mInd] # Default model parameters
            setGUIparams(modelParams[model], modelParamsList[model], pValArr[mInd][:]) # setGUImodelParams(paramSet,model)
            
        elif paramTabs.selected_index == TabGroups['Simulators']:
            sInd = simParamsTabs.selected_index
            simulator = simList[sInd]
            setGUIparams(simParams[simulator], simParamsList[simulator], sim_pValArr[sInd][:])
            
        elif paramTabs.selected_index == TabGroups['Protocols']: # Protocols tab selected
            pInd = protParamsTabs.selected_index
            #pInd = protIndDict[protocol]
            protocol = protList[pInd]
            #pSet = protParams[protocol] # Default protocol parameters
            setGUIparams(protParams[protocol], protParamsList[protocol], prot_pValArr[pInd][:]) #setGUIprotParams(paramSet,protocol) #modelParamsList[protocol]
            
        else: 
            raise ValueError('Unknown Tab Index!')
            
        return paramSet

        
    def onClickParamsReset(b): # on_button_clicked
        
        if paramTabs.selected_index == TabGroups['Models']:
            mInd = modelParamsTabs.selected_index
            model = statesArray[mInd]
            #mInd = statesDict[' '+str(nStates)]
            #pSet = modelParams[model] # Default model parameters
            #pSet = modelParams[statesArray[mInd]]
            setGUIparams(modelParams[model], modelParamsList[model], pValArr[mInd][:]) #setGUImodelParams(pSet,model)
            
        elif paramTabs.selected_index == TabGroups['Simulators']:
            sInd = simParamsTabs.selected_index
            simulator = simList[sInd]
            setGUIparams(simParams[simulator], simParamsList[simulator], sim_pValArr[sInd][:])
            
        elif paramTabs.selected_index == TabGroups['Protocols']:
            pInd = protParamsTabs.selected_index
            #pInd = protIndDict[protocol]
            protocol = protList[pInd]
            #pSet = protParams[protocol] # Default protocol parameters
            setGUIparams(protParams[protocol], protParamsList[protocol], prot_pValArr[pInd][:]) #setGUIprotParams(pSet,protocol)
            
        else: # Protocols tab selected
            raise ValueError('Unknown Tab Index!')
            
        return
    


    
    ### Utility functions 

    #getGUIparams(modelParams[model],pValArr[modelParams.keys().index(model)][:])
    #getGUIparams(simParams[simulator],sim_pValArr[simParams.keys().index(simulator)][:])
    #getGUIparams(protParams[protocol],prot_pValArr[protParams.keys().index(protocol)][:])
    
    # pSet = protParams[protocol]
    # widgetList = prot_pValArr[pInd][:]
    #getGUIparams(simParams[simulator],sim_pValArr[sInd][:])
    
    def getGUIparams(pSet, valueList, varyList=None, minList=None, maxList=None, exprList=None):
        # pSet must be the default parameters used to build the GUI
        userParams = copy.deepcopy(pSet) #Parameters()
        #print(userParams)
        #if isinstance(userParams, PyRhOparameters):
        for i, param in enumerate(userParams): #enumerate(pSet):
            #print(i, param)
            if isinstance(userParams[param].value, (list, tuple)): ################################## Change to not number!!!
                #userParams.add(param, value=ast.literal_eval(valueList[i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
                userParams[param].set(value=ast.literal_eval(valueList[i].value))
            else:
                #userParams.add(param, value=valueList[i].value)
                userParams[param].set(value=valueList[i].value)
            if varyList is not None:
                userParams[param].set(vary=varyList[i].value) #userParams[param].vary = varyList[i].value
            if minList is not None:
                userParams[param].set(min=minList[i].value)
            if maxList is not None:
                userParams[param].set(max=maxList[i].value)
            if exprList is not None:
                userParams[param].set(expr=exprList[i].value)
                
        return userParams
    
    def setGUIparams(pSet, paramList, valueList, varyList=None, minList=None, maxList=None, exprList=None):
        #for i, param in enumerate(pSet): # Changed to allow for pSet with arbitrary orders 
        for param in pSet:
            if param in paramList:
                i = paramList.index(param)
                
                if isinstance(pSet[param].value, (list, tuple)): ################################## Change to not number!!!
                    valueList[i].value = str(pSet[param].value)
                else:
                    valueList[i].value = pSet[param].value
                if varyList is not None:
                    varyList[i].value = pSet[param].vary
                if minList is not None:
                    minList[i].value = pSet[param].min
                if maxList is not None:
                    maxList[i].value = pSet[param].max
                if exprList is not None and pSet[param].expr is not None:
                    exprList[i].value = pSet[param].expr ### Change units handling
    
    
        
    ##### Main Run function #####
    
        #saveData = True
        #@interact(nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'delta', 'rectifier', 'shortPulse', 'recovery'), saveData=True, verbose=1)

    def runModel(model, protocol, simulator='Python', saveData=True, verbose=1): #verboseSlide.value
        """Main GUI function to create protocol, simulator and rhdopsin objects, set parameters, run and plot"""
        
        #verbose = verbose
        #nStates = int(model)
        # print("\nRunning Protocol '{}' on the {}-state model...".format(protocol,nStates))
        # print('--------------------------------------------------------------------------------\n')
        
        ### Choose between 3, 4, & 6 state models
        #if useFitCheck.value: # fitButton.value and 
        #    RhO = fitRhO
            ### To set GUI parameters, widget arrays should be dictionaries so subsets of parameters can be set in arbitrary order!!!
            #setGUIprotParams(expProtParams,'custom')            
        #else:
        RhO = selectModel(int(model))
        userModelParams = Parameters()
        mInd = statesDict[model] #statesDict[' '+str(nStates)]
        pSet = modelParams[modelList[mInd]]
        userModelParams = getGUIparams(pSet, pValArr[mInd][:])
        RhO.setParams(userModelParams)
        
        ### Get Protocol Parameters
        userProtParams = PyRhOparameters()
        pInd = protList.index(protocol) #protIndDict[protocol]
        userProtParams = getGUIparams(protParams[protocol], prot_pValArr[pInd][:])
        Prot = protocols[protocol](userProtParams)
        
        if protocol is 'custom': # and custPulseGenInput.value is not None:
            Prot.phi_ft = custPulseGenLoad(custPulseGenInput.value)
        
        ### Get Simulator Parameters
        userSimParams = PyRhOparameters()
        sInd = simList.index(simulator)
        userSimParams = getGUIparams(simParams[simulator], sim_pValArr[sInd][:])
        Sim = simulators[simulator](Prot, RhO, userSimParams)
        
        Sim.run(verbose)
        if verbose > 0: #saveData:
            Sim.plot()
        
        # print("\nFinished!")
        # print('================================================================================\n\n')
        
        return #Prot, RhO
    

    
    ##### Run bar #####
    
    ##### Run Bar Functions #####

    
    ### Protocol Dropdown
    def onChangeProtDropdown(change):#(name,value):
        paramTabs.selected_index = TabGroups['Protocols'] # Set to protocol parameters
        protParamsTabs.selected_index = protList.index(change['new'])#(value) # Get the index of the selected protocol #protIndDict[value]
    
    protDropdown = widgets.Dropdown(options=protList, value='step')
    #protDropdown.on_trait_change(onChangeProtDropdown, 'value')
    protDropdown.observe(onChangeProtDropdown, names='value')
    
    
    ### Run Button
    def onClickRunButton(b): # on_button_clicked
        if clearOutput.value:
            clear_output()
        runModel(stateButtons.value, protDropdown.value, simDropdown.value, saveButton.value, verboseSlide.value) #int(stateButtons.value)
        return
        
    runButton = widgets.Button(description="Run!")
    runButton.on_click(onClickRunButton)
    
    
    ### Model states button
    def onChangeModel(change): #(name,value):
        paramTabs.selected_index = TabGroups['Models'] # Set to model parameters
        modelParamsTabs.selected_index = statesDict[change['new']] #[value]
    
    stateButtons = widgets.ToggleButtons(description='Model states: ', options=statesArray, ) #https://github.com/ipython/ipython/issues/6469 #'Model states: '
    #stateButtons.on_trait_change(onChangeModel, 'value')#stateButtons.value_name  #'_view_name'
    stateButtons.observe(onChangeModel, names='value')
    
    ### Create states buttons link
    #modelParamsLink = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    
    ### Parameters Toggle
    def onClickParamsToggle(change): #(name, value): #paramsToggle
        if change['old'] is False: #value == True: # Set to current model and protocol tabs
            paramsControlBar.visible = True
            paramTabs.visible = True
            modelParamsTabs.selected_index = statesDict[stateButtons.value]
            protParamsTabs.selected_index = protList.index(protDropdown.value) #protIndDict[protDropdown.value]
        else:
            paramsControlBar.visible = False
            paramTabs.visible = False
        return
    paramsButton = widgets.ToggleButton(description='Parameters', value=False)
    #paramsButton.on_trait_change(onClickParamsToggle, 'value')
    paramsButton.observe(onClickParamsToggle, names='value')
    
    # ### Load NEURON parameters
    # hocFile = widgets.Text(description='HOC file: ',placeholder='<file name>')
    # hocLoad = widgets.Button(description='Load')
    # hocLoad.on_click(onHocLoad)
    
    # ### Load Brian parameters
    # brianFile = widgets.Text(description='Brian file: ',placeholder='<file name>')
    # brianLoad = widgets.Button(description='Load')
    # brianLoad.on_click(onBrianLoad)

    
    ### Output format Dropdown
    def onChangeSimDropdown(change): #(name,value):
        # if value == 'NEURON':
            # NEURONbox.visible=True
            # Brianbox.visible=False
        # elif value == 'Brian':
            # Brianbox.visible=True
            # NEURONbox.visible=False
        # else: # Set both to invisible
            # NEURONbox.visible=False
            # Brianbox.visible=False
        paramTabs.selected_index = TabGroups['Simulators'] # Set to simulator parameters
        simParamsTabs.selected_index = simList.index(change['new']) #(value) # Get the index of the selected protocol #protIndDict[value]
        return
    simDropdown = widgets.Dropdown(options=['Python', 'NEURON', 'Brian'], value='Python')
    #simDropdown.on_trait_change(onChangeSimDropdown, 'value')
    simDropdown.observe(onChangeSimDropdown, names='value')
    
    
    
    
    ### Plot and Save button
    saveButton = widgets.Checkbox(description='Plot & Save', value=True)
    saveButton.visible = False

    ### Verbosity slider
    def verboseChange(change): #(value):
        global verbose
        verbose = change['new'] #value
        return
    #verboseSlide = widgets.FloatProgressWidget(value=1, min=-1, max=3, step=1, description='Verbosity:')
    verboseSlide = widgets.IntSlider(value=1, min=0, max=3, description='Output:')
    #verboseSlide.on_trait_change(verboseChange, 'value')
    verboseSlide.observe(verboseChange, names='value')
    verboseSlide.visible = True #False

    ### Clear ouput checkbox
    clearOutput = widgets.Checkbox(description='Clear', value=True)
    
    
    #### Run Bar container
    runBar = widgets.HBox(children=[stateButtons,protDropdown,simDropdown,paramsButton,saveButton,verboseSlide,clearOutput,runButton]) #fitButton,
    #####display(runBar) # Commented to nest in GUI box
    

    
    
    ### Set Button styles after displaying the runBar
    # paramsButton.button_style = 'info'
    # stateButtons.button_style = 'success'
    # stateButtons.margin = '5px'
    # simDropdown.button_style = 'success'
    # simDropdown.margin = '5px'
    # protDropdown.button_style = 'success'
    # protDropdown.margin = '5px'
    ###fitButton.button_style = 'success'
    # runButton.button_style = 'danger'
    # runButton.margin = '5px'
    
    
    #stateButtons.add_class('btn-primary')
    #protDropdown.add_class('btn-success')
    #simDropdown.add_class('btn-warning')
    
    ### button bootstrap styles: default, primary, success, info, warning, danger, link
    

    # runBar.align = 'center'
    # verboseSlide.width = '60px'
    
    ##### Parameters control bar #####
    
    paramVar = widgets.Text(description='Parameter Set: ', placeholder='<variable name>')
    paramLoad = widgets.Button(description='Load')
    paramLoad.on_click(onClickParamsLoad)
    paramOutput = widgets.Output()
    
    paramBox = widgets.HBox(children=[paramVar])
    
    ### Reset Parameters Button
    paramsResetButton = widgets.Button(description="Reset")
    paramsResetButton.on_click(onClickParamsReset)
    
    ### Create Parameters control bar
    paramsControlBar = widgets.HBox(children=[paramBox,paramLoad,paramsResetButton,paramOutput]) #paramVar
    #####display(paramsControlBar) # Commented to nest in GUI box
    paramsControlBar.visible = False
    
    
    ### Set Parameters Bar style
    # paramVar.width = '150px' #set_css({'width': '150px'})
    # paramLoad.button_style = 'warning' #add_class('btn-warning')
    # paramsResetButton.button_style = 'warning' #add_class('btn-warning')
    
    # paramsControlBar.align = 'center'
    
    
    
    
    ##### Create Parameters Tab widgets #####
    
    ### lmfit Parameters() are ordereddict objects (https://github.com/lmfit/lmfit-py/blob/master/lmfit/parameter.py) but it might be better to change pValArr and prot_pValArr etc. to ordereddicts too...
    
    
    
    
    
    
    

    
    
    
    
    
    ##### Fit Data parameters #####
    
    ### Create Data set entry
    dataVar = widgets.Text(placeholder='<Data Set>') #description='Data Set: ',
    dataLoad = widgets.Button(description='Load')
    dataLoad.on_click(onDataLoad)
    
    dataBox = widgets.HBox(children=[dataVar]) #dataLoad
    
    #rArrow = widgets.Latex(value='$\\Rightarrow$')
    
    ### Create Fit Model States buttons
    def onChangeModelToFit(change): #(name,value):
        #paramTabs.selected_index = TabGroups['Models'] # Set to model parameters
        fitParamsTabs.selected_index = statesDict[change['new']] #[value]
    statesToFitButtons = widgets.ToggleButtons(description='Model: ', options=statesArray) #,value=u' 3') #https://github.com/ipython/ipython/issues/6469
    #statesToFitButtons.on_trait_change(onChangeModelToFit, 'value')
    statesToFitButtons.observe(onChangeModelToFit, names='value')
    
    #fitLabel = widgets.HTML(value='Models to fit: ')
    #fit3sCheck = widgets.Checkbox(description='3 state', value=True)
    #fit4sCheck = widgets.Checkbox(description='4 state', value=False)
    #fit6sCheck = widgets.Checkbox(description='6 state', value=False)
    
    
    ### Hyperparameters
    #p0fVvar = widgets.Text(value=str(p0fV), description='p0fV: ', width='150px')
    #p0IPIvar = widgets.Text(value=str(p0IPI), description='p0IPI: ', width='150px')
    
    #hyperTitle = widgets.HTML(value='<h4>General Parameters</h4>')
    #hyperParams = widgets.VBox(children=[hyperTitle, p0fVvar, p0IPIvar])#, border_color='blue')
    
    
    ### Create Checkboxes
    runSSAcheck = widgets.Checkbox(description='Characterise', value=False)
    #useFitCheck = widgets.Checkbox(description='Use fit', value=False)
    relaxFactWid = widgets.FloatSlider(value=2.0, min=1.0, max=10.0, step=0.1, description='Relax factor:')
    fitMethods = widgets.Dropdown(options=fitMethodsDict, value=fitMethodsDict[defMethod]) #, description='Method:')
    #plotExpData = widgets.Checkbox(description='Plot', value=True)
    def onTickChange(change):
        postOptFitMethods.visible = change['new']
    runPostOpt = widgets.Checkbox(description='Post-fit opt.', value=True)
    runPostOpt.observe(onTickChange, names='value')
    postOptFitMethods = widgets.Dropdown(options=fitMethodsDict, value=fitMethodsDict[defMethod]) #, description='Post-fit Meth:')
    
    ### Create Run Fit Button
    runFitButton = widgets.Button(description="Fit!")
    runFitButton.on_click(onClickFitButton)
    characteriseButton = widgets.Button(description="Characterise")
    characteriseButton.on_click(onClickCharacteriseButton)
    exportFitButton = widgets.Button(description="Export")
    exportFitButton.on_click(onClickExportFitButton)
    
    ### Create Fit Bar
    #fitBar = widgets.HBox(children=[dataVar, dataLoad, statesToFitButtons, runSSAcheck, useFitCheck, plotExpData, runFitButton]) #fit3sCheck, fit4sCheck, fit6sCheck
    #display(fitBar)
    #fitBar.visible=False
    
    # ### Set formatting
    # dataVar.width = '150px' #set_css({'width': '150px'})
    # #dataVar.margin = '0px'
    # dataLoad.button_style = 'warning' #add_class('btn-warning')
    # runFitButton.button_style = 'danger' #add_class('btn-danger')
    # #runFitButton.margin = '5px'
    # characteriseButton.button_style = 'success'
    # # Set Fit Bar formatting after display
    # dataLoad._dom_classes = ('margin-left','10px')
    # #fitBar.align = 'center'
    # exportFitButton.button_style = 'info'
    
    # runSSAcheck.margin = '5px'
    # #plotExpData.margin = '5px'
    # runPostOpt.margin = '5px'
    # runFitButton.margin = '5px'
    # #characteriseButton.margin = '5px'
    # exportFitButton.margin = '5px'
    
    #dataLoadBar = widgets.HBox(children=[dataVar, dataLoad, useFitCheck], align='center')
    ##statesBar = widgets.HBox(children=[fitLabel, fit3sCheck, fit4sCheck, fit6sCheck], align='center') #.align = 'center'
    #statesBar = widgets.HBox(children=[statesToFitButtons], align='center')
    

    
    #fitOutputBar = widgets.HBox(children=[runSSAcheck, plotExpData, runFitButton], align='center')
    #runSSAcheck
    fitBarMain = widgets.HBox(children=[dataBox, dataLoad, statesToFitButtons, runFitButton, characteriseButton, exportFitButton]) #, align='baseline')
    fitBarMeta = widgets.HBox(children=[fitMethods, runPostOpt, postOptFitMethods, relaxFactWid]) #, align='baseline')
    fitBar = widgets.VBox(children=[fitBarMain, fitBarMeta]) #, align='start') #align='center') #plotExpData #dataVar, dataLoad
    #fitParamsHead = widgets.VBox(children=[dataLoadBar, statesBar, hyperParams, fitOutputBar]) # Left side container
    fitParamsHead = widgets.VBox(children=[fitBar]) # , hyperParams
    #fitNotesBoxes = widgets.VBox(children=[]) # Right side container for figure and equations
    
    #fitParamsTabs = widgets.HBox(children=[fitBox]) # fitNotesBoxes # Hack to have a parent tab with no children
    #fitParamsBox = widgets.VBox(children=[fitBox,fitParamsTabs])
    
    ### Create parameter fitting entries
    
    modelFitBoxes = [None for m in range(len(modelParams))]
    modelFitParamBoxes = [None for m in range(len(modelParams))] # Left side container
    modelFitNotesBoxes = [None for m in range(len(modelParams))] # Right side container for figure and equations
    #figHTML = [None for m in range(len(modelParams))] 
    #eqBox = [None for m in range(len(modelParams))]
    pfLabArr = [[None for p in modelParams[m]] for m in modelList]
    pfMinArr = [[None for p in modelParams[m]] for m in modelList]
    pfValArr = [[None for p in modelParams[m]] for m in modelList] # Array of parameter values #range(len(modelParams))
    pfMaxArr = [[None for p in modelParams[m]] for m in modelList]
    fUnitArr = [[None for p in modelParams[m]] for m in modelList] # Array of units #[[None for p in modelParamSet] for modelParamSet in modelParams]
    fVaryArr = [[None for p in modelParams[m]] for m in modelList]
    fExprArr = [[None for p in modelParams[m]] for m in modelList]
    pfBoxArr = [[None for p in modelParams[m]] for m in modelList] # Array of parameter boxes ###[[None for p in modelParams[modelList[m]]] for m in range(len(modelParams))]
    
    spacer = widgets.HTML()
    ###spacer.width = '150px' # Equal to the width of a drop-drown menu
    
    #widgets.HTML(value='<strong>Initial Value</strong>')
    #colHeadings = widgets.HTML(value='<table><tr><th width="40px"></th><th width=60px align="right">Vary</th><div align="center"><th width=150px align="center">Minimum</th><th width=150px align="center">Initial</th><th width=150px align="center">Maximum</th><th width=150px align="center">Units</th><th width=150px align="center">Expression</th></div></tr></table>')
    # style="width:800px" <th>Parameter</th>
    
    colHeadings = widgets.HBox(children = [widgets.HTML(value='<p align="right">Vary</p>', width='80px'), widgets.HTML(value='<p align="center">Minimum</p>', width='150px'), widgets.HTML(value='<p align="center">Value</p>', width='150px'), widgets.HTML(value='<p align="center">Maximum</p>', width='150px'), widgets.HTML(value='<p align="center">Units</p>', width='150px'), widgets.HTML(value='<p align="center">Expression</p>', width='150px')])
    for model, pSet in modelParams.items(): #range(len(modelParams)):
        #pSet = modelParams[m]            # Set of n-state model parameters
        m = modelList.index(model)
        nParams = len(pSet)#.keys())
        
        i=0
        for key in pSet:#.keys(): #, value in pSet.items():
            # if isinstance(pSet[key].value, list): # list ==> Text
                # pluginValues[pluginInd][pInd] = widgets.Text(value=str(pSet[key].value), description=key) # np.asarray
            # if isinstance(pSet[key].value, bool): # Allow model features to be turned on or off
                # pValArr[m][i] = widgets.Dropdown(options=boolDict,value=pSet[key].value,description=key)
            # elif isinstance(pSet[key].value, numbers.Number): # Number: (int, long, float, complex)
            #if pSet[key].min == None or pSet[key].max == None:
            #    pValArr[m][i] = widgets.FloatText(value=pSet[key].value, description=key)#"{} [{}]:".format(key, pSet[key].expr))
            #else:
            
            if key in modelLabels:
                name = '$'+modelLabels[key]+'$'
            else:
                name = key
            
            pfLabArr[m][i] = widgets.Latex(value=name)
            #pfLabArr[m][i].width = '50px'#'20px'
            minVal = pSet[key].min if pSet[key].min != None else -np.inf
            pfMinArr[m][i] = widgets.FloatText(value=minVal)#,description='min') #str(minVal)
            maxVal = pSet[key].max if pSet[key].max != None else np.inf
            pfMaxArr[m][i] = widgets.FloatText(value=maxVal)#,description='max') #str(maxVal)
            pfValArr[m][i] = widgets.BoundedFloatText(value=pSet[key].value, min=minVal, max=maxVal)#, description='initial')#"{} [{}]:".format(key, pSet[key].expr))
            #pfValArr[m][i].width = '150px'
            
            fitExpr = pSet[key].expr if pSet[key].expr is not None else ''
            # if pSet[key].expr is None:
                # fExprArr[m][i] = widgets.Text(value='')#, description='expr')
            # else:
                # fExprArr[m][i] = widgets.Text(value=pSet[key].expr)
            fExprArr[m][i] = widgets.Text(value=fitExpr)
            ###fExprArr[m][i].width = '150px'
            #fVaryArr[m][i] = widgets.Checkbox(value=True)#, description='vary')
            fVaryArr[m][i] = widgets.Checkbox(value=pSet[key].vary)#, description='vary')
            ###fVaryArr[m][i].width = '30px'
            
            if key in unitLabels: #pSet[key].expr is not None:
                #fUnitArr[m][i] = widgets.Dropdown(options=[pSet[key].expr])#,value=pSet[key].expr) ### Change units handling
                fUnitArr[m][i] = widgets.Dropdown(options=[unitLabels[key]])#,value=pSet[key].expr) ### Change units handling
                pfBoxArr[m][i] = widgets.HBox(children=[pfLabArr[m][i],fVaryArr[m][i],pfMinArr[m][i],pfValArr[m][i],pfMaxArr[m][i],fUnitArr[m][i],fExprArr[m][i]], align='center')
            else:
                
                pfBoxArr[m][i] = widgets.HBox(children=[pfLabArr[m][i],fVaryArr[m][i],pfMinArr[m][i],pfValArr[m][i],pfMaxArr[m][i],spacer,fExprArr[m][i]], align='center')
            i+=1
        
        modelFitParamBoxes[m] = widgets.VBox(children=pfBoxArr[m]) #Box
        #modelFitNotesBoxes[m] = widgets.HBox(children=[])#[figHTML[m],eqBox[m]])
        
        #modelNotesBoxes[m].add_class('box-flex1')
        #modelBox = widgets.HBox(children=[modelParamBoxes,modelNotesBoxes])
        modelFitBoxes[m] = widgets.VBox(children=[colHeadings, modelFitParamBoxes[m]]) #, #,modelFitNotesBoxes[m]])#modelBox
        ###modelFitBoxes[m].margin = '5px'
        
    ### Linked parameters
    #E_Link = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    
    #modelParamsTabs = widgets.Tab(description='Parameter Settings', children=modelBoxes, values=statesArray)
    #modelParamsTabs.margin = '5px'
    #modelParamsTabs.on_trait_change(onChangeModelTab,'selected_index')
    #modelParamsTabs.on_displayed(setModelParamsTabs,'selected_index')
    #modelParamsTabs.on_displayed()
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    
    #fitParamsTabs = buildLayerTab(modelFitBoxes)
    
    # fitParamsTabs.set_title(0, 'Three-state model') # Moved to loop over modelTitles
    # fitParamsTabs.set_title(1, 'Four-state model')
    # fitParamsTabs.set_title(2, 'Six-state model')
    ###fitParamsTabs.width = '800px'
    
    #modelFitTitle = widgets.HTML(value='<h4>Model Specific Parameters</h4>')
    #fitBox = widgets.VBox(children=[fitParamsHead,modelFitTitle,fitParamsTabs]) #HBox
    
    ###fitBox.margin = '5px'
    #modelParamsTabs = widgets.Tab(description='Parameter Settings', children=modelBoxes, values=statesArray)
    #modelParamsTabs.margin = '5px'
        
    def onChangeFitTab(change):#(name,value):
        statesToFitButtons.value = statesArray[change['new']] #[value]
    fitParamsTabs = widgets.Tab(description='Parameter Settings', children=modelFitBoxes, values=statesArray)
    #fitParamsTabs.on_trait_change(onChangeFitTab, 'selected_index')
    fitParamsTabs.observe(onChangeFitTab, names='selected_index')
    fitBox = widgets.VBox(children=[fitParamsHead, fitParamsTabs]) #HBox
    
    #statesToFitButtons.button_style = 'info'
    #statesToFitButtons.margin = '5px'
    
    
    
    
    
    
    
    
    #import numbers
    
    ##### Model parameters #####
    modelBoxes = [None for m in range(len(modelParams))]
    modelParamBoxes = [None for m in range(len(modelParams))] # Left side container
    modelNotesBoxes = [None for m in range(len(modelParams))] # Right side container for figure and equations
    figHTML = [None for m in range(len(modelParams))] 
    eqBox = [None for m in range(len(modelParams))]
    pValArr = [[None for p in modelParams[m]] for m in modelList] # Array of parameter values #range(len(modelParams))
    unitArr = [[None for p in modelParams[m]] for m in modelList] # Array of units #[[None for p in modelParamSet] for modelParamSet in modelParams]
    pBoxArr = [[None for p in modelParams[m]] for m in modelList] # Array of parameter boxes ###[[None for p in modelParams[modelList[m]]] for m in range(len(modelParams))]
    for model, pSet in modelParams.items(): #range(len(modelParams)):
        #pSet = modelParams[m]            # Set of n-state model parameters
        m = modelList.index(model)
        nParams = len(pSet)#.keys())
        
        i=0
        for key in pSet:#.keys(): #, value in pSet.items():
            # if isinstance(pSet[key].value, list): # list ==> Text
                # pluginValues[pluginInd][pInd] = widgets.Text(value=str(pSet[key].value), description=key) # np.asarray
            # if isinstance(pSet[key].value, bool): # Allow model features to be turned on or off
                # pValArr[m][i] = widgets.Dropdown(options=boolDict,value=pSet[key].value,description=key)
            # elif isinstance(pSet[key].value, numbers.Number): # Number: (int, long, float, complex)
            if key in modelLabels:
                name = '$'+modelLabels[key]+'$'
            else:
                name = key
            if pSet[key].min == (None or -np.inf) or pSet[key].max == (None or np.inf):
                pValArr[m][i] = widgets.FloatText(value=pSet[key].value, description=name)#"{} [{}]:".format(key, pSet[key].expr))
            else:
                pValArr[m][i] = widgets.BoundedFloatText(value=pSet[key].value, min=pSet[key].min, max=pSet[key].max, description=name)#"{} [{}]:".format(key, pSet[key].expr))
            
            if key in unitLabels: #pSet[key].expr is not None:
                #unitArr[m][i] = widgets.Dropdown(options=[pSet[key].expr],value=pSet[key].expr) ### Change units handling
                unitArr[m][i] = widgets.Dropdown(options=[unitLabels[key]],value=unitLabels[key]) ### Change units handling
                pBoxArr[m][i] = widgets.HBox(children=[pValArr[m][i],unitArr[m][i]])
            else:
                pBoxArr[m][i] = widgets.HBox(children=[pValArr[m][i]])
            i+=1
        
        guiFigDir = GUIdir + '/' #'gui/'
        
        figHTML[m] = widgets.HTML()
        eqBox[m] = widgets.Latex()
        #figPrefix = os.path.join(pyrhoPath, 'gui')
        
        if int(statesArray[m]) == 3: # m==0
            #figHTML[m].value = '<img src="{}" alt="Three state model" width=150>'.format('file://'+os.path.join(figPrefix, '3state_model.png'))
            #figHTML[m].value = '<img src="3state_model.png" alt="Three state model" width=150>'
            figHTML[m].value = '<img src="{}{}state_model.png" alt="{}-state model" width=150>'.format(guiFigDir, statesArray[m], stateLabs[statesArray[m]])
            eqBox[m].value = models[model].equations
            
        elif int(statesArray[m]) == 4: # m==1
            figHTML[m].value = '<img src="{}{}state_model.png" alt="{}-state model" width=180>'.format(guiFigDir, statesArray[m], stateLabs[statesArray[m]])
            eqBox[m].value = models[model].equations
            
        else: #int(statesArray[m]) == 6:
            figHTML[m].value = '<img src="{}{}state_model.png" alt="{}-state model" width=210>'.format(guiFigDir, statesArray[m], stateLabs[statesArray[m]])
            eqBox[m].value = models[model].equations
            
        modelParamBoxes[m] = widgets.Box(children=pBoxArr[m])
        modelNotesBoxes[m] = widgets.HBox(children=[figHTML[m],eqBox[m]])
        #modelNotesBoxes[m].add_class('box-flex1')
        #modelBox = widgets.HBox(children=[modelParamBoxes,modelNotesBoxes])
        modelBoxes[m] = widgets.HBox(children=[modelParamBoxes[m],modelNotesBoxes[m]])#modelBox
        ###modelBoxes[m].margin = '5px'
    
    ### Linked parameters
    #E_Link = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    
    
    
    
    
    def onChangeModelTab(change): #(name,value):
        stateButtons.value = statesArray[change['new']] #[value] #' '+str(value)
    modelParamsTabs = widgets.Tab(description='Parameter Settings', children=modelBoxes, values=statesArray)
    ###modelParamsTabs.margin = '5px'
    #modelParamsTabs.on_trait_change(onChangeModelTab, 'selected_index')
    modelParamsTabs.observe(onChangeModelTab, names='selected_index')
    #modelParamsTabs.on_displayed(setModelParamsTabs,'selected_index')
    #modelParamsTabs.on_displayed()
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    
    
    
    
    ##### Simulator parameters #####
    
    simBoxes = [None for sim in simParams]
    simParamBoxes = [None for sim in simParams] # Left side container
    simNotesBoxes = [None for sim in simParams] # Right side container for figure and equations
    #simStimHTML = [None for sim in simParams]
    simFigHTML = [None for sim in simParams]
    #eqBox = [None for m in range(len(modelParams))]
    sim_pValArr = [[None for p in simParams[sim]] for sim in simParams] # Array of parameter values
    sim_unitArr = [[None for p in simParams[sim]] for sim in simParams] # Array of units
    sim_pBoxArr = [[None for p in simParams[sim]] for sim in simParams] # Array of parameter boxes
    #sim_notesArr = [[None for p in simParams[sim]] for sim in simParams] # Array of notes boxes
    
    for sInd, sim in enumerate(simParams):
        pSet = simParams[sim]
        #i=0
        for i, param in enumerate(pSet): #param in pSet:#.keys(): #, value in pSet.items():
            if isinstance(pSet[param].value, list):
                sim_pValArr[sInd][i] = widgets.Text(value=str(pSet[param].value), description=param) # np.asarray
            elif isinstance(pSet[param].value, str):
                sim_pValArr[sInd][i] = widgets.Text(value=str(pSet[param].value), description=param)
            elif isinstance(pSet[param].value, bool):
                sim_pValArr[sInd][i] = widgets.Dropdown(options=boolDict, value=pSet[param].value, description=param)
            else:
                if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf):
                    sim_pValArr[sInd][i] = widgets.FloatText(value=pSet[param].value, description=param)
                else:
                    sim_pValArr[sInd][i] = widgets.BoundedFloatText(value=pSet[param].value, min=pSet[param].min, max=pSet[param].max, description=param)
                ###sim_pValArr[sInd][i].width = '150px'
            if param in simUnitLabels:
                sim_unitArr[sInd][i] = widgets.Dropdown(options=[simUnitLabels[param]], value=simUnitLabels[param])
            else:
                sim_unitArr[sInd][i] = widgets.HTML(value='') # Spacer
                #sim_unitArr[sInd][i].width = '150px'
            #sim_notesArr[sInd][i] = widgets.HTML(value=' '+simParamNotes[param]) #widgets.Latex(value='$'+simParamNotes[param]+'$')
            # if param in simUnitLabels: #pSet[param].expr is not None:
                # sim_unitArr[sInd][i] = widgets.Dropdown(options=[simUnitLabels[param]],value=simUnitLabels[param]) #pSet[param].expr ### Change units handling
                # sim_pBoxArr[sInd][i] = widgets.HBox(children=[sim_pValArr[sInd][i],sim_unitArr[sInd][i]])
            # else:
                # sim_pBoxArr[sInd][i] = widgets.HBox(children=[sim_pValArr[sInd][i]])
            sim_pBoxArr[sInd][i] = widgets.HBox(children=[sim_pValArr[sInd][i], sim_unitArr[sInd][i]]) #, sim_notesArr[sInd][i]])
            
            #i+=1
            
        simFigHTML[sInd] = widgets.HTML()
        if sim is 'Brian': ### HACK!
            simFigHTML[sInd].value = 'Brian has not yet been implemented in the GUI - please use the interpreter to run network level simulations. '
        #exampleProt = '{}{}6s.{}'.format(fDir,prot,'png')#saveFigFormat)
        #if os.path.isfile(exampleProt):
        #    protFigHTML[pInd].value='<img src="{}" alt=Example {} width=200px>'.format(exampleProt,prot)
        #else:
        #    protFigHTML[pInd].value='Example Figure'
        simParamBoxes[sInd] = widgets.Box(children=sim_pBoxArr[sInd])
        simNotesBoxes[sInd] = widgets.HBox(children=[simFigHTML[sInd]])# simStimHTML[sInd]  , ])#[figHTML[prot],eqBox[prot]])
        
        
        simBoxes[sInd] = widgets.HBox(children=[simParamBoxes[sInd],simNotesBoxes[sInd]])#modelBox
        #display(protBoxes[pInd])
        #simBoxes[sInd].margin = '5px'
        
    
    ##### Simulator parameters tab #####
    def onChangeSimTab(change): #(name,value):
        simDropdown.value = simList[change['new']] #[value] #simList.index(value)# simDropdown.selected_label #simParamsTabs._titles[value] #[protParamsTabs.selected_index]
    simParamsTabs = widgets.Tab(description='Simulator Settings', children=simBoxes)# \
    ###simParamsTabs.margin = '5px'
    #simParamsTabs.on_trait_change(onChangeSimTab, 'selected_index')
    simParamsTabs.observe(onChangeSimTab, names='selected_index')
    
    
    
    ##### Protocol parameters #####
    
    protBoxes = [None for prot in protList]
    protParamBoxes = [None for prot in protList] # Left side container
    protNotesBoxes = [None for prot in protList] # Right side container for figure and equations
    protStimHTML = [None for prot in protList]
    protFigHTML = [None for prot in protList] 
    protHeaders = [None for prot in protList]
    protBoxesWrapper = [None for prot in protList]
    
    #eqBox = [None for m in range(len(modelParams))]
    prot_pValArr = [[None for p in protParams[prot]] for prot in protList] # Array of parameter values
    prot_unitArr = [[None for p in protParams[prot]] for prot in protList] # Array of units
    prot_pBoxArr = [[None for p in protParams[prot]] for prot in protList] # Array of parameter boxes
    ### These should be included in prot_pBoxArr but are put in th eprotNotesBoxes to avoid alignment issues
    protParamNotesWid = [[None for p in protParams[prot]] for prot in protList] # Array of parameter note boxes
    #protParamNotesBoxes = [None for prot in protList] # Array of boxes for parameter note boxes
    
    ### Create Irradiance <-> Flux converter
    def onClickIrradButton(b):
        fluxVal.value = irrad2flux(E=irradVal.value, lam=lamVal.value)
    
    def onClickFluxButton(b):
        irradVal.value = flux2irrad(phi=fluxVal.value, lam=lamVal.value)
    
    irradVal = widgets.FloatText(value=1) #description='Irradiance', 
    irradUnits = widgets.Dropdown(options=['mW/mm^2'], value='mW/mm^2')
    lamVal = widgets.FloatText(description='$\lambda$', value=470)
    lamUnits = widgets.Dropdown(options=['nm'], value='nm')
    irradButton = widgets.Button(description="==>")
    irradButton.on_click(onClickIrradButton)
    irradButton.tooltips = "Calculate Flux"
    fluxVal = widgets.FloatText(value=1e17) #description='Flux', 
    fluxUnits = widgets.Dropdown(options=['ph./mm^2/s'], value='ph./mm^2/s')
    fluxButton = widgets.Button(description="<==")
    fluxButton.on_click(onClickFluxButton)
    fluxButton.tooltips = "Calculate Irradiance"
    protParamsHead = widgets.HBox(children=[irradVal, irradUnits, fluxButton, lamVal, lamUnits, irradButton, fluxVal, fluxUnits])
    
    
    spacerProtDD = widgets.HTML(value='')
    spacerProtNum = widgets.HTML(value='')
    spacerNull = widgets.HTML(value='')
    
    
    def custPulseGenLoad(funcName):
        #global func
        print("Loading custom pulse generator function: '{}'... ".format(custPulseGenInput.value), end='')
        #import pdb
        #pdb.set_trace()
        # Use whos to list variables in IPython interactive workspace
        if custPulseGenInput.value in vars():
            func = vars()[custPulseGenInput.value]
            print('Successfully loaded from vars!')
            
        #elif custPulseGenInput.value in dir():
        #    func = dir()[custPulseGenInput.value]
        #    print('Successfully loaded from dir!')
            
        elif custPulseGenInput.value in locals():
            func = locals()[custPulseGenInput.value]
            print('Successfully loaded from locals!')
        
        elif custPulseGenInput.value in globals():
            func = globals()[custPulseGenInput.value]
            print('Successfully loaded from globals!')
            
        else:
            func = None
            #with paramOutput:
            print("Unable to find '{}'!".format(custPulseGenInput.value))
            #warnings.warn('Warning! Variable: {} not found!'.format(paramVar.value))
        
        return func
    
    
    for pInd, prot in enumerate(protList):
        pSet = protParams[prot]
        if prot is 'custom': ### Exception here for loading pulse generator function
            custPulseGenInput = widgets.Text(placeholder='<Pulse generator function>', description='$\phi(t)$') # value=str(pSet[param].value)
            protHeaders[pInd] = widgets.HBox(children = [custPulseGenInput, widgets.HTML(value=protParamNotes[prot]['phi_ft'])])
        else:
            protHeaders[pInd] = widgets.HTML(value='')
        for i, param in enumerate(pSet):
            label = '$' + protParamLabels[param] + '$'
            if isinstance(pSet[param].value, list): # or isinstance(pSet[param].value, str):
                prot_pValArr[pInd][i] = widgets.Text(value=str(pSet[param].value), description=label) # np.asarray #param
                spacerProtParams = spacerNull
            elif isinstance(pSet[param].value, bool):
                prot_pValArr[pInd][i] = widgets.Dropdown(options=boolDict, value=pSet[param].value, description=label) #param
                spacerProtParams = spacerProtDD
            
            #elif pSet[param].value is None: # Used for phi_ft
            #    if pSet[param].name is "phi_ft":
            #        prot_pValArr[pInd][i] = widgets.Text(value=str(pSet[param].value), description=label)
            #        spacerProtParams = spacerNull
            
            else:
                if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf):
                    prot_pValArr[pInd][i] = widgets.FloatText(value=pSet[param].value, description=label)
                else:
                    prot_pValArr[pInd][i] = widgets.BoundedFloatText(value=pSet[param].value, min=pSet[param].min, max=pSet[param].max, description=label)
                ###prot_pValArr[pInd][i].width = '150px'
                spacerProtParams = spacerProtNum
            
            
            if param in protUnitLabels: #pSet[param].expr is not None:
                prot_unitArr[pInd][i] = widgets.Dropdown(options=[protUnitLabels[param]], value=protUnitLabels[param])#pSet[param].expr ### Change units handling
            else:
                prot_unitArr[pInd][i] = widgets.HTML(value='')
            
            
            #prot_unitArr[pInd][i] = widgets.Latex(value='$[\mathrm{'+protUnitLabels[param]+'}]$')
            protParamNotesWid[pInd][i] = widgets.HTML(value=protParamNotes[prot][param]) # Latex is ok too
            prot_pBoxArr[pInd][i] = widgets.HBox(children=[prot_pValArr[pInd][i], spacerProtParams, prot_unitArr[pInd][i], protParamNotesWid[pInd][i]])
            
            """
            if param in protUnitLabels: #pSet[param].expr is not None:
                print(param, '-->', protUnitLabels[param])
                prot_unitArr[pInd][i] = widgets.Dropdown(options=[protUnitLabels[param]], value=protUnitLabels[param])#pSet[param].expr ### Change units handling
                prot_pBoxArr[pInd][i] = widgets.HBox(children=[prot_pValArr[pInd][i], prot_unitArr[pInd][i]])
            else:
                prot_pBoxArr[pInd][i] = widgets.HBox(children=[prot_pValArr[pInd][i]])
            """
            
            #i+=1
        
        """
        protStimHTML[pInd] = widgets.HTML()
        ### Add stimulus figures by printing to file
        # IPython.core.pylabtools.print_figure(fig, fmt='png', bbox_inches='tight', **kwargs)
        stimFig = plt.figure()
        x = np.linspace(0, 3*np.pi, 500)
        #plt.figure(Ifig.number)
        plt.plot(x, np.sin(x**2))
        #fdata64 = base64.b64encode(print_figure(stimFig))
        title = 'Stimulus for {}'.format(prot)
        ###print(title,fdata64)
        html_fig = ''#'<img alt="{}" src="data:image/png;base64,{}">'.format(title,fdata64)
        #return html_tpl.format(**locals())
        plt.close(stimFig)
        protStimHTML[pInd].value = html_fig
        """
        
        protFigHTML[pInd] = widgets.HTML()
        exampleProt = '{}{}6s.{}'.format(guiFigDir,prot,'png')#saveFigFormat) # fDir
        if os.path.isfile(exampleProt):
            protFigHTML[pInd].value='<img src="{}" alt=Example {} width=200px>'.format(exampleProt,prot)
        else:
            protFigHTML[pInd].value='Example Figure'
        protParamBoxes[pInd] = widgets.Box(children=prot_pBoxArr[pInd])
        
        #protParamNotesBoxes[pInd] = widgets.VBox(children=protParamNotesWid[pInd])
        
        '''#protNotesBoxes[pInd] = widgets.HBox(children=[protStimHTML[pInd], protFigHTML[pInd]])#[figHTML[prot],eqBox[prot]])'''
        #protNotesBoxes[pInd] = widgets.HBox(children=[ protParamNotesBoxes[pInd], protFigHTML[pInd] ])#[figHTML[prot],eqBox[prot]])
        protNotesBoxes[pInd] = widgets.HBox(children=[ protFigHTML[pInd] ])#[figHTML[prot],eqBox[prot]])
        
        
        protBoxes[pInd] = widgets.HBox(children=[ protParamBoxes[pInd], protNotesBoxes[pInd] ])#modelBox
        #protBoxes[pInd].margin = '5px'
        #display(protBoxes[pInd])
        
        protBoxesWrapper[pInd] = widgets.VBox(children=[protHeaders[pInd], protBoxes[pInd]])
        
        
    
    ##### Protocol parameters tab #####
    def onChangeProtTab(change): #(name,value):
        protDropdown.value = protList[change['new']] #[value] #protDropdown.value_name = protParamsTabs._titles[value] #[protParamsTabs.selected_index]
    protParamsTabs = widgets.Tab(description='Parameter Settings', children=protBoxesWrapper) # protBoxes)# \
    ###protParamsTabs.margin = '5px'
    #protParamsTabs.on_trait_change(onChangeProtTab, 'selected_index')
    protParamsTabs.observe(onChangeProtTab, names='selected_index')
    protBox = widgets.VBox(children=[protParamsHead, protParamsTabs]) #HBox
    
    ##### Configure tabs for abstraction layers #####
    paramTabs = widgets.Tab(description='Parameter Settings', children=[fitBox, modelParamsTabs, protBox, simParamsTabs]) # protParamsTabs #fitParamsTabs #,values=['Model', 'Protocol']) #E_box,k_box
    #####display(paramTabs) # Commented to nest in GUI box
    paramTabs.selected_index = TabGroups['Models'] # Set to show model parameters initially
    paramTabs.visible = False
    ###paramTabs.margin = '5px'
    
    
    
    
    
    # ### Create NEURON Control Bar
    # NEURONbox = widgets.HBox(children=[hocFile, hocLoad])
    # display(NEURONbox)
    # NEURONbox.visible=False
    # NEURONbox.align = 'center'
    
    # hocLoad.button_style = 'warning' #add_class('btn-warning')
    
    
    # ### Create Brian Control Bar
    # Brianbox = widgets.HBox(children=[brianFile, brianLoad])
    # display(Brianbox)
    # Brianbox.visible=False
    # Brianbox.align = 'center'
    
    # brianLoad.button_style = 'warning' #add_class('btn-warning')
    
    
    
    ### Container box for the entire GUI
    GUI = widgets.VBox()
    GUI.children = [runBar, paramsControlBar, paramTabs] #[fitBar,runBar,paramsControlBar,paramTabs,NEURONbox,Brianbox]
    display(GUI)
    
    
    
    ##### Add formatting here after display #####
    
    ### align := ['start', 'center', 'end', 'baseline', 'stretch']
    ### border_color := 'Valid HTML colour'
    ### border_radius := Specify curvature in pixels
    ### border_style := ['none', 'hidden', 'dotted', 'dashed', 'solid', 'double', 'groove', 'ridge', 'inset', 'outset', 'initial', 'inherit', '']
    
    ### Set Button styles after displaying the runBar
    paramsButton.button_style = 'info'
    #stateButtons.button_style = 'success' ### Removed due to widgets bug
    stateButtons.margin = '5px'
    simDropdown.button_style = 'success'
    simDropdown.margin = '5px'
    protDropdown.button_style = 'success'
    #protDropdown.margin = '5px'
    #fitButton.button_style = 'success'
    runButton.button_style = 'danger'
    runButton.margin = '5px'
    runBar.align = 'center'
    verboseSlide.width = '48px' # '60px'
    
    ### Set paramTabs formatting
    
    ### Set Parameters Bar style
    paramVar.width = '150px' #set_css({'width': '150px'})
    paramBox.width = '250px'
    paramOutput.width = '300px'
    paramLoad.button_style = 'warning' #add_class('btn-warning')
    paramsResetButton.button_style = 'warning' #add_class('btn-warning')
    paramsControlBar.align = 'center'
    
    #lamVal.width = '75px'
    #lamUnits.width = '75px'
    
    
    ### Set fitBar formatting
    dataVar.width = '150px' #set_css({'width': '150px'})
    #dataVar.margin = '0px'
    dataBox.width = '150px'
    dataLoad.button_style = 'warning' #add_class('btn-warning')
    # statesToFitButtons.button_style = 'info' ### Removed due to widgets bug
    runFitButton.button_style = 'danger' #add_class('btn-danger')
    #runFitButton.margin = '5px'
    characteriseButton.button_style = 'success'
    # Set Fit Bar formatting after display
    dataLoad._dom_classes = ('margin-left','10px')
    #fitBar.align = 'center'
    exportFitButton.button_style = 'info'
    runSSAcheck.margin = '5px'
    #plotExpData.margin = '5px'
    runPostOpt.margin = '5px'
    runFitButton.margin = '5px'
    #characteriseButton.margin = '5px'
    exportFitButton.margin = '5px'
    fitMethods.margin = '5px'
    postOptFitMethods.margin = '5px'
    fitBarMain.align = 'center'
    fitBarMeta.align = 'center'
    fitBar.align = 'start'
    relaxFactWid.width = '100px'
    relaxFactWid.margin = '5px'
    fitParamsTabs.width = '800px'
    fitBox.margin = '5px'
    
    
    ### Set Tab titles ### # Titles must be set after display of tabs
    paramTabs.set_title(0, 'Fitting Parameters')
    paramTabs.set_title(1, 'Model Parameters')
    paramTabs.set_title(2, 'Protocol Parameters')
    paramTabs.set_title(3, 'Simulator Parameters')
    
    stateButtons.tooltips = modelTitles
    statesToFitButtons.tooltips = modelTitles
            
    for ind, title in enumerate(modelTitles):
        modelParamsTabs.set_title(ind, title)
        fitParamsTabs.set_title(ind, title)
        
    # modelParamsTabs.set_title(0, 'Three-state model')
    # modelParamsTabs.set_title(1, 'Four-state model')
    # modelParamsTabs.set_title(2, 'Six-state model')
    modelParamsTabs.width = '800px' #set_css({'width': '800px'}) # 800
    
    ### Set Parameter Tabs style
    for pInd, prot in enumerate(protList): # Change to protParams #for prot, pSet in protParams.items():
        protParamsTabs.set_title(pInd, prot)
        protNotesBoxes[pInd].flex = True #protNotesBoxes[pInd].add_class('box-flex1')
        protFigHTML[pInd].width = '200px' #set_css({'width': '200px', 'margin-left': '20px'})
        protFigHTML[pInd].margin = '20px'
        # protParamBoxes[pInd].width = '200px' # Does not work
        for i, param in enumerate(protParams[prot]):
            if isinstance(protParams[prot][param].value, list):
                prot_pValArr[pInd][i].width = '300px' #'150px'#'300px' # Hack since surrounding space does not contract '150px' #set_css({'width': '150px'})
            elif isinstance(protParams[prot][param].value, bool):
                prot_pValArr[pInd][i].width = '278px' 
            else:
                prot_pValArr[pInd][i].width = '300px' #'150px'
                if param in protUnitLabels:
                    prot_unitArr[pInd][i].description = ' ' # Hack to correct spacing
            protParamNotesWid[pInd][i].margin = '5px'
        protBoxes[pInd].margin = '5px'
    
    ### Hacks to fix the spacing irregularities betweeen e.g. Text and Float widgets by adjusting HTML widget spacers
    spacerProtNum.width = '128px'
    spacerProtDD.width = '124px'
    spacerNull.width = '0px'
    
    
    ### Set Simulator Tabs style
    for sInd, sim in enumerate(simList):
        simParamsTabs.set_title(sInd, sim)
        simNotesBoxes[sInd].flex = True #protNotesBoxes[pInd].add_class('box-flex1')
        simFigHTML[sInd].width = '200px' #set_css({'width': '200px', 'margin-left': '20px'})
        simFigHTML[sInd].margin = '20px'
        for i, param in enumerate(simParams[sim]):
        #for i in range(len(simParams[sim])): # pSet.keys(): #for i in range(len(protParams[prot])):
            sim_pValArr[sInd][i].width = '150px' #set_css({'width': '150px'})
            if param in simUnitLabels:
                sim_unitArr[sInd][i].description = ' ' # Hack to correct spacing
        simBoxes[sInd].margin = '5px'
    
    
    ### Hack to tile parameter fields horizontally - must come after displaying the parent
    #for model, pSet in modelParams.items():
    for model in modelParams:
        m = modelList.index(model)
        eqBox[m].width = '150px'
        eqBox[m].margin = '10px'
        modelNotesBoxes[m].width = '300px'
        modelNotesBoxes[m].margin = '20px'
        modelBoxes[m].margin = '5px'
        #for i, param in enumerate(pSet):
        #    pValArr[m][i].width = '300px'
    
    spacer.width = '150px' # Equal to the width of a drop-drown menu
    
    for model, pSet in modelParams.items(): #range(len(modelParams)):
        #pSet = modelParams[m]            # Set of n-state model parameters
        m = modelList.index(model)
        #nParams = len(pSet)#.keys())
        i=0
        for key in pSet:#.keys(): #, value in pSet.items():
            pfLabArr[m][i].width = '50px'#'20px'
            #pfValArr[m][i].width = '150px'
            fExprArr[m][i].width = '150px'
            fVaryArr[m][i].width = '30px'
            i+=1
        modelFitBoxes[m].margin = '5px'
        
        # Loop over Model Parameters Tab
        #for i, param in enumerate(pSet):
            #pValArr[m][i].width = '150px'
            # #if param in protUnitLabels:
            #unitArr[m][i].width = '150px'
            #unitArr[m][i].description = ' ' # Hack to correct spacing
    
    modelParamsTabs.margin = '5px'
    simParamsTabs.margin = '5px'
    protParamsTabs.margin = '5px'
    paramTabs.margin = '5px'
    
    
    return #GUI
    
#GUI = widgets.Box()
#interact(runModel, nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'delta', 'rectifier', 'shortPulse', 'recovery'), saveData=True, verbose=1);

if __name__ == '__main__':
    loadGUI()
    
