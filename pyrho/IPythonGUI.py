from IPython.html.widgets import interact
from IPython.html import widgets
from IPython.display import display
from IPython.display import clear_output
from IPython.utils.traitlets import link
from IPython.core.pylabtools import print_figure
import base64

from .models import *
from .simulators import *
from .protocols import *
from .loadData import *
import warnings
import time
from collections import OrderedDict
import os.path
import ast


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

#%run -i modProtocols.py
#%run -i models.py

# State keys must be padded with a leading ' ' to avoid a widgets bug: https://github.com/ipython/ipython/issues/6469
#statesDict = OrderedDict([(' '+s,i) for i,s in enumerate(list(modelParams.keys()))]) # enumerate(modelList)
#statesDict = OrderedDict([(' 3',0), (' 4',1), (' 6',2)]) 
statesDict = OrderedDict([(s,i) for i,s in enumerate(list(modelParams))]) #.keys()
statesArray = list(statesDict) #.keys() #[' 3', ' 4', ' 6'] # [u' 3',u' 4',u' 6']

TabGroups = {'Fit':0, 'Models':1, 'Simulators':2, 'Protocols':3}
#TabGroups = {'Models':0, 'Simulators':1, 'Protocols':2}

#clearDelay = 1.5 # Pause [s] before clearing text entry fields

    
    


    
    

def loadGUI():
    
    
    ##### Model fitting bar Functions #####

    def fitToggle(name, value):
        if value == True:
            fitBar.visible = True
            #time.sleep(clearDelay) # Pause then clear the input field
            #dataVar.value = ''
        else:
            fitBar.visible = False
            #dataVar.value='<variable name>'
        return
        
    def onLoad(name):
        global dataSet
        print('Loading: "',dataVar.value, '"...', end=' ')
        if dataVar.value in vars(): ### locals()? # http://stackoverflow.com/questions/7969949/whats-the-difference-between-globals-locals-and-vars
            dataSet = vars()[dataVar.value]
            useFitCheck.value = True
            print('Successfully loaded from vars!')
        elif dataVar.value in globals():
            dataSet = globals()[dataVar.value] # eval(dataVar.value) ##### Is this safe?
            useFitCheck.value = True
            print('Successfully loaded from globals!')
        else:
            dataSet = None
            useFitCheck.value = False
            warnings.warn('Warning! Variable: {} not found!'.format(dataVar.value))
        return dataSet
        
    def runFitButton_on_click(b): # on_button_clicked
        #global fitParamsPopup
        global clearOutput
        if clearOutput.value:
            clear_output()
            #if 'fitParamsPopup' in vars() or 'fitParamsPopup' in globals():
            #    fitParamsPopup.close()
        global fitRhO
        fitRhO = fitModels(dataSet, int(statesToFitButtons.value))
        #fitParamReport = widgets.TextareaWidget(description='Report:',value=fitRhO.reportParams())
        #fitParamsPopup = widgets.PopupWidget(children=[fitParamReport],button_text='Fitted Parameters',description='Fitted {} state model parameters from: {}'.format(int(statesToFitButtons.value),dataVar.value))
        #display(fitParamsPopup)
        if runSSAcheck.value == True:
            characterise(fitRhO)    
        return
    
    
    
    ##### Run Bar Functions #####

    def protDropdownChange(name,value):
        paramTabs.selected_index = TabGroups['Protocols'] # Set to protocol parameters
        protParamsTabs.selected_index = protList.index(value) # Get the index of the selected protocol #protIndDict[value]


    def runButton_on_click(b): # on_button_clicked
        if clearOutput.value:
            clear_output()
        runModel(stateButtons.value, protDropdown.value, simDropdown.value, saveButton.value, verboseSlide.value) #int(stateButtons.value)
        return        
        
    def changeModel(name,value):
        paramTabs.selected_index = TabGroups['Models'] # Set to model parameters
        modelParamsTabs.selected_index = statesDict[value]    

    def paramsToggle(name, value):
        
        if value == True: # Set to current model and protocol tabs
            paramsControlBar.visible = True
            paramTabs.visible = True
            modelParamsTabs.selected_index = statesDict[stateButtons.value]
            protParamsTabs.selected_index = protList.index(protDropdown.value) #protIndDict[protDropdown.value]
        else:
            paramsControlBar.visible = False
            paramTabs.visible = False
        return

    def simDropdownChange(name,value):
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
        simParamsTabs.selected_index = simList.index(value) # Get the index of the selected protocol #protIndDict[value]
        return
        
        
    ##### NEURON bar functions #####
    def onHocLoad(name):
        print('Loading: "',hocFile.value, '"...', end=' ')
        try: # Load mechanism and set some appropriate parameters
            h = hoc.HocObject()
            h.xopen(hocFile.value)
            h.load_mechanisms() #h('nrn_load_dll("libnrnmech.so")')
            
            #from neuron import gui # Load NEURON GUI
            
        except:
            print('Error! File: {} not found!'.format(dataVar.value))
        return

        
    ##### Brian bar functions #####
    def onBrianLoad(name):
        print('Loading: "',brianFile.value, '"...', end=' ')
        try: # Load mechanism and set some appropriate parameters
            print('Finish me!!!')
            
        except:
            print('Error! File: {} not found!'.format(dataVar.value))
        return
        

        
    ##### Parameter Bar and Tabs functions #####
    
    def paramSetOnLoad(name):
        global paramSet
        print('Loading: "',paramVar.value, '"...', end=' ')
        
        if paramVar.value in vars():
            paramSet = vars()[paramVar.value]
            print('Successfully loaded from vars!')
            
        elif paramVar.value in globals():
            paramSet = globals()[paramVar.value]
            print('Successfully loaded from globals!')
            
        else:
            paramSet = None
            warnings.warn('Warning! Variable: {} not found!'.format(paramVar.value))
            
        if paramTabs.selected_index == TabGroups['Models']:
            mInd = modelParamsTabs.selected_index
            model = statesArray[mInd]
            #mInd = statesDict[' '+str(nStates)]
            #pSet = modelParams[mInd] # Default model parameters
            setGUIparams(modelParams[model],pValArr[mInd][:]) # setGUImodelParams(paramSet,model)
            
        elif paramTabs.selected_index == TabGroups['Simulators']:
            sInd = simParamsTabs.selected_index
            simulator = simList[sInd]
            setGUIparams(simParams[simulator],sim_pValArr[sInd][:])
            
        elif paramTabs.selected_index == TabGroups['Protocols']: # Protocols tab selected
            pInd = protParamsTabs.selected_index
            #pInd = protIndDict[protocol]
            protocol = protList[pInd]
            #pSet = protParams[protocol] # Default protocol parameters
            setGUIparams(protParams[protocol],prot_pValArr[pInd][:]) #setGUIprotParams(paramSet,protocol)
            
        else: 
            raise ValueError('Unknown Tab Index!')
            
        return paramSet

        
    def resetButton_on_click(b): # on_button_clicked
        
        if paramTabs.selected_index == TabGroups['Models']:
            mInd = modelParamsTabs.selected_index
            model = statesArray[mInd]
            #mInd = statesDict[' '+str(nStates)]
            #pSet = modelParams[model] # Default model parameters
            #pSet = modelParams[statesArray[mInd]]
            setGUIparams(modelParams[model],pValArr[mInd][:]) #setGUImodelParams(pSet,model)
            
        elif paramTabs.selected_index == TabGroups['Simulators']:
            sInd = simParamsTabs.selected_index
            simulator = simList[sInd]
            setGUIparams(simParams[simulator],sim_pValArr[sInd][:])
            
        elif paramTabs.selected_index == TabGroups['Protocols']:
            pInd = protParamsTabs.selected_index
            #pInd = protIndDict[protocol]
            protocol = protList[pInd]
            #pSet = protParams[protocol] # Default protocol parameters
            setGUIparams(protParams[protocol],prot_pValArr[pInd][:]) #setGUIprotParams(pSet,protocol)
            
        else: # Protocols tab selected
            raise ValueError('Unknown Tab Index!')
            
        return
    
    def modelTabChange(name,value):
        stateButtons.value = statesArray[value] #' '+str(value)
    
    def simTabChange(name,value):
        simDropdown.value = simList[value] #simList.index(value)# simDropdown.selected_label #simParamsTabs._titles[value] #[protParamsTabs.selected_index]
        
    def protTabChange(name,value):
        protDropdown.value = protList[value] #protDropdown.value_name = protParamsTabs._titles[value] #[protParamsTabs.selected_index]
    
    
    ### Utility functions 
    # def getGUImodelParams(model):
        # userParams = Parameters()
        # mInd = statesDict[model] #statesDict[' '+str(nStates)]
        # pSet = modelParams[modelList[mInd]]           
        # i=0
        # for key in pSet:  # Set of n-state model pararmeters
            # userParams.add(key, value=pValArr[mInd][i].value)
            # i+=1
        # return userParams

    # def setGUImodelParams(pSet,model):
        # mInd = statesDict[model] #statesDict[' '+str(nStates)]
        # i=0
        # for key in pSet: #, value in pSet.items():
            # pValArr[mInd][i].value = pSet[key].value
            # i+=1
        # return
    
    
    # def getGUIsimParams(simulator):
        # userParams = Parameters()
        # sInd = simList.index(simulator) #simIndDict[simulator]
        # pSet = simParams[simulator]
        # #pluginInd = paramGroup.keys().index(plugin)
        # i=0
        # for param in pSet: #.keys():
            # if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                # userParams.add(param, value=ast.literal_eval(sim_pValArr[sInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            # else:
                # userParams.add(param, value=sim_pValArr[sInd][i].value)
            # i+=1
        # return userParams
        
    # def setGUIsimParams(pSet,simulator):
        # sInd = simList.index(simulator) #simIndDict[simulator]
        # i=0
        # for param in pSet: #.keys(): #, value in pSet.items():
            # if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                # sim_pValArr[sInd][i].value = str(pSet[param].value)
            # else:
                # sim_pValArr[sInd][i].value = pSet[param].value
            # i+=1
        # return
        
    
    # def getGUIprotParams(protocol):
        # userParams = Parameters()
        # pInd = protList.index(protocol) #protIndDict[protocol]
        # pSet = protParams[protocol]
        # i=0
        # for param in pSet: #.keys():
            # if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                # userParams.add(param, value=ast.literal_eval(prot_pValArr[pInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            # else:
                # userParams.add(param, value=prot_pValArr[pInd][i].value)
            # i+=1
        # return userParams
        
    # def setGUIprotParams(pSet,protocol):
        # pInd = protList.index(protocol) #protIndDict[protocol]
        # i=0
        # for param in pSet: #.keys(): #, value in pSet.items():
            # if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                # prot_pValArr[pInd][i].value = str(pSet[param].value)
            # else:
                # prot_pValArr[pInd][i].value = pSet[param].value
            # i+=1
        # return
    
    
    
    #getGUIparams(modelParams[model],pValArr[modelParams.keys().index(model)][:])
    #getGUIparams(simParams[simulator],sim_pValArr[simParams.keys().index(simulator)][:])
    #getGUIparams(protParams[protocol],prot_pValArr[protParams.keys().index(protocol)][:])
    
    # pSet = protParams[protocol]
    # widgetList = prot_pValArr[pInd][:]
    def getGUIparams(pSet, widgetList):
        userParams = Parameters()
        for i, param in enumerate(pSet):
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                userParams.add(param, value=ast.literal_eval(widgetList[i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            else:
                userParams.add(param, value=widgetList[i].value)
        return userParams
    
    def setGUIparams(pSet, widgetList):
        for i, param in enumerate(pSet):
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                widgetList[i].value = str(pSet[param].value)
            else:
                widgetList[i].value = pSet[param].value
    
    
    
        
    ##### Main Run function #####
    
        #saveData = True
        #@interact(nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'saturate', 'inwardRect', 'varyPL', 'varyIPI'), saveData=True, verbose=1)

    def runModel(model, protocol, simulator='Python', saveData=True, verbose=1): #verboseSlide.value
        
        verbose = verbose
        nStates = int(model)
        print("\nRunning Protocol '{}' on the {}-state model...".format(protocol,nStates))
        print('--------------------------------------------------------------------------------\n')
        
        ### Choose between 3, 4, & 6 state models
        if useFitCheck.value: # fitButton.value and 
            RhO = fitRhO
            
            ### To set GUI parameters, widget arrays should be dictionaries so subsets of parameters can be set in arbitrary order!!!
            #setGUIprotParams(expProtParams,'custom')
            
        else:
            RhO = selectModel(int(model))
            if verbose > 0:
                print(RhO)
            userModelParams = Parameters()
            mInd = statesDict[model] #statesDict[' '+str(nStates)]
            pSet = modelParams[modelList[mInd]]
            userModelParams = getGUIparams(pSet,pValArr[mInd][:])
            # nParams = len(pSet)#.keys())
            # i=0
            # for key in pSet:#.keys():  # Set of n-state model pararmeters
                # userModelParams.add(key,value=pValArr[mInd][i].value)
                # i+=1
            RhO.setParams(userModelParams)
        
        
        
        ### Get Simulator Parameters
        if verbose > 0:
            print('Simulating with {}...'.format(simulator))
        #if simulator != 'Python': ########################################### Hack!!!!!
        userSimParams = Parameters()
        sInd = simList.index(simulator)
        userSimParams = getGUIparams(simParams[simulator],sim_pValArr[sInd][:])
        Sim = simulators[simulator](RhO, userSimParams)
        
        
        
        ### Get Protocol Parameters
        userProtParams = Parameters()
        pInd = protList.index(protocol) #protIndDict[protocol]
        userProtParams = getGUIparams(protParams[protocol],prot_pValArr[pInd][:])
        # pSet = protParams[protocol]
        # i=0
        # for param in pSet:#.keys(): #, value in pSet.items():
            # if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                # userProtParams.add(param, value=ast.literal_eval(prot_pValArr[pInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            # else:
                # userProtParams.add(param, value=prot_pValArr[pInd][i].value)
            # i+=1
        Prot = protocols[protocol](userProtParams)
        Prot.runProtocol(Sim,RhO,verbose)
        if verbose > 0: #saveData:
            Prot.plotProtocol(Sim,RhO,verbose)
        
        print("\nFinished!")
        print('================================================================================\n\n')
        
        return #Prot, RhO
        
        
        

        
    
    ##### Fit Bar #####
    
    # ### Create Fit Bar toggle
    # fitButton = widgets.ToggleButton(description='Fit data', value=False)
    # fitButton.on_trait_change(fitToggle,'value') #paramsButton.get_state()

    # ### Create Data set entry
    # dataVar = widgets.Text(description='Data Set: ',placeholder='<variable name>')
    # dataLoad = widgets.Button(description='Load')
    # dataLoad.on_click(onLoad)
    
    # ### Create Fit Model States buttons
    # statesToFitButtons = widgets.ToggleButtons(description='Highest model to fit: ',options=statesArray)#,value=u' 3') #https://github.com/ipython/ipython/issues/6469
    # #fit3sCheck = widgets.Checkbox(description='Fit 3 state model', value=True)
    # #fit4sCheck = widgets.Checkbox(description='Fit 4 state model', value=False)
    # #fit6sCheck = widgets.Checkbox(description='Fit 6 state model', value=False)
    
    # ### Create Checkboxes
    # runSSAcheck = widgets.Checkbox(description='Characterise', value=True)
    # useFitCheck = widgets.Checkbox(description='Use fit', value=False)
    # plotExpData = widgets.Checkbox(description='Plot data: ', value=True)
    
    # ### Create Run Button
    # runFitButton = widgets.Button(description="Fit!")
    # runFitButton.on_click(runFitButton_on_click)
    
    # ### Create Fit Bar
    # fitBar = widgets.HBox(children=[dataVar, dataLoad, statesToFitButtons, runSSAcheck, useFitCheck, plotExpData, runFitButton]) #fit3sCheck, fit4sCheck, fit6sCheck
    # display(fitBar)
    # fitBar.visible=False
    
    # ### Set formatting
    # dataVar.width = '150px' #set_css({'width': '150px'})
    # dataLoad.button_style = 'warning' #add_class('btn-warning')
    # runFitButton.button_style = 'danger' #add_class('btn-danger')
    # # Set Fit Bar formating after display
    # fitBar.align = 'center'
    
    
    ##### Run bar #####
    
    ### Protocol Dropdown
    protDropdown = widgets.Dropdown(options=protList,value='saturate')###########,value='saturate') #protDict
    protDropdown.on_trait_change(protDropdownChange,'value')
    
    ### Run Button
    runButton = widgets.Button(description="Run!")
    
    runButton.on_click(runButton_on_click)
    
    ### Model states button
    
    stateButtons = widgets.ToggleButtons(description='Model states: ',options=statesArray,) #https://github.com/ipython/ipython/issues/6469
    stateButtons.on_trait_change(changeModel,'value')#stateButtons.value_name  #'_view_name'
    
    ### Create states buttons link
    #modelParamsLink = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    ### Parameters Toggle
    paramsButton = widgets.ToggleButton(description='Parameters', value=False)
    paramsButton.on_trait_change(paramsToggle,'value')
    
    # ### Load NEURON parameters
    # hocFile = widgets.Text(description='HOC file: ',placeholder='<file name>')
    # hocLoad = widgets.Button(description='Load')
    # hocLoad.on_click(onHocLoad)
    
    # ### Load Brian parameters
    # brianFile = widgets.Text(description='Brian file: ',placeholder='<file name>')
    # brianLoad = widgets.Button(description='Load')
    # brianLoad.on_click(onBrianLoad)
    
    ### Output format Dropdown
    simDropdown = widgets.Dropdown(options=['Python', 'NEURON', 'Brian'],value='Python')
    simDropdown.on_trait_change(simDropdownChange,'value')
    
    
    
    
    ### Plot and Save button
    saveButton = widgets.Checkbox(description='Plot & Save', value=True)
    saveButton.visible = False

    ### Verbosity slider
    #verboseSlide = widgets.FloatProgressWidget(value=1, min=-1, max=3, step=1, description='Verbosity:')
    verboseSlide = widgets.IntSlider(value=1, min=0, max=3, description='Output:')
    verboseSlide.visible=True #False

    ### Clear ouput checkbox
    clearOutput = widgets.Checkbox(description='Clear', value=True)
    
    
    #### Run Bar container
    runBar = widgets.HBox(children=[stateButtons,simDropdown,protDropdown,paramsButton,saveButton,verboseSlide,clearOutput,runButton]) #fitButton,
    display(runBar)

    
    
    ### Set Button styles after displaying the runBar
    paramsButton.button_style = 'info'
    stateButtons.button_style = 'success'
    stateButtons.margin = '5px'
    simDropdown.button_style = 'success'
    simDropdown.margin = '5px'
    protDropdown.button_style = 'success'
    protDropdown.margin = '5px'
    #fitButton.button_style = 'success'
    runButton.button_style = 'danger'
    runButton.margin = '5px'
    
    #stateButtons.add_class('btn-primary')
    #protDropdown.add_class('btn-success')
    #simDropdown.add_class('btn-warning')
    
    ### button bootstrap styles: default, primary, success, info, warning, danger, link
    

    runBar.align = 'center'
    verboseSlide.width = '60px'
    
    ##### Parameters control bar #####
    
    paramVar = widgets.Text(description='Parameter Set: ',placeholder='<variable name>')
    paramLoad = widgets.Button(description='Load')
    paramLoad.on_click(paramSetOnLoad)
    
    
    ### Reset Parameters Button
    paramsResetButton = widgets.Button(description="Reset")
    paramsResetButton.on_click(resetButton_on_click)
    
    ### Create Parameters control bar
    paramsControlBar = widgets.HBox(children=[paramVar,paramLoad,paramsResetButton])
    display(paramsControlBar)
    paramsControlBar.visible = False
    
    
    ### Set Parameters Bar style
    paramVar.width = '150px' #set_css({'width': '150px'})
    paramLoad.button_style = 'warning' #add_class('btn-warning')
    paramsResetButton.button_style = 'warning' #add_class('btn-warning')
    
    paramsControlBar.align = 'center'
    
    
    
    
    ##### Create Parameters Tab widgets #####
    
    ### lmfit Parameters() are ordereddict objects (https://github.com/lmfit/lmfit-py/blob/master/lmfit/parameter.py) but it might be better to change pValArr and prot_pValArr etc. to ordereddicts too...
    
    
    boolDict = OrderedDict([('True',True), ('False',False)])
    
    
    ##### Fit Data parameters #####
    
    ### Create Data set entry
    dataVar = widgets.Text(description='Data Set: ',placeholder='<variable name>')
    dataLoad = widgets.Button(description='Load')
    dataLoad.on_click(onLoad)
    
    ### Create Fit Model States buttons
    #statesToFitButtons = widgets.ToggleButtons(description='Highest model to fit: ',options=statesArray)#,value=u' 3') #https://github.com/ipython/ipython/issues/6469
    fitLabel = widgets.HTML(value='Models to fit: ')
    fit3sCheck = widgets.Checkbox(description='3 state', value=True)
    fit4sCheck = widgets.Checkbox(description='4 state', value=False)
    fit6sCheck = widgets.Checkbox(description='6 state', value=False)
    
    ### Hyperparameters
    p0fVvar = widgets.Text(value=str(p0fV), description='p0fV: ')
    p0IPIvar = widgets.Text(value=str(p0IPI), description='p0IPI: ')
    
    
    ### Create Checkboxes
    runSSAcheck = widgets.Checkbox(description='Characterise', value=True)
    useFitCheck = widgets.Checkbox(description='Use fit', value=False)
    plotExpData = widgets.Checkbox(description='Plot data: ', value=True)
    
    ### Create Run Button
    runFitButton = widgets.Button(description="Fit!")
    runFitButton.on_click(runFitButton_on_click)
    
    ### Create Fit Bar
    #fitBar = widgets.HBox(children=[dataVar, dataLoad, statesToFitButtons, runSSAcheck, useFitCheck, plotExpData, runFitButton]) #fit3sCheck, fit4sCheck, fit6sCheck
    #display(fitBar)
    #fitBar.visible=False
    
    ### Set formatting
    dataVar.width = '150px' #set_css({'width': '150px'})
    #dataVar.margin = '0px'
    dataLoad.button_style = 'warning' #add_class('btn-warning')
    runFitButton.button_style = 'danger' #add_class('btn-danger')
    runFitButton.margin = '5px'
    # Set Fit Bar formating after display
    dataLoad._dom_classes = ('margin-left','10px')
    #fitBar.align = 'center'
    
    dataLoadBar = widgets.HBox(children=[dataVar, dataLoad, useFitCheck], align='center')
    statesBar = widgets.HBox(children=[fitLabel, fit3sCheck, fit4sCheck, fit6sCheck], align='center') #.align = 'center'
    hyperParams = widgets.VBox(children=[p0fVvar, p0IPIvar])
    fitOutputBar = widgets.HBox(children=[runSSAcheck, plotExpData, runFitButton], align='center')
    fitParamBoxes = widgets.VBox(children=[dataLoadBar, statesBar, hyperParams, fitOutputBar]) # Left side container
    fitNotesBoxes = widgets.VBox(children=[]) # Right side container for figure and equations
    fitBox = widgets.HBox(children=[fitParamBoxes])
    fitBox.margin = '5px'
    fitParamsTabs = widgets.HBox(children=[fitBox]) # fitNotesBoxes # Hack to have a parent tab with no children
    
    
    
    
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
                if pSet[param].expr == None: # No units
                    pluginParams[pluginInd][pInd] = widgets.HBox(children=[pluginValues[pluginInd][pInd]])
                else:
                    pluginUnits[pluginInd][pInd] = widgets.Dropdown(options=[pSet[param].expr],value=pSet[param].expr)
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
        
        #pluginParamsTabs.on_trait_change(simTabChange,'selected_index') # External
    
    
    
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
        #pSet = modelParams[m]            # Set of n-state model pararmeters
        m = modelList.index(model)
        nParams = len(pSet)#.keys())
        
        i=0
        for key in pSet:#.keys(): #, value in pSet.items():
            if pSet[key].min == None or pSet[key].max == None:
                pValArr[m][i] = widgets.FloatText(value=pSet[key].value, description=key)#"{} [{}]:".format(key, pSet[key].expr))
            else:
                pValArr[m][i] = widgets.BoundedFloatText(value=pSet[key].value, min=pSet[key].min, max=pSet[key].max, description=key)#"{} [{}]:".format(key, pSet[key].expr))
            if not pSet[key].expr == None:
                unitArr[m][i] = widgets.Dropdown(options=[pSet[key].expr],value=pSet[key].expr)
                pBoxArr[m][i] = widgets.HBox(children=[pValArr[m][i],unitArr[m][i]])
            else:
                pBoxArr[m][i] = widgets.HBox(children=[pValArr[m][i]])
            i+=1
        
        figHTML[m] = widgets.HTML()
        eqBox[m] = widgets.Latex()
        if int(statesArray[m]) == 3: # m==0
            figHTML[m].value = '<img src="3-state-model.png" alt="Three state model" width=150>'
            eqBox[m].value = """
                    $$\\dot{C} = G_rD - \\epsilon F C$$
                    $$\\dot{O} = \epsilon FC -G_{d}O$$
                    $$\\dot{D} = G_{d}O-G_{r}D$$
                    $$C+O+D=1$$
                    $$\\epsilon F = \\phi\\frac{\\epsilon \\sigma_{ret}}{w_{loss}} = k\\phi$$
                    $$G_r = G_{r0} + G_{r1}(\\phi)$$
                    $$I_{\\phi} = g O (v-E)$$
                    """ #$$ $$
                    #G_r = G_{r,d} + \\mathcal{H}(\\phi) \\cdot G_{r,l}
                    #I_{\\phi} = \\bar{g} O \\cdot (v-E)
        elif int(statesArray[m]) == 4: # m==1
            figHTML[m].value = '<img src="4-state-model.png" alt="Four state model" width=180>' #width=200>'
            eqBox[m].value = """
                    $$\\dot{C_1} = G_rC_2 + G_{d1}(\\phi)O_1 - G_{a1}(\\phi)C_1$$
                    $$\\dot{O_1} = G_{a1}(\\phi)C_1 - (G_{d1}+e_{12}(\\phi))O_1 + e_{21}(\\phi)O_2$$
                    $$\\dot{O_2} = G_{a2}(\\phi)C_2 + e_{12}(\\phi)O_1 - (G_{d2}+e_{21}(\\phi))O_2$$
                    $$\\dot{C_2} = G_{d2}O_2 - (G_{a2}(\\phi)+G_r)C_2$$
                    $$C_1+O_1+O_2+C_2=1$$
                    $$$$
                    $$G_{a1}(\\phi) = \\phi\\frac{\\epsilon_1 \\sigma_{ret}}{w_{loss}} = k_1\\phi$$
                    $$G_{a2}(\\phi) = \\phi\\frac{\\epsilon_2 \\sigma_{ret}}{w_{loss}} = k_2\\phi$$
                    $$e_{12}(\\phi) = e_{12, d} + c_1 log(1+\\phi / \\phi_0)$$
                    $$e_{21}(\\phi) = e_{21, d} + c_2 log(1+\\phi / \\phi_0)$$
                    $$$$
                    $$I_{\\phi} = g (O_1+\\gamma O_2) (v-E)$$
                    """ #\\frac{\\phi}{\\phi_0}
        else: #int(statesArray[m]) == 6:
            figHTML[m].value = '<img src="http://link.springer.com/static-content/images/46/art%253A10.1007%252Fs10827-012-0431-7/MediaObjects/10827_2012_431_Fig1_HTML.gif" width=220>'
            eqBox[m].value = """
                    $$\\dot{C_1} = -a_1(\\phi)C_1 + b_1O_1 + a_6C_2$$
                    $$\dot{I_1} = a_1(\\phi)C_1 - a_2I_1$$
                    $$\\dot{O_1} = a_2I_1 - (b_1 + a_3(\\phi))O_1 + b_2(\\phi)O_2$$
                    $$\\dot{O_2} = a_3(\\phi)O_1 - (b_2(\\phi) + a_4)O_2 + b_3I_2$$
                    $$\dot{I_2} = -b_3I_2 + b_4(\\phi)C_2$$
                    $$\\dot{C_2} = a_4O_2 - (b_4(\\phi)+a_6)C_2$$
                    $$C_1+I_1+O_1+O_2+I_2+C_2=1$$
                    $$$$
                    $$a_1(\\phi) = a_{10}(\\phi / \\phi_0)$$
                    $$a_3(\\phi) = a_{30} + a_{31} \\ln(1 + \\phi / \\phi_0)$$
                    $$b_2(\\phi) = b_{20} + b_{21} \\ln(1 + \\phi / \\phi_0)$$
                    $$b_4(\\phi) = b_{40} (\\phi / \\phi_0)$$
                    $$$$
                    $$f(v) = \\frac{1-\\exp({-(v-E)/v_0})}{(v-E)/v_1}$$
                    $$I_{\\phi} = g (O_1+\\gamma O_2) f(v)(v-E)$$
                    """
        modelParamBoxes[m] = widgets.Box(children=pBoxArr[m])
        modelNotesBoxes[m] = widgets.HBox(children=[figHTML[m],eqBox[m]])
        #modelNotesBoxes[m].add_class('box-flex1')
        #modelBox = widgets.HBox(children=[modelParamBoxes,modelNotesBoxes])
        modelBoxes[m] = widgets.HBox(children=[modelParamBoxes[m],modelNotesBoxes[m]])#modelBox
        modelBoxes[m].margin = '5px'
    
    ### Linked parameters
    #E_Link = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    
    modelParamsTabs = widgets.Tab(description='Parameter Settings', children=modelBoxes, values=statesArray)
    modelParamsTabs.margin = '5px'
    modelParamsTabs.on_trait_change(modelTabChange,'selected_index')
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
    
    
    for sInd, sim in enumerate(simParams):
        pSet = simParams[sim]
        i=0
        for param in pSet:#.keys(): #, value in pSet.items():
            if isinstance(pSet[param].value, list):
                sim_pValArr[sInd][i] = widgets.Text(value=str(pSet[param].value), description=param) # np.asarray
            elif isinstance(pSet[param].value, str):
                sim_pValArr[sInd][i] = widgets.Text(value=str(pSet[param].value), description=param)
            elif isinstance(pSet[param].value, bool):
                sim_pValArr[sInd][i] = widgets.Dropdown(options=boolDict,value=pSet[param].value,description=param)
            else:
                if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf):
                    sim_pValArr[sInd][i] = widgets.FloatText(value=pSet[param].value, description=param)
                else:
                    sim_pValArr[sInd][i] = widgets.BoundedFloatText(value=pSet[param].value, min=pSet[param].min, max=pSet[param].max, description=param)
                sim_pValArr[sInd][i].width = '150px'
            if not pSet[param].expr == None:
                sim_unitArr[sInd][i] = widgets.Dropdown(options=[pSet[param].expr],value=pSet[param].expr)
                sim_pBoxArr[sInd][i] = widgets.HBox(children=[sim_pValArr[sInd][i],sim_unitArr[sInd][i]])
            else:
                sim_pBoxArr[sInd][i] = widgets.HBox(children=[sim_pValArr[sInd][i]])
            
            i+=1
        
        simFigHTML[sInd] = widgets.HTML()
        #exampleProt = '{}{}6s.{}'.format(fDir,prot,'png')#saveFigFormat)
        #if os.path.isfile(exampleProt):
        #    protFigHTML[pInd].value='<img src="{}" alt=Example {} width=200px>'.format(exampleProt,prot)
        #else:
        #    protFigHTML[pInd].value='Example Figure'
        simParamBoxes[sInd] = widgets.Box(children=sim_pBoxArr[sInd])
        simNotesBoxes[sInd] = widgets.HBox(children=[simFigHTML[sInd]])# simStimHTML[sInd]  , ])#[figHTML[prot],eqBox[prot]])
        
        
        simBoxes[sInd] = widgets.HBox(children=[simParamBoxes[sInd],simNotesBoxes[sInd]])#modelBox
        #display(protBoxes[pInd])
        simBoxes[sInd].margin = '5px'
    
    
    ##### Simulator parameters tab #####
    simParamsTabs = widgets.Tab(description='Simulator Settings', children=simBoxes)# \
    simParamsTabs.margin = '5px'
    simParamsTabs.on_trait_change(simTabChange,'selected_index')
    
    
    
    
    
    
    ##### Protocol parameters #####
    
    protBoxes = [None for prot in protList]
    protParamBoxes = [None for prot in protList] # Left side container
    protNotesBoxes = [None for prot in protList] # Right side container for figure and equations
    protStimHTML = [None for prot in protList]
    protFigHTML = [None for prot in protList] 
    #eqBox = [None for m in range(len(modelParams))]
    prot_pValArr = [[None for p in protParams[prot]] for prot in protList] # Array of parameter values
    prot_unitArr = [[None for p in protParams[prot]] for prot in protList] # Array of units
    prot_pBoxArr = [[None for p in protParams[prot]] for prot in protList] # Array of parameter boxes
    
    
    for pInd, prot in enumerate(protList):
        pSet = protParams[prot]
        i=0
        for param in pSet:#.keys(): #, value in pSet.items():
            if isinstance(pSet[param].value, list):
                prot_pValArr[pInd][i] = widgets.Text(value=str(pSet[param].value), description=param) # np.asarray
            else:
                if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf):
                    prot_pValArr[pInd][i] = widgets.FloatText(value=pSet[param].value, description=param)
                else:
                    prot_pValArr[pInd][i] = widgets.BoundedFloatText(value=pSet[param].value, min=pSet[param].min, max=pSet[param].max, description=param)
                prot_pValArr[pInd][i].width = '150px'
            if not pSet[param].expr == None:
                prot_unitArr[pInd][i] = widgets.Dropdown(options=[pSet[param].expr],value=pSet[param].expr)
                prot_pBoxArr[pInd][i] = widgets.HBox(children=[prot_pValArr[pInd][i],prot_unitArr[pInd][i]])
            else:
                prot_pBoxArr[pInd][i] = widgets.HBox(children=[prot_pValArr[pInd][i]])
            
            i+=1
        
        protStimHTML[pInd] = widgets.HTML()
        ### Add stimulus figures by printing to file
        # IPython.core.pylabtools.print_figure(fig, fmt='png', bbox_inches='tight', **kwargs)
        stimFig = plt.figure()
        x = np.linspace(0, 3*np.pi, 500)
        #plt.figure(Ifig.number)
        plt.plot(x, np.sin(x**2))
        fdata64 = base64.b64encode(print_figure(stimFig))
        title = 'Stimulus for {}'.format(prot)
        #print(title,fdata64)
        html_fig = ''#'<img alt="{}" src="data:image/png;base64,{}">'.format(title,fdata64)
        #return html_tpl.format(**locals())
        plt.close(stimFig)
        protStimHTML[pInd].value = html_fig
        
        protFigHTML[pInd] = widgets.HTML()
        exampleProt = '{}{}6s.{}'.format(fDir,prot,'png')#saveFigFormat)
        if os.path.isfile(exampleProt):
            protFigHTML[pInd].value='<img src="{}" alt=Example {} width=200px>'.format(exampleProt,prot)
        else:
            protFigHTML[pInd].value='Example Figure'
        protParamBoxes[pInd] = widgets.Box(children=prot_pBoxArr[pInd])
        protNotesBoxes[pInd] = widgets.HBox(children=[protStimHTML[pInd], protFigHTML[pInd]])#[figHTML[prot],eqBox[prot]])
        
        
        protBoxes[pInd] = widgets.HBox(children=[protParamBoxes[pInd],protNotesBoxes[pInd]])#modelBox
        protBoxes[pInd].margin = '5px'
        #display(protBoxes[pInd])




    ##### Protocol parameters tab #####
    protParamsTabs = widgets.Tab(description='Parameter Settings', children=protBoxes)# \
    protParamsTabs.margin = '5px'
    protParamsTabs.on_trait_change(protTabChange,'selected_index')
    
    
    ##### Configure tabs for abstraction layers #####
    paramTabs = widgets.Tab(description='Parameter Settings', children=[fitParamsTabs,modelParamsTabs,simParamsTabs,protParamsTabs])#,values=['Model', 'Protocol']) #E_box,k_box
    display(paramTabs)
    paramTabs.selected_index = TabGroups['Models'] # Set to show model parameters initially
    paramTabs.visible=False
    paramTabs.margin = '5px'
    
    
    
    ### Set Tab titles ### # Titles must be set after display of tabs
    paramTabs.set_title(0, 'Fit Parameters')
    paramTabs.set_title(1, 'Model Parameters')
    paramTabs.set_title(2, 'Simulator Parameters')
    paramTabs.set_title(3, 'Protocol Parameters')

    modelParamsTabs.set_title(0, 'Three-state model')
    modelParamsTabs.set_title(1, 'Four-state model')
    modelParamsTabs.set_title(2, 'Six-state model')
    modelParamsTabs.width = '800px' #set_css({'width': '800px'}) # 800
    
    ### Set Parameter Tabs style
    for pInd, prot in enumerate(protList): # Change to protParams #for prot, pSet in protParams.items():
        protParamsTabs.set_title(pInd, prot)
        protNotesBoxes[pInd].flex = True #protNotesBoxes[pInd].add_class('box-flex1')
        protFigHTML[pInd].width = '200px' #set_css({'width': '200px', 'margin-left': '20px'})
        protFigHTML[pInd].margin = '20px'
        for i in range(len(protParams[prot])): # pSet.keys(): #for i in range(len(protParams[prot])):
            prot_pValArr[pInd][i].width = '150px' #set_css({'width': '150px'})
    
    ### Set Simulator Tabs style
    for sInd, sim in enumerate(simList):
        simParamsTabs.set_title(sInd, sim)
        simNotesBoxes[sInd].flex = True #protNotesBoxes[pInd].add_class('box-flex1')
        simFigHTML[sInd].width = '200px' #set_css({'width': '200px', 'margin-left': '20px'})
        simFigHTML[sInd].margin = '20px'
        for i in range(len(simParams[sim])): # pSet.keys(): #for i in range(len(protParams[prot])):
            sim_pValArr[sInd][i].width = '150px' #set_css({'width': '150px'})
    
    
    ### Hack to tile parameter fields horizontally - must come after displaying the parent
    for model in modelParams: #for m in range(len(modelParams)):
        m = modelList.index(model)
        eqBox[m].width = '150px'
        eqBox[m].margin = '10px'
        modelNotesBoxes[m].width = '300px'
        modelNotesBoxes[m].margin = '20px'
    
    
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
    
    
    
    #GUI.children=[fitBar,runBar,paramsControlBar,paramTabs,NEURONbox,Brianbox]
    #display(GUI)
    #GUI.remove_class('vbox')
    #GUI.add_class('hbox')
    

    
    
    return

GUI = widgets.Box()
#interact(runModel, nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'saturate', 'inwardRect', 'varyPL', 'varyIPI'), saveData=True, verbose=1);

if __name__ == '__main__':
    loadGUI()