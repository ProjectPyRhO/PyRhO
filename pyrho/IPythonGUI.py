from IPython.html.widgets import interact
from IPython.html import widgets
from IPython.display import display
from IPython.display import clear_output
from IPython.utils.traitlets import link
from IPython.core.pylabtools import print_figure
import base64

from .models import *
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

paramsGroups = {'Models':0, 'Protocols':1}

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
        paramTabs.selected_index = paramsGroups['Protocols'] # Set to protocol parameters
        protParamsTabs.selected_index = protList.index(value) # Get the index of the selected protocol #protIndDict[value]


    def runButton_on_click(b): # on_button_clicked
        if clearOutput.value:
            clear_output()
        runModel(stateButtons.value, protDropdown.value, outputDropdown.value, saveButton.value, verboseSlide.value) #int(stateButtons.value)
        return        
        
    def changeModel(name,value):
        paramTabs.selected_index = paramsGroups['Models'] # Set to model parameters
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

    def outputDropdownChange(name,value):
        if value == 'NEURON':
            NEURONbox.visible=True
            Brianbox.visible=False
        elif value == 'Brian':
            Brianbox.visible=True
            NEURONbox.visible=False
        else: # Set both to invisible
            NEURONbox.visible=False
            Brianbox.visible=False
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
            
        if paramTabs.selected_index == paramsGroups['Models']:
            mInd = modelParamsTabs.selected_index
            model = statesArray[mInd]
            #mInd = statesDict[' '+str(nStates)]
            #pSet = modelParams[mInd] # Default model parameters
            setGUImodelParams(paramSet,model)
        else: # Protocols tab selected
            pInd = protParamsTabs.selected_index
            #pInd = protIndDict[protocol]
            protocol = protList[pInd]
            #pSet = protParams[protocol] # Default protocol parameters
            setGUIprotParams(paramSet,protocol)
        return paramSet
        
    def resetButton_on_click(b): # on_button_clicked
        if paramTabs.selected_index == paramsGroups['Models']:
            mInd = modelParamsTabs.selected_index
            model = statesArray[mInd]
            #mInd = statesDict[' '+str(nStates)]
            pSet = modelParams[mInd] # Default model parameters
            setGUImodelParams(pSet,model)
        else: # Protocols tab selected
            pInd = protParamsTabs.selected_index
            #pInd = protIndDict[protocol]
            protocol = protList[pInd]
            pSet = protParams[protocol] # Default protocol parameters
            setGUIprotParams(pSet,protocol)
        return
    
    def modelTabChange(name,value):
        stateButtons.value = statesArray[value] #' '+str(value)

    def protTabChange(name,value):
        protDropdown.value_name = protParamsTabs._titles[value] #[protParamsTabs.selected_index]
    
    
    ### Utility functions 
    def getGUImodelParams(model):
        userModelParams = Parameters()
        mInd = statesDict[model] #statesDict[' '+str(nStates)]
        pSet = modelParams[modelList[mInd]]           
        i=0
        for key in pSet:  # Set of n-state model pararmeters
            userModelParams.add(key, value=pValArr[mInd][i].value)
            i+=1
        return userModelParams

    def setGUImodelParams(pSet,model):
        mInd = statesDict[model] #statesDict[' '+str(nStates)]
        i=0
        for key in pSet: #, value in pSet.items():
            pValArr[mInd][i].value = pSet[key].value
            i+=1
        return
        
    def getGUIprotParams(protocol):
        userProtParams = Parameters()
        pInd = protList.index(protocol) #protIndDict[protocol]
        pSet = protParams[protocol]
        i=0
        for param in pSet: #.keys():
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                userProtParams.add(param, value=ast.literal_eval(prot_pValArr[pInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            else:
                userProtParams.add(param, value=prot_pValArr[pInd][i].value)
            i+=1
        return userProtParams
        
    def setGUIprotParams(pSet,protocol):
        pInd = protList.index(protocol) #protIndDict[protocol]
        i=0
        for param in pSet: #.keys(): #, value in pSet.items():
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                prot_pValArr[pInd][i].value = str(pSet[param].value)
            else:
                prot_pValArr[pInd][i].value = pSet[param].value
            i+=1
        return

        
        
    ##### Main Run function #####
    
        #saveData = True
        #@interact(nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'saturate', 'inwardRect', 'varyPL', 'varyIPI'), saveData=True, verbose=1)

    def runModel(model, protocol, output='Python', saveData=True, verbose=1): #verboseSlide.value
        
        verbose = verbose
        nStates = int(model)
        print("\nRunning Protocol '{}' on the {}-state model...".format(protocol,nStates))
        print('--------------------------------------------------------------------------------\n')
        
        ### Choose between 3, 4, & 6 state models
        if fitButton.value and useFitCheck.value:
            RhO = fitRhO
            
            ### To set GUI parameters, widget arrays should be dictionaries so subsets of parameters can be set in arbitrary order!!!
            #setGUIprotParams(expProtParams,'custom')
            
        else:
            RhO = selectModel(int(model))
            print(RhO)
            userModelParams = Parameters()
            mInd = statesDict[model] #statesDict[' '+str(nStates)]
            pSet = modelParams[modelList[mInd]]           
            nParams = len(pSet)#.keys())
            i=0
            for key in pSet:#.keys():  # Set of n-state model pararmeters
                userModelParams.add(key,value=pValArr[mInd][i].value)
                i+=1
            RhO.setParams(userModelParams)
        
        
        if output == 'Python':
            pInd = protList.index(protocol) #protIndDict[protocol]
            userProtParams = Parameters()
            pSet = protParams[protocol]
            i=0
            for param in pSet:#.keys(): #, value in pSet.items():
                if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                    userProtParams.add(param, value=ast.literal_eval(prot_pValArr[pInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
                else:
                    userProtParams.add(param, value=prot_pValArr[pInd][i].value)
                i+=1
            Prot = protocols[protocol](userProtParams)
            Prot.runProtocol(RhO,verbose)
            if verbose > 0: #saveData:
                Prot.plotProtocol(RhO,verbose)
        elif output == 'NEURON':
            print('Simulating with NEURON...')
        elif output == 'Brian':
            print('Simulating with Brian...')
        else:
            print('Error! Unknown output selection!')
        print("\nFinished!")
        print('================================================================================\n\n')
        
        return #Prot, RhO
        
        
        

        
    
    ##### Fit Bar #####
    
    ### Create Fit Bar toggle
    fitButton = widgets.ToggleButton(description='Fit data', value=False)
    fitButton.on_trait_change(fitToggle,'value') #paramsButton.get_state()

    ### Create Data set entry
    dataVar = widgets.Text(description='Data Set: ',placeholder='<variable name>')
    dataLoad = widgets.Button(description='Load')
    dataLoad.on_click(onLoad)
    
    ### Create Fit Model States buttons
    statesToFitButtons = widgets.ToggleButtons(description='Highest model to fit: ',options=statesArray)#,value=u' 3') #https://github.com/ipython/ipython/issues/6469
    #fit3sCheck = widgets.Checkbox(description='Fit 3 state model', value=True)
    #fit4sCheck = widgets.Checkbox(description='Fit 4 state model', value=False)
    #fit6sCheck = widgets.Checkbox(description='Fit 6 state model', value=False)
    
    ### Create Checkboxes
    runSSAcheck = widgets.Checkbox(description='Characterise', value=True)
    useFitCheck = widgets.Checkbox(description='Use fit', value=False)
    plotExpData = widgets.Checkbox(description='Plot data: ', value=True)
    
    ### Create Run Button
    runFitButton = widgets.Button(description="Fit!")
    runFitButton.on_click(runFitButton_on_click)
    
    ### Create Fit Bar
    fitBar = widgets.HBox(children=[dataVar, dataLoad, statesToFitButtons, runSSAcheck, useFitCheck, plotExpData, runFitButton]) #fit3sCheck, fit4sCheck, fit6sCheck
    display(fitBar)
    fitBar.visible=False
    
    ### Set formatting
    dataVar.width = '150px' #set_css({'width': '150px'})
    dataLoad.button_style = 'warning' #add_class('btn-warning')
    runFitButton.button_style = 'danger' #add_class('btn-danger')
    # Set Fit Bar formating after display
    fitBar.align = 'center'
    
    
    ##### Run bar #####
    
    ### Protocol Dropdown
    protDropdown = widgets.Dropdown(options=protList,value='saturate')###########,value='saturate') #protDict
    protDropdown.on_trait_change(protDropdownChange,'value')
    
    ### Run Button
    runButton = widgets.Button(description="Run!")
    
    runButton.on_click(runButton_on_click)
    
    ### Model states button
    
    stateButtons = widgets.ToggleButtons(description='Model states: ',options=statesArray) #https://github.com/ipython/ipython/issues/6469
    stateButtons.on_trait_change(changeModel,'value')#stateButtons.value_name  #'_view_name'
    
    ### Create states buttons link
    modelParamsLink = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    ### Parameters Toggle
    paramsButton = widgets.ToggleButton(description='Parameters', value=False)
    paramsButton.on_trait_change(paramsToggle,'value')
    
    ### Load NEURON parameters
    hocFile = widgets.Text(description='HOC file: ',placeholder='<file name>')
    hocLoad = widgets.Button(description='Load')
    hocLoad.on_click(onHocLoad)
    
    ### Load Brian parameters
    brianFile = widgets.Text(description='Brian file: ',placeholder='<file name>')
    brianLoad = widgets.Button(description='Load')
    brianLoad.on_click(onBrianLoad)
    
    ### Output format Dropdown
    outputDropdown = widgets.Dropdown(options=['Python', 'NEURON', 'Brian'],value='Python')
    outputDropdown.on_trait_change(outputDropdownChange,'value')
    
    
    
    
    ### Plot and Save button
    saveButton = widgets.Checkbox(description='Plot & Save', value=True)
    saveButton.visible = False

    ### Verbosity slider
    #verboseSlide = widgets.FloatProgressWidget(value=1, min=-1, max=3, step=1, description='Verbosity:')
    verboseSlide = widgets.IntSlider(value=1, min=0, max=3, description='Output:')
    verboseSlide.visible=True #False

    ### Clear ouput checkbox
    clearOutput = widgets.Checkbox(description='Clear Output', value=True)
    
    
    #### Run Bar container
    runBar = widgets.HBox(children=[fitButton,stateButtons,protDropdown,paramsButton,outputDropdown,saveButton,verboseSlide,clearOutput,runButton])
    display(runBar)

    
    
    ### Set Button styles after displaying the runBar
    paramsButton.button_style = 'info'
    fitButton.button_style = 'success'
    runButton.button_style = 'danger'
    
    #stateButtons.add_class('btn-primary')
    #protDropdown.add_class('btn-success')
    #outputDropdown.add_class('btn-warning')
    
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
        else: #int(statesArray[m]) == 6:
            figHTML[m].value = '<img src="http://link.springer.com/static-content/images/46/art%253A10.1007%252Fs10827-012-0431-7/MediaObjects/10827_2012_431_Fig1_HTML.gif" width=220>'
        
        modelParamBoxes[m] = widgets.Box(children=pBoxArr[m])
        modelNotesBoxes[m] = widgets.HBox(children=[figHTML[m],eqBox[m]])
        #modelNotesBoxes[m].add_class('box-flex1')
        #modelBox = widgets.HBox(children=[modelParamBoxes,modelNotesBoxes])
        modelBoxes[m] = widgets.HBox(children=[modelParamBoxes[m],modelNotesBoxes[m]])#modelBox
    
    ### Linked parameters
    #E_Link = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    
    modelParamsTabs = widgets.Tab(description='Parameter Settings', children=modelBoxes, values=statesArray)
    modelParamsTabs.margin = '5px'
    modelParamsTabs.on_trait_change(modelTabChange,'selected_index')
    #modelParamsTabs.on_displayed(setModelParamsTabs,'selected_index')
    #modelParamsTabs.on_displayed()
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    
    
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
        #display(protBoxes[pInd])




    ##### Protocol parameters tab #####
    protParamsTabs = widgets.Tab(description='Parameter Settings', children=protBoxes)# \
    protParamsTabs.margin = '5px'
    protParamsTabs.on_trait_change(protTabChange,'selected_index')
    
    paramTabs = widgets.Tab(description='Parameter Settings', children=[modelParamsTabs,protParamsTabs])#,values=['Model', 'Protocol']) #E_box,k_box
    display(paramTabs)
    paramTabs.selected_index = paramsGroups['Models'] # Set to show model parameters initially
    paramTabs.visible=False
    
    
    ### Set Tab titles ### # Titles must be set after display of tabs
    paramTabs.set_title(0, 'Model Parameters')
    paramTabs.set_title(1, 'Protocol Parameters')

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
        
    ### Hack to tile parameter fields horizontally - must come after displaying the parent
    for model in modelParams: #for m in range(len(modelParams)):
        m = modelList.index(model)
        eqBox[m].width = '150px'
        eqBox[m].margin = '10px'
        modelNotesBoxes[m].width = '300px'
        modelNotesBoxes[m].margin = '20px'
    
    
    ### Create NEURON Control Bar
    NEURONbox = widgets.HBox(children=[hocFile, hocLoad])
    display(NEURONbox)
    NEURONbox.visible=False
    NEURONbox.align = 'center'
    
    hocLoad.button_style = 'warning' #add_class('btn-warning')
    
    
    ### Create Brian Control Bar
    Brianbox = widgets.HBox(children=[brianFile, brianLoad])
    display(Brianbox)
    Brianbox.visible=False
    Brianbox.align = 'center'

    brianLoad.button_style = 'warning' #add_class('btn-warning')
    
    #GUI.children=[fitBar,runBar,paramsControlBar,paramTabs,NEURONbox,Brianbox]
    #display(GUI)
    #GUI.remove_class('vbox')
    #GUI.add_class('hbox')
    

    
    
    return

GUI = widgets.Box()
#interact(runModel, nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'saturate', 'inwardRect', 'varyPL', 'varyIPI'), saveData=True, verbose=1);