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
statesDict = OrderedDict([(' '+s,i) for i,s in enumerate(list(modelParams.keys()))]) # enumerate(modelList)
#statesDict = OrderedDict([(' 3',0), (' 4',1), (' 6',2)]) 
statesArray = list(statesDict.keys()) #[' 3', ' 4', ' 6'] # [u' 3',u' 4',u' 6']

paramsGroups = {'Models':0, 'Protocols':1}

clearDelay = 1.5 # Pause [s] before clearing text entry fields

#verbose = 1

# Plot settings
#plotPeakRecovery = False
#plotStateVars = False
#plotKinetics = False

#p0fV = (40,4,1)#25000)#,1)      # v0,v1,E,G
#p0IPI = (-1e-8,400,-1e-7) # a*exp(-t/b)+c





    

    
    


    
    

def loadGUI():
    
    
    ##### Model fitting bar Functions #####

    def fitToggle(name, value):
        if value == True:
            fitBar.visible = True
            time.sleep(clearDelay) # Pause then clear the input field
            dataVar.value = ''
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
        global fitParamsPopup
        global clearOutput
        if clearOutput.value:
            clear_output()
            if 'fitParamsPopup' in vars() or 'fitParamsPopup' in globals():
                fitParamsPopup.close()
        global fitRhO
        fitRhO = fitModels(dataSet, int(statesToFitButtons.value))
        fitParamReport = widgets.TextareaWidget(description='Report:',value=fitRhO.reportParams())
        fitParamsPopup = widgets.PopupWidget(children=[fitParamReport],button_text='Fitted Parameters',description='Fitted {} state model parameters from: {}'.format(int(statesToFitButtons.value),dataVar.value))
        display(fitParamsPopup)
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

    #firstCall = True
    def paramsToggle(name, value):
        #global firstCall
        #if firstCall == True: # Workaround to fix tab titles not displaying correctly: https://github.com/ipython/ipython/issues/5743
            #protParamsTabs.selected_index = 1
            #protParamsTabs.selected_index = 0
        #    paramTabs.selected_index = paramsGroups['Models'] # Start with model parameters
            #modelParamsTabs.selected_index = 1
            #modelParamsTabs.selected_index = 0
        #    firstCall = False
        
        if value == True: # Set to current model and protocol tabs
            paramsControlBar.visible = True
            paramTabs.visible = True
            modelParamsTabs.selected_index = statesDict[stateButtons.value]
            protParamsTabs.selected_index = protList.index(protDropdown.value) #protIndDict[protDropdown.value]
            time.sleep(clearDelay) # Pause then clear the input field
            paramVar.value = ''
        else:
            paramsControlBar.visible = False
            paramTabs.visible = False
        return

    def outputDropdownChange(name,value):
        if value == 'NEURON':
            NEURONbox.visible=True
            Brianbox.visible=False
            time.sleep(clearDelay) # Pause then clear the input field
            hocFile.value = ''
        elif value == 'Brian':
            Brianbox.visible=True
            NEURONbox.visible=False
            time.sleep(clearDelay) # Pause then clear the input field
            brianFile.value = ''
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
        #RhO = selectModel(nStates)
        #if paramsButton.value == True:
        #print(RhO)
        userModelParams = Parameters()
        mInd = statesDict[model] #statesDict[' '+str(nStates)]
        pSet = modelParams[modelList[mInd]]           
        #nParams = len(pSet.keys())
        i=0
        for key in pSet.keys():  # Set of n-state model pararmeters
            userModelParams.add(key, value=pValArr[mInd][i].value)
            i+=1
        #RhO.setParams(userModelParams)
        return userModelParams

    def setGUImodelParams(pSet,model):
        mInd = statesDict[model] #statesDict[' '+str(nStates)]
        #pSet = modelParams[mInd]            # Set of n-state model pararmeters
        #nParams = len(pSet.keys())
        i=0
        for key in pSet.keys(): #, value in pSet.items():
            pValArr[mInd][i].value = pSet[key].value
            i+=1
        return
        
    def getGUIprotParams(protocol):
        userProtParams = Parameters()
        pInd = protList.index(protocol) #protIndDict[protocol]
        pSet = protParams[protocol]
        i=0
        for param in pSet.keys():
            if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                userProtParams.add(param, value=ast.literal_eval(prot_pValArr[pInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
            else:
                userProtParams.add(param, value=prot_pValArr[pInd][i].value)
            i+=1
        #Prot = protocols[protocol](userProtParams)
        return userProtParams
        
    def setGUIprotParams(pSet,protocol):
        #pSet = protParams[protocol]
        pInd = protList.index(protocol) #protIndDict[protocol]
        i=0
        for param in pSet.keys(): #, value in pSet.items():
            #print(param)
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
            #if paramsButton.value == True:
            print(RhO)
            userModelParams = Parameters()
            mInd = statesDict[model] #statesDict[' '+str(nStates)]
            pSet = modelParams[modelList[mInd]]           
            nParams = len(pSet.keys())
            i=0
            for key in pSet.keys():  # Set of n-state model pararmeters
                userModelParams.add(key,value=pValArr[mInd][i].value)
                #print(key,pValArr[mInd][i].value,type(pValArr[mInd][i].value))
                i+=1
            RhO.setParams(userModelParams)
        
        #dataTag = str(nStates)+"s"
        
        #%run -i modProtocols.py
        #%run -i models.py
        
        if output == 'Python':
            #Prot = protocols[protocol]()
            #Prot = selectProtocol(protocol)
            
            pInd = protList.index(protocol) #protIndDict[protocol]
            userProtParams = Parameters()
            pSet = protParams[protocol]
            i=0
            for param in pSet.keys(): #, value in pSet.items():
                #print(protocol,pInd,param,i)
                if isinstance(pSet[param].value, list): ################################## Change to not number!!!
                    #print(prot_pValArr[pInd][i].value, type(prot_pValArr[pInd][i].value))
                    userProtParams.add(param, value=ast.literal_eval(prot_pValArr[pInd][i].value)) # or http://stackoverflow.com/questions/5387208/convert-a-string-to-an-array
                
                #    print('TextBox: ',param,' = ',str(pSet[param].value))
                else:
                    userProtParams.add(param, value=prot_pValArr[pInd][i].value)
                i+=1
            Prot = protocols[protocol](userProtParams)
            #Prot.setParams(userProtParams)
            #Prot.printParams()
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
    fitButton = widgets.ToggleButtonWidget(description='Fit data', value=False)
    fitButton.on_trait_change(fitToggle,'value') #paramsButton.get_state()

    ### Create Data set entry
    dataVar = widgets.TextWidget(description='Data Set: ',value='<variable name>')
    #dataVar.set_css({'width': '200px'})
    #def dataVarClick(name):
    #    dataVar.value = ''
    #dataVar.on_click(dataVarClick)
    dataLoad = widgets.ButtonWidget(description='Load')
    dataLoad.on_click(onLoad)
    
    ### Create Fit Model States buttons
    statesToFitButtons = widgets.ToggleButtonsWidget(description='Highest model to fit: ',values=statesArray)#,value=u' 3') #https://github.com/ipython/ipython/issues/6469
    #fit3sCheck = widgets.CheckboxWidget(description='Fit 3 state model', value=True)
    #fit4sCheck = widgets.CheckboxWidget(description='Fit 4 state model', value=False)
    #fit6sCheck = widgets.CheckboxWidget(description='Fit 6 state model', value=False)
    
    ### Create Checkboxes
    runSSAcheck = widgets.CheckboxWidget(description='Characterise', value=True)
    useFitCheck = widgets.CheckboxWidget(description='Use fit', value=False)
    plotExpData = widgets.CheckboxWidget(description='Plot data: ', value=True)
    
    ### Create Run Button
    runFitButton = widgets.ButtonWidget(description="Fit!")
    runFitButton.on_click(runFitButton_on_click)
    
    ### Create Fit Bar
    fitBar = widgets.ContainerWidget(children=[dataVar, dataLoad, statesToFitButtons, runSSAcheck, useFitCheck, plotExpData, runFitButton]) #fit3sCheck, fit4sCheck, fit6sCheck
    display(fitBar)
    fitBar.visible=False
    
    ### Set formatting
    dataVar.set_css({'width': '150px'})
    dataLoad.add_class('btn-warning')
    runFitButton.add_class('btn-danger')
    # Set Fit Bar formating after display
    fitBar.remove_class('vbox')
    fitBar.add_class('hbox')
    fitBar.remove_class('align-start')
    fitBar.add_class('align-center')
    
    
    
    ##### Run bar #####
    
    ### Protocol Dropdown
    protDropdown = widgets.DropdownWidget(values=protList,value='saturate') #protDict
    protDropdown.on_trait_change(protDropdownChange,'value')
    
    ### Run Button
    runButton = widgets.ButtonWidget(description="Run!")
    
    runButton.on_click(runButton_on_click)
    
    ### Model states button
    
    stateButtons = widgets.ToggleButtonsWidget(description='Model states: ',values=statesArray) #https://github.com/ipython/ipython/issues/6469
    stateButtons.on_trait_change(changeModel,'value')#stateButtons.value_name  #'_view_name'
    
    ### Create states buttons link
    modelParamsLink = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    #modelParamsLink = link((stateButtons, 'value'), (modelParamsTabs, 'value'))
    
    ### Parameters Toggle
    paramsButton = widgets.ToggleButtonWidget(description='Parameters', value=False)
    paramsButton.on_trait_change(paramsToggle,'value')
    
    ### Load NEURON parameters
    hocFile = widgets.TextWidget(description='HOC file: ',value='<file name>')
    hocLoad = widgets.ButtonWidget(description='Load')
    hocLoad.on_click(onHocLoad)
    
    ### Load Brian parameters
    brianFile = widgets.TextWidget(description='Brian file: ',value='<file name>')
    brianLoad = widgets.ButtonWidget(description='Load')
    brianLoad.on_click(onBrianLoad)
    
    ### Output format Dropdown
    outputDropdown = widgets.DropdownWidget(values=['Python', 'NEURON', 'Brian'],value='Python')
    outputDropdown.on_trait_change(outputDropdownChange,'value')
    
    
    
    ### Model figure popup

    #figButton = widgets.ButtonWidget(description="Figure")
    #counter = widgets.IntTextWidget(description='Counter:')
    #popup = widgets.PopupWidget(children=[counter], description='Popup Demo', button_text='Popup Button')
    #display(popup)

    #popup.set_css({
    #    'width': '420px',
    #    'height': '350px'}, selector='modal')
    #graph_widget = ForceDirectedGraphWidget(graph)
    #popup.children = [graph_widget]


    #figHTML = widgets.HTMLWidget()

    #figPopup = widgets.PopupWidget(children=[figHTML])
    #figPopup.button_text = "Figure"
    #def figButton_on_click(b): # on_button_clicked
    #    #stateButtons.value = statesArray[value]
    #    if stateButtons.value == statesArray[0] or stateButtons.value == statesArray[1]: # 3 or 4 state
    #        figPopup.description = "Three and Four state model figures"
    #        figHTML.value = '<img src="http://api.onlinelibrary.wiley.com/asset/v1/doi/10.1111%2Fj.1751-1097.2008.00460.x/asset/image_n%2FPHP_460_f1.gif?l=j6%2BNsqLlmq%2F8Bqdug7r6wFbpY%2B373Yz%2Fg76IMstuTvGHfO6q%2BFbQB8c%2ByMjlpZ6tKRQQk3vBGWG3%0AN2tSYhODvw%3D%3D&o=e8c686ea-f58f-49bb-a923-369073ce3542">'
    #    else:
    #        figPopup.description = "Six state model figure"
    #        figHTML.value = '<img src="http://link.springer.com/static-content/images/46/art%253A10.1007%252Fs10827-012-0431-7/MediaObjects/10827_2012_431_Fig1_HTML.gif">'

    #    #display(figPopup)
    #figPopup.on_displayed(figButton_on_click)
    #figPopup.close()
    
    
    
    ### Plot and Save button
    saveButton = widgets.CheckboxWidget(description='Plot & Save', value=True)
    saveButton.visible = False

    ### Verbosity slider
    #verboseSlide = widgets.FloatProgressWidget(value=1, min=-1, max=3, step=1, description='Verbosity:')
    verboseSlide = widgets.IntSliderWidget(value=1, min=0, max=3, description='Output:')
    verboseSlide.visible=True #False

    ### Clear ouput checkbox
    clearOutput = widgets.CheckboxWidget(description='Clear Output', value=True)
    
    
    #### Run Bar container
    runBar = widgets.ContainerWidget(children=[fitButton,stateButtons,protDropdown,paramsButton,outputDropdown,saveButton,verboseSlide,clearOutput,runButton]) #???
    display(runBar)
    #container.border_color = 'red'
    #container.border_style = 'dotted'
    #container.border_width = 3
    #display(container)
    
    
    ### Set Button styles after displaying the runBar
    paramsButton.add_class('btn-info')
    fitButton.add_class('btn-success')
    runButton.add_class('btn-danger')
    
    #stateButtons.add_class('btn-primary')
    #protDropdown.add_class('btn-success')
    #outputDropdown.add_class('btn-warning')
    
    ### button bootstrap styles: default, primary, success, info, warning, danger, link
    
    runBar.remove_class('vbox')
    runBar.add_class('hbox')
    runBar.remove_class('align-start')
    runBar.add_class('align-center')
    #outputDropdown.set_css({'width': '60px'})
    verboseSlide.set_css({'width': '60px'})
    
    ##### Parameters control bar #####
    
    paramVar = widgets.TextWidget(description='Parameter Set: ',value='<variable name>')
    paramLoad = widgets.ButtonWidget(description='Load')
    paramLoad.on_click(paramSetOnLoad)
    
    
    ### Reset Parameters Button
    paramsResetButton = widgets.ButtonWidget(description="Reset")
    paramsResetButton.on_click(resetButton_on_click)
    
    ### Create Parameters control bar
    paramsControlBar = widgets.ContainerWidget(children=[paramVar,paramLoad,paramsResetButton])
    display(paramsControlBar)
    paramsControlBar.visible = False
    
    
    ### Set Parameters Bar style
    paramVar.set_css({'width': '150px'})
    paramLoad.add_class('btn-warning')
    paramsResetButton.add_class('btn-warning')
    
    paramsControlBar.remove_class('vbox')
    paramsControlBar.add_class('hbox')
    paramsControlBar.remove_class('align-start')
    paramsControlBar.add_class('align-center')
    
    
    
    
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
        nParams = len(pSet.keys())
        
        i=0
        for key in pSet.keys(): #, value in pSet.items():
            if pSet[key].min == None or pSet[key].max == None:
                pValArr[m][i] = widgets.FloatTextWidget(value=pSet[key].value, description=key)#"{} [{}]:".format(key, pSet[key].expr))
            else:
                pValArr[m][i] = widgets.BoundedFloatTextWidget(value=pSet[key].value, min=pSet[key].min, max=pSet[key].max, description=key)#"{} [{}]:".format(key, pSet[key].expr))
            if not pSet[key].expr == None:
                unitArr[m][i] = widgets.DropdownWidget(values=[pSet[key].expr],value=pSet[key].expr)
                pBoxArr[m][i] = widgets.ContainerWidget(children=[pValArr[m][i],unitArr[m][i]])
            else:
                pBoxArr[m][i] = widgets.ContainerWidget(children=[pValArr[m][i]])
            
            pBoxArr[m][i].remove_class('vbox')
            pBoxArr[m][i].add_class('hbox')
            #pBoxArr[m][i].visible = False
            i+=1
        
        figHTML[m] = widgets.HTMLWidget()
        eqBox[m] = widgets.LatexWidget()
        if int(statesArray[m]) == 3: # m==0
            ### Popup widget
            
            #figPopup = widgets.PopupWidget(children=[figHTML])
            #figPopup.button_text = "Figure"
            #figPopup.description = "Three and Four state model figures"
            figHTML[m].value = '<img src="3-state-model.png" alt="Three state model" width=135>' #width=150>' # height="42" width="42"
            #'<img src="http://api.onlinelibrary.wiley.com/asset/v1/doi/10.1111%2Fj.1751-1097.2008.00460.x/asset/image_n%2FPHP_460_f1.gif?l=j6%2BNsqLlmq%2F8Bqdug7r6wFbpY%2B373Yz%2Fg76IMstuTvGHfO6q%2BFbQB8c%2ByMjlpZ6tKRQQk3vBGWG3%0AN2tSYhODvw%3D%3D&o=e8c686ea-f58f-49bb-a923-369073ce3542">'
            
            #def figButton_on_click(b): # on_button_clicked
                #stateButtons.value = statesArray[value]
            #    if stateButtons.value == statesArray[0] or stateButtons.value == statesArray[1]: # 3 or 4 state
            #        figPopup.description = "Three and Four state model figures"
            #        figHTML.value = '<img src="http://api.onlinelibrary.wiley.com/asset/v1/doi/10.1111%2Fj.1751-1097.2008.00460.x/asset/image_n%2FPHP_460_f1.gif?l=j6%2BNsqLlmq%2F8Bqdug7r6wFbpY%2B373Yz%2Fg76IMstuTvGHfO6q%2BFbQB8c%2ByMjlpZ6tKRQQk3vBGWG3%0AN2tSYhODvw%3D%3D&o=e8c686ea-f58f-49bb-a923-369073ce3542">'
            #    else:
            #        figPopup.description = "Six state model figure"
            #        figHTML.value = '<img src="http://link.springer.com/static-content/images/46/art%253A10.1007%252Fs10827-012-0431-7/MediaObjects/10827_2012_431_Fig1_HTML.gif">'

                #display(figPopup)
            #figPopup.on_displayed(figButton_on_click)
            #figPopup.close()
            
            #modelpArrBox = widgets.ContainerWidget(children=pBoxArr[m])
            #modelBox = widgets.ContainerWidget(children=[modelpArrBox,figPopup])
            #RhO = selectModel(3)
            #eqBox[m].value = RhO_3states().equations
            eqBox[m].value = """$$ $$
                    $$\\dot{C} = G_rD - \\epsilon F C$$
                    $$\\dot{O} = \epsilon FC -G_{d}O$$
                    $$\\dot{D} = G_{d}O-G_{r}D$$
                    $$C+O+D=1$$
                    $$\\epsilon F = \\phi\\frac{\\epsilon \\sigma_{ret}}{w_{loss}} = k\\phi$$
                    $$G_r = G_{r0} + G_{r1}(\\phi)$$
                    $$I_{\\phi} = g O (v-E)$$
                    """
                    #G_r = G_{r,d} + \\mathcal{H}(\\phi) \\cdot G_{r,l}
                    #I_{\\phi} = \\bar{g} O \\cdot (v-E)
        elif int(statesArray[m]) == 4: # m==1
            figHTML[m].value = '<img src="4-state-model.png" alt="Four state model" width=180>' #width=200>'
        else: #int(statesArray[m]) == 6:
            figHTML[m].value = '<img src="http://link.springer.com/static-content/images/46/art%253A10.1007%252Fs10827-012-0431-7/MediaObjects/10827_2012_431_Fig1_HTML.gif" width=220>'
            #modelBox = widgets.ContainerWidget(children=pBoxArr[m])
        
        
        
        modelParamBoxes[m] = widgets.ContainerWidget(children=pBoxArr[m])
        modelNotesBoxes[m] = widgets.ContainerWidget(children=[figHTML[m],eqBox[m]])
        #modelNotesBoxes[m].add_class('box-flex1')
        #modelBox = widgets.ContainerWidget(children=[modelParamBoxes,modelNotesBoxes])
        modelBoxes[m] = widgets.ContainerWidget(children=[modelParamBoxes[m],modelNotesBoxes[m]])#modelBox
    
    ### Linked parameters
    #E_Link = link((stateButtons, 'value'), (statesToFitButtons, 'value'))
    
    modelParamsTabs = widgets.TabWidget(description='Parameter Settings', children=modelBoxes, values=statesArray)
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
    
    #pInd = 0
    #for prot, pSet in protParams.items():
    for pInd, prot in enumerate(protList):
        pSet = protParams[prot]
        i=0
        for param in pSet.keys(): #, value in pSet.items():
            #print(param)
            if isinstance(pSet[param].value, list):
                prot_pValArr[pInd][i] = widgets.TextWidget(value=str(pSet[param].value), description=param) # np.asarray
            #    print('TextBox: ',param,' = ',str(pSet[param].value))
            else:
                if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf):
                    prot_pValArr[pInd][i] = widgets.FloatTextWidget(value=pSet[param].value, description=param)
                else:
                    #print(param,pSet[param].value,'[',pSet[param].min,pSet[param].max,']')
                    prot_pValArr[pInd][i] = widgets.BoundedFloatTextWidget(value=pSet[param].value, min=pSet[param].min, max=pSet[param].max, description=param)
            if not pSet[param].expr == None:
                prot_unitArr[pInd][i] = widgets.DropdownWidget(values=[pSet[param].expr],value=pSet[param].expr)
                prot_pBoxArr[pInd][i] = widgets.ContainerWidget(children=[prot_pValArr[pInd][i],prot_unitArr[pInd][i]])
            else:
                prot_pBoxArr[pInd][i] = widgets.ContainerWidget(children=[prot_pValArr[pInd][i]])
            
            i+=1
        
        protStimHTML[pInd] = widgets.HTMLWidget()
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
        
        protFigHTML[pInd] = widgets.HTMLWidget()
        exampleProt = '{}{}6s.{}'.format(fDir,prot,'png')#saveFigFormat)
        if os.path.isfile(exampleProt):
            protFigHTML[pInd].value='<img src="{}" alt=Example {}>'.format(exampleProt,prot)
        else:
            protFigHTML[pInd].value='Example Figure'
        protParamBoxes[pInd] = widgets.ContainerWidget(children=prot_pBoxArr[pInd])
        protNotesBoxes[pInd] = widgets.ContainerWidget(children=[protStimHTML[pInd], protFigHTML[pInd]])#[figHTML[prot],eqBox[prot]])
        
        
        protBoxes[pInd] = widgets.ContainerWidget(children=[protParamBoxes[pInd],protNotesBoxes[pInd]])#modelBox
        #display(protBoxes[pInd])
        #pInd+=1




    ##### Protocol parameters tab #####
    protParamsTabs = widgets.TabWidget(description='Parameter Settings', children=protBoxes)# \
    protParamsTabs.on_trait_change(protTabChange,'selected_index')
    
    paramTabs = widgets.TabWidget(description='Parameter Settings', children=[modelParamsTabs,protParamsTabs])#,values=['Model', 'Protocol']) #E_box,k_box
    display(paramTabs)
    paramTabs.selected_index = paramsGroups['Models'] # Set to show model parameters initially
    paramTabs.visible=False
    
    
    ### Set Tab titles ### # Titles must be set after display of tabs
    paramTabs.set_title(0, 'Model Parameters')
    paramTabs.set_title(1, 'Protocol Parameters')

    modelParamsTabs.set_title(0, 'Three-state model')
    modelParamsTabs.set_title(1, 'Four-state model')
    modelParamsTabs.set_title(2, 'Six-state model')
    modelParamsTabs.set_css({'width': '800px'}) # 800
    
    ### Set Parameter Tabs style
    for pInd, prot in enumerate(protList): # Change to protParams #for prot, pSet in protParams.items():
        protParamsTabs.set_title(pInd, prot)
        protNotesBoxes[pInd].add_class('box-flex1')
        protFigHTML[pInd].set_css({'width': '200px', 'margin-left': '20px'})
        for i in range(len(protParams[prot])): # pSet.keys(): #for i in range(len(protParams[prot])):
            #print(pInd,prot,i)
            prot_pBoxArr[pInd][i].remove_class('vbox')
            prot_pBoxArr[pInd][i].add_class('hbox')
            #if (pSet[param].min == None or pSet[param].min == -np.inf) or (pSet[param].max == None or pSet[param].max == np.inf):
            prot_pValArr[pInd][i].set_css({'width': '150px'})
        protBoxes[pInd].remove_class('vbox')
        protBoxes[pInd].add_class('hbox')

    ### Hack to tile parameter fields horizontally - must come after displaying the parent
    for model in modelParams: #for m in range(len(modelParams)):
        m = modelList.index(model)
        #modelParamsTabs.selected_index = m
        for p in range(len(modelParams[model])):
            pBoxArr[m][p].remove_class('vbox')
            pBoxArr[m][p].add_class('hbox')
            #if not unitArr[m][p] == None:
            #    unitArr[m][p].set_css({'width': '50px'})
            #pValArr[m][p].set_css({'width': '100px'})
        eqBox[m].set_css({'width': '150px', 'margin-left': '20px'}) #remove_class('vbox')
        #eqBox[m].add_class('hbox')
        modelBoxes[m].remove_class('vbox')
        modelBoxes[m].add_class('hbox')
        #modelNotesBoxes[m].add_class('box-flex1')
        modelNotesBoxes[m].set_css({'width': '300px', 'margin-left': '20px'})
        modelNotesBoxes[m].remove_class('vbox')
        modelNotesBoxes[m].add_class('hbox')
    
    
    ### Create NEURON Control Bar
    NEURONbox = widgets.ContainerWidget(children=[hocFile, hocLoad])
    display(NEURONbox)
    NEURONbox.visible=False

    NEURONbox.remove_class('vbox')
    NEURONbox.add_class('hbox')
    NEURONbox.remove_class('align-start')
    NEURONbox.add_class('align-center')

    hocLoad.add_class('btn-warning')
    
    
    ### Create Brian Control Bar
    Brianbox = widgets.ContainerWidget(children=[brianFile, brianLoad])
    display(Brianbox)
    Brianbox.visible=False

    Brianbox.remove_class('vbox')
    Brianbox.add_class('hbox')
    Brianbox.remove_class('align-start')
    Brianbox.add_class('align-center')

    brianLoad.add_class('btn-warning')
    
    #GUI.children=[fitBar,runBar,paramsControlBar,paramTabs,NEURONbox,Brianbox]
    #display(GUI)
    #GUI.remove_class('vbox')
    #GUI.add_class('hbox')
    

    
    
    return

GUI = widgets.ContainerWidget()
#interact(runModel, nStates={'Three-state':3,'Four-state':4,'Six-state':6}, protocol=('custom', 'step', 'sinusoid', 'ramp', 'saturate', 'inwardRect', 'varyPL', 'varyIPI'), saveData=True, verbose=1);