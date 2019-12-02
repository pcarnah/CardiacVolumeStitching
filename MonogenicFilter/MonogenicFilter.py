#  /*====================================================================
#  Copyright (C) 2019  Patrick Carnahan <pcarnah@uwo.ca>
#  Unauthorized copying of this file, via any medium is strictly prohibited
#  Proprietary and confidential
#  ====================================================================*/

import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
from MonogenicSignal import MonogenicSignal
import sitkUtils
import SimpleITK as sitk


#
# MonogenicSignal
#

class MonogenicFilter(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Monogenic Signal"
        self.parent.categories = ["Examples"]
        self.parent.dependencies = []
        self.parent.contributors = ["Patrick Carnahan (Robarts Research Institute)"]
        self.parent.helpText = """
This module computes the monogenic signal fro a given input volume.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""  # replace with organization, grant and thanks.


#
# MonogenicFilterWidget
#

class MonogenicFilterWidget(ScriptedLoadableModuleWidget):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # Instantiate and connect widgets ...

        #
        # Parameters Area
        #
        parametersCollapsibleButton = ctk.ctkCollapsibleButton()
        parametersCollapsibleButton.text = "Parameters"
        self.layout.addWidget(parametersCollapsibleButton)

        # Layout within the dummy collapsible button
        parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

        #
        # input volume selector
        #
        self.inputSelector = slicer.qMRMLNodeComboBox()
        self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.inputSelector.selectNodeUponCreation = True
        self.inputSelector.addEnabled = False
        self.inputSelector.removeEnabled = False
        self.inputSelector.noneEnabled = False
        self.inputSelector.showHidden = False
        self.inputSelector.showChildNodeTypes = False
        self.inputSelector.setMRMLScene(slicer.mrmlScene)
        self.inputSelector.setToolTip("Pick the input to the algorithm.")
        parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

        #
        # output volume selector
        #
        self.outputSelector = slicer.qMRMLNodeComboBox()
        self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.outputSelector.selectNodeUponCreation = True
        self.outputSelector.addEnabled = True
        self.outputSelector.removeEnabled = True
        self.outputSelector.noneEnabled = True
        self.outputSelector.showHidden = False
        self.outputSelector.showChildNodeTypes = False
        self.outputSelector.setMRMLScene(slicer.mrmlScene)
        self.outputSelector.setToolTip("Pick the output to the algorithm.")
        parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

        #
        # Mode selector
        #
        modeHBox = qt.QHBoxLayout()

        fsButton = qt.QRadioButton("Feature Symmetry")
        fsButton.checked = True
        leButton = qt.QRadioButton("Local Energy")
        loButton = qt.QRadioButton("Local Orientation")
        lpButton = qt.QRadioButton("Local Phase")
        osButton = qt.QRadioButton("Oriented Symmetry")

        self.modeSelector = qt.QButtonGroup()
        self.modeSelector.addButton(fsButton)
        self.modeSelector.setId(fsButton, MonogenicFilterLogic.FEATURE_SYMMETRY)
        self.modeSelector.addButton(leButton)
        self.modeSelector.setId(leButton, MonogenicFilterLogic.LOCAL_ENERGY)
        self.modeSelector.addButton(loButton)
        self.modeSelector.setId(loButton, MonogenicFilterLogic.LOCAL_ORIENTATION)
        self.modeSelector.addButton(lpButton)
        self.modeSelector.setId(lpButton, MonogenicFilterLogic.LOCAL_PHASE)
        self.modeSelector.addButton(osButton)
        self.modeSelector.setId(osButton, MonogenicFilterLogic.ORIENTED_SYMMETRY)

        modeHBox.addWidget(fsButton)
        modeHBox.addWidget(leButton)
        modeHBox.addWidget(loButton)
        modeHBox.addWidget(lpButton)
        modeHBox.addWidget(osButton)

        parametersFormLayout.addRow(modeHBox)

        #
        # Apply Button
        #
        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.toolTip = "Run the algorithm."
        self.applyButton.enabled = False
        parametersFormLayout.addRow(self.applyButton)

        # connections
        self.applyButton.connect('clicked(bool)', self.onApplyButton)
        self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        pass

    def onSelect(self):
        self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

    def onApplyButton(self):
        logic = MonogenicFilterLogic()
        logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), self.modeSelector.checkedId())


#
# MonogenicFilterLogic
#

class MonogenicFilterLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    FEATURE_SYMMETRY = 0
    LOCAL_ENERGY = 1
    LOCAL_ORIENTATION = 2
    LOCAL_PHASE = 3
    ORIENTED_SYMMETRY = 4

    def run(self, inputVolume, outputVolume, mode=FEATURE_SYMMETRY):
        """
        Run the actual algorithm
        """

        im = sitkUtils.PullVolumeFromSlicer(inputVolume)
        im2 = sitk.SmoothingRecursiveGaussian(im, 1.0)

        volume = sitk.GetArrayFromImage(im2)

        # volume = slicer.util.arrayFromVolume(inputVolume)
        volume = volume.transpose([2,1,0])

        filters = MonogenicSignal.MonogenicFilters(volume.shape, im.GetSpacing(), [6, 12, 24, 32])
        monogenic = filters.getMonogenicSignal(volume)

        if mode == MonogenicFilterLogic.FEATURE_SYMMETRY:
            outVolume, _ = monogenic.featureSymmetry()
        elif mode == MonogenicFilterLogic.LOCAL_ENERGY:
            outVolume = monogenic.localEnergy()[:,:,:,2]
        elif mode == MonogenicFilterLogic.LOCAL_ORIENTATION:
            outVolume = monogenic.localOrientation()[:,:,:,2]
        elif mode == MonogenicFilterLogic.LOCAL_PHASE:
            outVolume = monogenic.localPhase()[:,:,:,2]
        elif mode == MonogenicFilterLogic.ORIENTED_SYMMETRY:
            _, _, _, outVolume = monogenic.orientedSymmetry()
        else:
            logging.warning("Invalid mode.")
            return False

        outVolume = outVolume.transpose([2, 1, 0])

        outputVolume.CopyOrientation(inputVolume)
        slicer.util.updateVolumeFromArray(outputVolume, outVolume)

        return True


class MonogenicFilterTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """ Do whatever is needed to reset the state - typically a scene clear will be enough.
        """
        slicer.mrmlScene.Clear(0)

    def runTest(self):
        """Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_MonogenicFilter1()

    def test_MonogenicFilter1(self):
        """ Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")
        #
        # first, get some data
        #
        import SampleData
        SampleData.downloadFromURL(
            nodeNames='FA',
            fileNames='FA.nrrd',
            uris='http://slicer.kitware.com/midas3/download?items=5767',
            checksums='SHA256:12d17fba4f2e1f1a843f0757366f28c3f3e1a8bb38836f0de2a32bb1cd476560')
        self.delayDisplay('Finished with download and loading')

        volumeNode = slicer.util.getNode(pattern="FA")
        logic = MonogenicFilterLogic()
        self.assertIsNotNone(logic.hasImageData(volumeNode))
        self.delayDisplay('Test passed!')
