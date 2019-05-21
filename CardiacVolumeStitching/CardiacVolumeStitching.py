import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np


#
# CardiacVolumeStitching
#

class CardiacVolumeStitching(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Cardiac Volume Stitching"
        self.parent.categories = ["Cardiac"]
        self.parent.dependencies = []
        self.parent.contributors = [
            "Patrick Carnahan (Robarts Research Institute)"]  # replace with "Firstname Lastname (Organization)"
        self.parent.helpText = """
This module facilitates the registration and stitching or transgastric and TEE volumes.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""  # replace with organization, grant and thanks.


#
# CardiacVolumeStitchingWidget
#

class CardiacVolumeStitchingWidget(ScriptedLoadableModuleWidget):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModuleWidget.__init__(self, parent)

        # To store list of input nodes (volume or sequence)
        self.inputVolumeList = []

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        self._logic = CardiacVolumeStitchingLogic()

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
        # input volume folder selector
        #
        # self.inputDirSelector = ctk.ctkPathLineEdit()
        # self.inputDirSelector.filters = ctk.ctkPathLineEdit.Dirs
        # self.inputDirSelector.settingKey = 'Philips4dUsDicomPatcherInputDir'
        #
        # parametersFormLayout.addRow("Input Volume directory:", self.inputDirSelector)

        volumeSelectorHBox = qt.QHBoxLayout()

        #
        # input volume selector
        #
        self.inputSelector = slicer.qMRMLNodeComboBox()
        self.inputSelector.nodeTypes = ["vtkMRMLSequenceNode"]
        self.inputSelector.selectNodeUponCreation = True
        self.inputSelector.addEnabled = False
        self.inputSelector.removeEnabled = False
        self.inputSelector.noneEnabled = False
        self.inputSelector.showHidden = False
        self.inputSelector.showChildNodeTypes = False
        self.inputSelector.setMRMLScene(slicer.mrmlScene)
        self.inputSelector.setToolTip("Select input node to volume stitching.")
        volumeSelectorHBox.addWidget(self.inputSelector)

        self.addVolumeButton = qt.QPushButton('+')
        self.addVolumeButton.enabled = True
        self.addVolumeButton.toolTip = "Add the selected volume to the list."
        self.addVolumeButton.setMaximumWidth(40)
        volumeSelectorHBox.addWidget(self.addVolumeButton)

        self.removeVolumeButton = qt.QPushButton('-')
        self.removeVolumeButton.enabled = True
        self.removeVolumeButton.toolTip = "Remove the selected volume from the list"
        self.removeVolumeButton.setMaximumWidth(40)
        volumeSelectorHBox.addWidget(self.removeVolumeButton)

        parametersFormLayout.addRow(volumeSelectorHBox)

        #
        # Table with nodes to register and stitch
        #
        self.inputTable = qt.QTableWidget(0, 2)
        self.inputTable.setHorizontalHeaderLabels(('Volume', 'Position'))
        header = self.inputTable.horizontalHeader()
        header.setSectionResizeMode(0, qt.QHeaderView.Stretch)
        header.setSectionResizeMode(1, qt.QHeaderView.Fixed)
        header.resizeSection(1, 170)
        parametersFormLayout.addRow(self.inputTable)

        #
        # master sequence selector
        #
        self.masterSelector = slicer.qMRMLNodeComboBox()
        self.masterSelector.nodeTypes = ["vtkMRMLSequenceNode"]
        self.masterSelector.selectNodeUponCreation = True
        self.masterSelector.addEnabled = False
        self.masterSelector.removeEnabled = False
        self.masterSelector.renameEnabled = False
        self.masterSelector.noneEnabled = False
        self.masterSelector.showHidden = False
        self.masterSelector.showChildNodeTypes = False
        self.masterSelector.setMRMLScene(slicer.mrmlScene)
        self.masterSelector.setToolTip("Select master node for output sequence browser.")
        parametersFormLayout.addRow('Master Sequence', self.masterSelector)

        #
        # output volume selector
        #
        self.outputSelector = slicer.qMRMLNodeComboBox()
        self.outputSelector.nodeTypes = ["vtkMRMLSequenceNode"]
        self.outputSelector.selectNodeUponCreation = True
        self.outputSelector.addEnabled = True
        self.outputSelector.removeEnabled = True
        self.outputSelector.renameEnabled = True
        self.outputSelector.noneEnabled = True
        self.outputSelector.showHidden = False
        self.outputSelector.showChildNodeTypes = False
        self.outputSelector.baseName = 'Stitched-Output'
        self.outputSelector.setMRMLScene(slicer.mrmlScene)
        self.outputSelector.setToolTip("Select output node to volume stitching.")
        parametersFormLayout.addRow('Output Sequence', self.outputSelector)

        #
        # Apply Button
        #
        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.toolTip = "Run the algorithm."
        self.applyButton.enabled = False
        parametersFormLayout.addRow(self.applyButton)

        # connections
        self.addVolumeButton.connect('clicked(bool)', self.onAddVolumeButton)
        self.removeVolumeButton.connect('clicked(bool)', self.onRemoveVolumeButton)
        self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.masterSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.applyButton.connect('clicked(bool)', self.onApplyButton)

        # observers
        self.nodeRemovedObserver = slicer.mrmlScene.AddObserver(slicer.vtkMRMLScene.NodeAboutToBeRemovedEvent,
                                                                self.onNodeRemoved)

        # Add vertical spacer
        self.layout.addStretch(1)

    def cleanup(self):
        slicer.mrmlScene.RemoveObserver(self.nodeRemovedObserver)
        pass

    def onSelect(self):
        self.applyButton.enabled = self.outputSelector.currentNode() and self.inputVolumeList and len(
            self.inputVolumeList) > 1 and self.masterSelector.currentNode()

    def onAddVolumeButton(self):
        # Add volume to list and update table
        self.addVolume(self.inputSelector.currentNode())
        self.onSelect()

        return

    def addVolume(self, volume):
        if volume and volume.GetID() not in self.inputVolumeList:
            self.inputVolumeList.append(volume.GetID())

            node = self.inputSelector.currentNode()

            # Append row to end of table
            i = self.inputTable.rowCount
            self.inputTable.insertRow(i)

            # Set the node name
            self.inputTable.setCellWidget(i, 0, qt.QLabel(node.GetName()))

            volumeTypeDropdown = qt.QComboBox()
            presets = self._logic.getInitializationPresets()
            for j, preset in enumerate(presets):
                volumeTypeDropdown.addItem("{0} ({1})".format(*preset))

                if any([n in node.GetName() for n in preset]):
                    volumeTypeDropdown.currentIndex = j

            self.inputTable.setCellWidget(i, 1, volumeTypeDropdown)

    def onRemoveVolumeButton(self):
        # Remove volume from list and update table
        self.removeVolume(self.inputSelector.currentNode())
        self.onSelect()

        return

    def removeVolume(self, volume):
        if volume and volume.GetID() in self.inputVolumeList:
            i = self.inputVolumeList.index(volume.GetID())
            self.inputVolumeList.remove(volume.GetID())

            self.inputTable.removeRow(i)

    def onApplyButton(self):

        presets = [self.inputTable.cellWidget(i, 1).currentIndex for i in range(self.inputTable.rowCount)]
        self._logic.stitchSequences(self.inputVolumeList, presets, self.outputSelector.currentNode(),
                                    self.masterSelector.currentNode())

    @vtk.calldata_type(vtk.VTK_OBJECT)
    def onNodeRemoved(self, caller, event, calldata):
        node = calldata
        if isinstance(node, slicer.vtkMRMLSequenceNode):
            self.removeVolume(node)
            self.onSelect()


#
# CardiacVolumeStitchingLogic
#

class CardiacVolumeStitchingLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        self._parameterSets = [['Par0001Rigid.txt'], ['Par0001NonRigid.txt']]

        self._initializationPresets = [['Transesophageal', 'TEE'], ['Transgastric', 'TG'], ['Transthoracic', 'TTE']]

        # Create initialization transforms as identity
        self._initializationPresetTransforms = [None] * 3

        # TEE Transform
        tr = vtk.vtkTransform()
        tr.Identity()
        tr.RotateWXYZ(98.4210581181494, 0.35740674433659325, -0.8628562094610169, -0.35740674433659314)
        # Rotates the en face view to the correct orientation with respect to transgastric views
        self._initializationPresetTransforms[InitializationPresets_TEE] = tr

        # TG Transform
        # Uses identity
        tr = vtk.vtkTransform()
        tr.Identity()
        self._initializationPresetTransforms[InitializationPresets_TG] = tr

        # TODO TTE Transform
        tr = vtk.vtkTransform()
        tr.Identity()
        self._initializationPresetTransforms[InitializationPresets_TTE] = tr

    def getInitializationPresets(self):
        return self._initializationPresets

    def run(self):
        # Replace with list from ui
        nodes = [None] * 5
        nodes[0] = slicer.util.getFirstNodeByClassByName('vtkMRMLScalarVolumeNode', '011 Transgastric 1 Cartesian')
        nodes[1] = slicer.util.getFirstNodeByClassByName('vtkMRMLScalarVolumeNode', '011 Transgastric 2 Cartesian')
        nodes[2] = slicer.util.getFirstNodeByClassByName('vtkMRMLScalarVolumeNode', '011 Transgastric 3 Cartesian')
        nodes[3] = slicer.util.getFirstNodeByClassByName('vtkMRMLScalarVolumeNode', '011 Transgastric 4 Cartesian')

        tr = vtk.vtkTransform()
        # Rotates the en face view to the correct orientation with respect to transgastric views
        tr.RotateWXYZ(98.4210581181494, 0.35740674433659325, -0.8628562094610169, -0.35740674433659314)
        trnode = slicer.util.getFirstNodeByClassByName('vtkMRMLTransformNode', 'TEEProbeToRAS')
        if not trnode or not trnode.GetName() == 'TEEProbeToRAS':
            trnode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode', 'TEEProbeToRAS')

        trnode.SetAndObserveTransformToParent(tr)

        ef = slicer.util.getFirstNodeByClassByName('vtkMRMLScalarVolumeNode', '011 Cartesian DICOM')
        ef.SetAndObserveTransformNodeID(trnode.GetID())

        temp_node = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode')
        temp_node.Copy(ef)
        # temp_node.SetName('temp_volume')
        temp_node.SetAndObserveTransformNodeID(trnode.GetID())
        temp_node.HardenTransform()
        nodes[4] = temp_node

        masks = [None] * 5
        for i, n in enumerate(nodes):
            masks[i] = self.generateMask(n)

        self.registerVolumes(nodes, masks)
        # self.mergeVolumes(nodes)

        slicer.mrmlScene.RemoveNode(temp_node)
        slicer.mrmlScene.RemoveNode(trnode)

    def stitchSequences(self, sequences, presets, outputSequence, masterSequence):
        # Check for input parameters
        if not sequences or not presets:
            logging.debug("stitchSequences: Missing argument")
            return

        nodes = [slicer.mrmlScene.GetNodeByID(id) for id in sequences]

        # Clear output sequence
        outputSequence.RemoveAllDataNodes()
        outputSequence.SetIndexType(masterSequence.GetIndexType())
        outputSequence.SetIndexName(masterSequence.GetIndexName())
        outputSequence.SetIndexUnit(masterSequence.GetIndexUnit())

        # Set up browser
        seqBrowser = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSequenceBrowserNode")
        for n in nodes:
            slicer.modules.sequencebrowser.logic().AddSynchronizedNode(n, None, seqBrowser)

        if masterSequence not in nodes:
            slicer.modules.sequencebrowser.logic().AddSynchronizedNode(masterSequence, None, seqBrowser)
        seqBrowser.SetAndObserveMasterSequenceNodeID(masterSequence.GetID())

        # Create temporary volume for output
        outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")

        # Get proxy nodes and generate masks
        proxyNodes = [seqBrowser.GetProxyNode(n) for n in nodes]
        masks = [self.generateMask(n) for n in proxyNodes]
        masksTransformed = [slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode') for _ in masks]

        # Create temporary transform nodes
        initialTrNodes = []
        for tr in self._initializationPresetTransforms:
            trNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
            trNode.SetAndObserveTransformToParent(tr)
            initialTrNodes.append(trNode)

        # Create transform nodes to store rigid registrations
        rigidTrNodes = []
        for i in range(len(proxyNodes)):
            tr = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
            tr.SetAndObserveTransformToParent(vtk.vtkTransform())
            rigidTrNodes.append(tr)

        # Create final transform nodes
        nonRigidTrNodes = []
        for i in range(len(proxyNodes)):
            tr = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLTransformNode')
            tr.SetAndObserveTransformToParent(vtk.vtkTransform())
            nonRigidTrNodes.append(tr)

        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            numberOfDataNodes = masterSequence.GetNumberOfDataNodes()

            seqBrowser.SetSelectedItemNumber(2)
            slicer.modules.sequencebrowser.logic().UpdateProxyNodesFromSequences(seqBrowser)

            # Harden initial transforms for proxy nodes so that the transformed volume is used
            for i, n in enumerate(proxyNodes):
                n.SetAndObserveTransformNodeID(initialTrNodes[presets[i]].GetID())
                n.HardenTransform()
            # Transform masks for use in registration
            for i, n in enumerate(masksTransformed):
                n.Copy(masks[i])
                n.SetAndObserveTransformNodeID(initialTrNodes[presets[i]].GetID())
                n.HardenTransform()

            # Rigid registration on single frame
            self.registerVolumes(proxyNodes, rigidTrNodes, masksTransformed, 0)

            # Convert displacement fields to linear transforms
            for n in rigidTrNodes:
                self.getElastixRigidFromDisplacementGrid(n, n)

            # Concatenate additional transforms
            for i in range(len(rigidTrNodes)):
                tr = rigidTrNodes[i].GetTransformToParent()
                if tr:
                    tr.PreMultiply()
                    tr.Concatenate(self._initializationPresetTransforms[presets[i]])
                    tr.PostMultiply()
                    if i > 0:
                        tr.Concatenate(rigidTrNodes[i-1].GetTransformToParent())
                    tr.Update()

            # Transform masks for use in non rigid phase of registration
            for i, n in enumerate(masksTransformed):
                n.Copy(masks[i])
                n.SetAndObserveTransformNodeID(rigidTrNodes[i].GetID())
                n.HardenTransform()

            # Loop through sequence browser
            #
            # Performs non-rigid registration to account for
            # small deviations using rigid registration for initialization
            for seqItemNumber in range(numberOfDataNodes):
                print('Started frame: {0} of {1}'.format(seqItemNumber + 1, numberOfDataNodes))
                slicer.app.processEvents(qt.QEventLoop.ExcludeUserInputEvents)
                seqBrowser.SetSelectedItemNumber(seqItemNumber)
                slicer.modules.sequencebrowser.logic().UpdateProxyNodesFromSequences(seqBrowser)

                # Register Volumes

                # Harden rigid transforms for this iteration
                for i, n in enumerate(proxyNodes):
                    n.SetAndObserveTransformNodeID(rigidTrNodes[i].GetID())
                    n.HardenTransform()

                self.registerVolumes(proxyNodes, maskList=masksTransformed, parSet=1)

                # Merge Volumes
                self.mergeVolumes(proxyNodes, outputVolume)

                # Saved stitched result
                outputSequence.SetDataNodeAtValue(outputVolume,
                                                  seqBrowser.GetMasterSequenceNode().GetNthIndexValue(seqItemNumber))

                print('Completed frame: {0} of {1}'.format(seqItemNumber + 1, numberOfDataNodes))

        finally:
            qt.QApplication.restoreOverrideCursor()

            # Clean up temporary nodes
            slicer.mrmlScene.RemoveNode(outputVolume)
            seqBrowser.RemoveAllSequenceNodes()
            slicer.mrmlScene.RemoveNode(seqBrowser)

            for trNode in initialTrNodes:
                slicer.mrmlScene.RemoveNode(trNode)

            for trNode in rigidTrNodes:
                slicer.mrmlScene.RemoveNode(trNode)

            for trNode in nonRigidTrNodes:
                slicer.mrmlScene.RemoveNode(trNode)

            for mask in masks:
                if mask and mask.GetDisplayNode():
                    slicer.mrmlScene.RemoveNode(mask.GetDisplayNode().GetColorNode())
                slicer.mrmlScene.RemoveNode(mask)

            for mask in masksTransformed:
                if mask and mask.GetDisplayNode():
                    slicer.mrmlScene.RemoveNode(mask.GetDisplayNode().GetColorNode())
                slicer.mrmlScene.RemoveNode(mask)


            # Create sequence browser for output if it does not already exist
            if not slicer.modules.sequencebrowser.logic().GetFirstBrowserNodeForSequenceNode(outputSequence):
                outputSequenceBrowser = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSequenceBrowserNode",
                                                                           outputSequence.GetName() + ' browser')
                slicer.modules.sequencebrowser.logic().AddSynchronizedNode(outputSequence, None, outputSequenceBrowser)

        return

    def registerVolumes(self, volumeList, trNodes=None, maskList=None, parSet=None):

        if not volumeList:
            logging.debug('registerVolumes: missing parameter')
            return

        import Elastix
        elastix = Elastix.ElastixLogic()
        # elastix.logStandardOutput = True
        # elastix.deleteTemporaryFiles = False

        moduleDir = os.path.dirname(os.path.abspath(__file__))
        registrationResourcesDir = os.path.abspath(os.path.join(moduleDir, 'Resources', 'RegistrationParameters'))

        elastix.registrationParameterFilesDir = registrationResourcesDir

        for i in range(len(volumeList) - 1):

            outputVolume = None

            if trNodes and not isinstance(trNodes, (list, tuple, np.ndarray)):
                trNode = trNodes
            elif trNodes and len(trNodes) > i + 1:
                trNode = trNodes[i + 1]
            else:
                trNode = None
                outputVolume = volumeList[i+1]

            if maskList and not isinstance(maskList, (list, tuple, np.ndarray)):
                fixedMask = maskList
                movingMask = None
            elif maskList and len(maskList) == 1:
                fixedMask = maskList[0]
                movingMask = None
            elif maskList and len(maskList) > i + 1:
                fixedMask = maskList[i]
                movingMask = maskList[i + 1]
            else:
                fixedMask = None
                movingMask = None

            if parSet and not isinstance(parSet, (list, tuple, np.ndarray)):
                parIdx = parSet
            elif parSet and len(parSet) == 1:
                parIdx = parSet[0]
            elif parSet and len(parSet) > i:
                parIdx = parSet[i]
            else:
                parIdx = 0

            elastix.registerVolumes(volumeList[i], volumeList[i + 1],
                                    self._parameterSets[parIdx],
                                    outputTransformNode=trNode,
                                    outputVolumeNode=outputVolume,
                                    fixedVolumeMaskNode=fixedMask,
                                    movingVolumeMaskNode=movingMask)


    def mergeVolumes(self, volumeList, out):
        # Get bounds of final volume
        bounds = np.zeros((len(volumeList), 6))
        for i in range(len(volumeList)):
            volumeList[i].GetSliceBounds(bounds[i, :], None)

        min = bounds.min(0)
        max = bounds.max(0)

        volumeBounds_ROI = np.array([min[0], max[1], min[2], max[3], min[4], max[5]])

        outputSpacing = 0.4

        roiDim = np.zeros(3)

        for i in range(3):
            roiDim[i] = (volumeBounds_ROI[i * 2 + 1] - volumeBounds_ROI[i * 2]) / outputSpacing;

        roiDim = np.ceil(roiDim).astype('uint16')

        out.SetOrigin(volumeBounds_ROI[0], volumeBounds_ROI[2], volumeBounds_ROI[4])
        out.SetSpacing([outputSpacing] * 3)

        blend = vtk.vtkImageBlend()
        blend.SetBlendModeToCompound()

        # for each volume, resample into output space
        # will need to do for each frame in sequence as well
        for n in volumeList:
            # Get transforms for input and output volumes
            inputIJK2RASMatrix = vtk.vtkMatrix4x4()
            n.GetIJKToRASMatrix(inputIJK2RASMatrix)
            referenceRAS2IJKMatrix = vtk.vtkMatrix4x4()
            out.GetRASToIJKMatrix(referenceRAS2IJKMatrix)
            inputRAS2RASTransform = vtk.vtkGeneralTransform()
            if n.GetTransformNodeID():
                slicer.mrmlScene.GetNodeByID(n.GetTransformNodeID()).GetTransformToWorld(inputRAS2RASTransform)

            # Create resample transform from reference volume to input
            resampleTransform = vtk.vtkGeneralTransform()
            resampleTransform.Identity()
            resampleTransform.PostMultiply()
            resampleTransform.Concatenate(inputIJK2RASMatrix)
            resampleTransform.Concatenate(inputRAS2RASTransform)
            resampleTransform.Concatenate(referenceRAS2IJKMatrix)
            resampleTransform.Inverse()

            # Resample the image to the output space using transform
            resampler = vtk.vtkImageReslice()
            resampler.SetInputConnection(n.GetImageDataConnection())
            resampler.SetOutputOrigin(0, 0, 0)
            resampler.SetOutputSpacing(1, 1, 1)
            resampler.SetOutputExtent(0, roiDim[0], 0, roiDim[1], 0, roiDim[2])
            resampler.SetResliceTransform(resampleTransform)
            resampler.SetInterpolationModeToCubic()
            resampler.SetOutputScalarType(vtk.VTK_DOUBLE)
            resampler.Update()

            # Get mask using threshold
            thresh = vtk.vtkImageThreshold()
            thresh.ThresholdByUpper(2)
            thresh.SetOutputScalarTypeToDouble()
            thresh.SetOutValue(0)
            thresh.SetInValue(1)
            thresh.ReplaceInOn()
            thresh.ReplaceOutOn()
            thresh.SetInputConnection(resampler.GetOutputPort())
            thresh.Update()

            # Erode mask to eliminate boundary artifacts
            dilateErode = vtk.vtkImageDilateErode3D()
            dilateErode.SetDilateValue(0)
            dilateErode.SetErodeValue(1)
            dilateErode.SetKernelSize(7, 7, 7)
            dilateErode.SetInputConnection(thresh.GetOutputPort())
            dilateErode.Update()

            # Append mask to image as alpha component
            app = vtk.vtkImageAppendComponents()
            app.AddInputData(resampler.GetOutput())
            app.AddInputData(dilateErode.GetOutput())
            app.Update()

            blend.AddInputData(app.GetOutput())

        # Set output volume
        blend.Update()

        # Extract first component to get image and cast to uchar
        extract = vtk.vtkImageExtractComponents()
        extract.SetComponents(0)
        extract.SetInputConnection(blend.GetOutputPort())
        extract.Update()

        cast = vtk.vtkImageCast()
        cast.SetOutputScalarTypeToUnsignedChar()
        cast.SetInputConnection(extract.GetOutputPort())
        cast.Update()

        outImage = vtk.vtkImageData()
        outImage.DeepCopy(cast.GetOutput())
        out.SetAndObserveImageData(outImage)

    def getElastixRigidFromDisplacementGrid(self, gridTransformNode, outputNode):

        gridTransform = gridTransformNode.GetTransformFromParentAs('vtkOrientedGridTransform')
        if not gridTransform:
            # Not stored as a displacement grid from parent
            return

        # Get center and dimensions of displacement grid
        center = np.zeros(3)
        bounds = np.zeros(6)

        displacementGrid = gridTransform.GetDisplacementGrid()
        displacementGrid.GetCenter(center)
        displacementGrid.GetBounds(bounds)

        dims = np.zeros(3)
        for i in range(3):
            dims[i] = (bounds[i * 2 + 1] - bounds[i * 2]);

        # Create sample points in grid
        pointSrc = vtk.vtkPointSource()
        pointSrc.SetCenter(center)
        pointSrc.SetRadius(np.min((dims) / 2))
        pointSrc.SetNumberOfPoints(200)
        pointSrc.Update()

        # Transform a copy of points using displacement grid
        toPoints = pointSrc.GetOutput().GetPoints()
        fromPointsGrid = vtk.vtkPoints()

        gridTransform.TransformPoints(toPoints, fromPointsGrid)

        # Find corresponding rigid transform between point sets
        landmark = vtk.vtkLandmarkTransform()
        landmark.SetModeToRigidBody()
        landmark.SetSourceLandmarks(fromPointsGrid)
        landmark.SetTargetLandmarks(toPoints)
        landmark.Update()

        error = 0
        fromPointsRigid = vtk.vtkPoints()
        landmark.TransformPoints(fromPointsGrid, fromPointsRigid)
        for i in range(toPoints.GetNumberOfPoints()):
            toPoint = np.zeros(3)
            fromPoint = np.zeros(3)
            toPoints.GetPoint(i, toPoint)
            fromPointsRigid.GetPoint(i, fromPoint)

            error += np.linalg.norm(toPoint - fromPoint) ** 2

        error = np.sqrt(error / toPoints.GetNumberOfPoints())

        linTransform = vtk.vtkTransform()
        linTransform.Identity()
        if error > 0.5:
            logging.debug("Cannot compute rigid transform from given displacement field")
        else:
            linTransform.SetMatrix(landmark.GetMatrix())

        outputNode.SetAndObserveTransformToParent(linTransform)

    def generateMask(self, refVolumeNode):

        segmentationNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', 'VolumeMask')
        segmentationNode.GetSegmentation().AddEmptySegment('VolumeMask')

        # Create segment editor to get access to effects
        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
        segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
        segmentEditorNode.SetOverwriteMode(segmentEditorNode.OverwriteNone)
        segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
        segmentEditorWidget.setSegmentationNode(segmentationNode)
        segmentEditorWidget.setMasterVolumeNode(refVolumeNode)
        segmentEditorWidget.setCurrentSegmentID('VolumeMask')

        # Threshold to get a mask of the ultrasound cone
        segmentEditorWidget.setActiveEffectByName("Threshold")
        effect = segmentEditorWidget.activeEffect()
        masterImageData = effect.masterVolumeImageData()
        lo, hi = masterImageData.GetScalarRange()
        effect.setParameter("MinimumThreshold", 1)
        effect.setParameter("MaximumThreshold", hi)
        effect.self().onApply()

        # Margin Shrink to exclude borders
        segmentEditorWidget.setActiveEffectByName("Margin")
        effect = segmentEditorWidget.activeEffect()
        effect.setParameter("MarginSizeMm", -6)
        effect.self().onApply()

        # Create label map node to store mask
        labelNode = slicer.util.getFirstNodeByClassByName('vtkMRMLLabelMapVolumeNode',
                                                          refVolumeNode.GetName() + '-label')
        if not labelNode:
            labelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode',
                                                           refVolumeNode.GetName() + '-label')
        labelNode.SetAndObserveTransformNodeID(refVolumeNode.GetTransformNodeID())

        segmentationIds = vtk.vtkStringArray()
        segmentationIds.InsertNextValue('VolumeMask')
        slicer.modules.segmentations.logic().ExportSegmentsToLabelmapNode(segmentationNode, segmentationIds,
                                                                          labelNode, refVolumeNode)
        # Clean up
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segmentEditorNode)
        slicer.mrmlScene.RemoveNode(segmentationNode)

        return labelNode


class CardiacVolumeStitchingTest(ScriptedLoadableModuleTest):
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


InitializationPresets_TEE = 0
InitializationPresets_TG = 1
InitializationPresets_TTE = 2
