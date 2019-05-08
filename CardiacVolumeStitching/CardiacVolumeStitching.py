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
        self.parent.contributors = ["Patrick Carnahan (Robarts Research Institute)"] # replace with "Firstname Lastname (Organization)"
        self.parent.helpText = """
This module facilitates the registration and stitching or transgastric and TEE volumes.
"""
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# CardiacVolumeStitchingWidget
#

class CardiacVolumeStitchingWidget(ScriptedLoadableModuleWidget):
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
        # input volume folder selector
        #
        # self.inputDirSelector = ctk.ctkPathLineEdit()
        # self.inputDirSelector.filters = ctk.ctkPathLineEdit.Dirs
        # self.inputDirSelector.settingKey = 'Philips4dUsDicomPatcherInputDir'
        #
        # parametersFormLayout.addRow("Input Volume directory:", self.inputDirSelector)

        #
        # output volume selector
        #
        self.outputSelector = slicer.qMRMLNodeComboBox()
        self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.outputSelector.selectNodeUponCreation = True
        self.outputSelector.addEnabled = True
        self.outputSelector.removeEnabled = True
        self.outputSelector.noneEnabled = False
        self.outputSelector.showHidden = False
        self.outputSelector.showChildNodeTypes = False
        self.outputSelector.setMRMLScene(slicer.mrmlScene)
        self.outputSelector.setToolTip("Pick the output of the algorithm.")
        parametersFormLayout.addRow("Output Volume", self.outputSelector)


        #
        # Apply Button
        #
        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.toolTip = "Run the algorithm."
        self.applyButton.enabled = True
        parametersFormLayout.addRow(self.applyButton)

        # connections
        self.applyButton.connect('clicked(bool)', self.onApplyButton)

        # Add vertical spacer
        self.layout.addStretch(1)


    def cleanup(self):
        pass



    def onApplyButton(self):
        logic = CardiacVolumeStitchingLogic()
        logic.run()


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
        self.parameterSets = [['Par0001Rigid.txt', 'Par0001NonRigid.txt'], ['Par0002Rotation.txt', 'Par0001Rigid.txt']]

    def run(self):
        # Replace with list from ui
        nodes = [None]*5
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
        temp_node.SetName('temp_volume')
        temp_node.SetAndObserveTransformNodeID(trnode.GetID())
        temp_node.HardenTransform()
        nodes[4] = temp_node

        masks = [None]*5
        for i,n in enumerate(nodes):
            masks[i] = self.generateMask(n)

        self.registerVolumes(nodes, masks)
        self.mergeVolumes(nodes)

        slicer.mrmlScene.RemoveNode(temp_node)
        slicer.mrmlScene.RemoveNode(trnode)

    def mergeVolumes(self, volumeList):
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

        out = slicer.util.getFirstNodeByName('Stitched-Output')
        if not out or not out.GetName() == 'Stitched-Output':
            out = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode', 'Stitched-Output')

        out.SetOrigin(volumeBounds_ROI[0], volumeBounds_ROI[2], volumeBounds_ROI[4])
        out.SetSpacing([outputSpacing] * 3)

        # Create accumulator image initialized to 0 to populate with final max values
        maxImage = vtk.vtkImageData()
        maxImage.SetOrigin(0,0,0)
        maxImage.SetSpacing(1, 1, 1)
        maxImage.SetExtent(0, roiDim[0], 0, roiDim[1], 0, roiDim[2])
        maxImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,0)

        # for each volume, resample into output space
        # will need to do for each frame in sequence as well
        for i in range(len(volumeList)):
            n = volumeList[i]

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
            resampler.SetOutputScalarType(vtk.VTK_UNSIGNED_CHAR)
            resampler.Update()

            # Take maximum value
            mathFilter = vtk.vtkImageMathematics()
            mathFilter.SetOperationToMax()
            mathFilter.SetInput1Data(maxImage)
            mathFilter.SetInput2Data(resampler.GetOutput())
            mathFilter.Update()

            maxImage.DeepCopy(mathFilter.GetOutput())

        # Set output volume
        out.SetAndObserveImageData(maxImage)

    def registerVolumes(self, volumeList, maskList = None, parSet = None):
        import Elastix
        elastix = Elastix.ElastixLogic()
        # elastix.logStandardOutput = True
        # elastix.deleteTemporaryFiles = False

        moduleDir = os.path.dirname(os.path.abspath(__file__))
        registrationResourcesDir = os.path.abspath(os.path.join(moduleDir, 'Resources', 'RegistrationParameters'))

        elastix.registrationParameterFilesDir = registrationResourcesDir

        for i in range(len(volumeList)-1):
            trName = 'Volume{0}To{1}'.format(i+2,i+1)
            trNode = slicer.util.getFirstNodeByName(trName)
            if not trNode or not trNode.GetName() == trName:
                trNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLGridTransformNode', trName)

            if maskList and not isinstance(maskList, (list, tuple, np.ndarray)):
                fixedMask = maskList
                movingMask = None
            elif maskList and len(maskList) == 1:
                fixedMask = maskList[0]
                movingMask = None
            elif maskList and len(maskList) > i+1:
                fixedMask = maskList[i]
                movingMask = maskList[i+1]
            else:
                fixedMask = None
                movingMask = None

            if parSet and len(parSet) == 1:
                parIdx = parSet[0]
            elif parSet and len(parSet) > i:
                parIdx = parSet[i]
            else:
                parIdx = 0


            elastix.registerVolumes(volumeList[i], volumeList[i+1],
                                    self.parameterSets[parIdx],
                                    outputTransformNode = trNode,
                                    fixedVolumeMaskNode = fixedMask,
                                    movingVolumeMaskNode = movingMask)

            volumeList[i+1].SetAndObserveTransformNodeID(trNode.GetID())

            trNode.SetAndObserveTransformNodeID(volumeList[i].GetTransformNodeID())


    def getElastixRigidFromDisplacementGrid(self, gridTransformNode, outputNode):

        gridTransform = gridTransformNode.GetTransformFromParent()

        # Get center and dimensions of displacement grid
        center = np.zeros(3)
        bounds = np.zeros(6)

        gridTransform.GetDisplacementGrid().GetCenter(center)
        gridTransform.GetDisplacementGrid().GetBounds(bounds)

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
        labelNode = slicer.util.getFirstNodeByClassByName('vtkMRMLLabelMapVolumeNode', refVolumeNode.GetName()+'-label')
        if not labelNode:
            labelNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', refVolumeNode.GetName()+'-label')
        labelNode.SetAndObserveTransformNodeID(refVolumeNode.GetTransformNodeID())

        segmentationIds = vtk.vtkStringArray()
        segmentationIds.InsertNextValue('VolumeMask')
        slicer.modules.segmentations.logic().ExportSegmentsToLabelmapNode(segmentationNode, segmentationIds,
                                                                          labelNode, refVolumeNode)
        # Clean up
        segmentEditorWidget = None
        slicer.mrmlScene.RemoveNode(segmentEditorNode)
        slicer.mrmlScene.RemoveNode(segmentationNode)

        return  labelNode


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


