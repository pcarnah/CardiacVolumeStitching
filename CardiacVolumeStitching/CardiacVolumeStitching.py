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
        self.inputDirSelector = ctk.ctkPathLineEdit()
        self.inputDirSelector.filters = ctk.ctkPathLineEdit.Dirs
        self.inputDirSelector.settingKey = 'Philips4dUsDicomPatcherInputDir'

        parametersFormLayout.addRow("Input Volume directory:", self.inputDirSelector)


        #
        # Apply Button
        #
        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.toolTip = "Run the algorithm."
        self.applyButton.enabled = False
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

    def run(self):
        # still need to include running elastix from code

        # Replace with list from ui
        nodes = slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode')

        for n in nodes:
            if not "Transgastric" in n.GetName():
                nodes.remove(n)

        # Will be used, get bounds of final volume
        bounds = np.zeros((4, 6))
        for i in xrange(4):
            nodes[i].GetSliceBounds(bounds[i, :], None)

        min = bounds.min(0)
        max = bounds.max(0)

        volumeBounds_ROI = np.array([min[0], max[1], min[2], max[3], min[4], max[5]])

        outputSpacing = 0.4

        roiDim = np.zeros(3)

        for i in range(3):
            roiDim[i] = (volumeBounds_ROI[i * 2 + 1] - volumeBounds_ROI[i * 2]) / outputSpacing;

        roiDim = np.ceil(roiDim).astype('uint16')

        final = np.zeros([4, roiDim[2] + 1, roiDim[1] + 1, roiDim[0] + 1], 'uint8')

        # for each volume, resample into output space
        # will need to do for each frame in sequence as well
        for i in xrange(4):
            n = nodes[i]

            out = slicer.util.getFirstNodeByName(n.GetName() + '-Resampled')
            if not out:
                out = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode', n.GetName() + '-Resampled')

            out.SetOrigin(volumeBounds_ROI[0], volumeBounds_ROI[2], volumeBounds_ROI[4])
            out.SetSpacing([outputSpacing] * 3)

            inputIJK2RASMatrix = vtk.vtkMatrix4x4()
            n.GetIJKToRASMatrix(inputIJK2RASMatrix)
            referenceRAS2IJKMatrix = vtk.vtkMatrix4x4()
            out.GetRASToIJKMatrix(referenceRAS2IJKMatrix)
            inputRAS2RASMatrix = vtk.vtkGeneralTransform()
            if n.GetTransformNodeID():
                slicer.mrmlScene.GetNodeByID(n.GetTransformNodeID()).GetTransformToWorld(inputRAS2RASMatrix)

            resampleTransform = vtk.vtkGeneralTransform()
            resampleTransform.Identity()
            resampleTransform.PostMultiply()
            resampleTransform.Concatenate(inputIJK2RASMatrix)
            resampleTransform.Concatenate(inputRAS2RASMatrix)
            resampleTransform.Concatenate(referenceRAS2IJKMatrix)
            resampleTransform.Inverse()

            resampler = vtk.vtkImageReslice()
            resampler.SetInputConnection(n.GetImageDataConnection())
            resampler.SetOutputOrigin(0, 0, 0)
            resampler.SetOutputSpacing(1, 1, 1)
            resampler.SetOutputExtent(0, roiDim[0], 0, roiDim[1], 0, roiDim[2])
            resampler.SetResliceTransform(resampleTransform)
            resampler.SetInterpolationModeToCubic()
            resampler.Update()

            out.SetAndObserveImageData(resampler.GetOutput())

        # replace this with running max and only one output volume for space efficiency
        for i in xrange(4):
            n = nodes[i]

            out = slicer.util.getFirstNodeByName(n.GetName() + '-Resampled')
            final[i, :] = slicer.util.arrayFromVolume(out)

        finalMax = final.max(0)
        out = slicer.util.getFirstNodeByName('Stitched-Output')
        if not out:
            out = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode', 'Stitched-Output')

        out.SetOrigin(volumeBounds_ROI[0], volumeBounds_ROI[2], volumeBounds_ROI[4])
        out.SetSpacing([outputSpacing] * 3)
        slicer.util.updateVolumeFromArray(out, finalMax)


    def getElastixRigidFromDisplacementGrid(self, gridTransformNode, outputNode):

        gridTransform = gridTransformNode.GetTransformFromParent()

        # Get center and dimensions of displacement grid
        center = np.zeros(3)
        bounds = np.zeros(6)

        gridTransform.GetDisplacementGrid().GetCenter(center)
        gridTransform.GetDisplacementGrid().GetBounds(bounds)

        dims = np.zeros(3)
        for i in xrange(3):
            dims[i] = (bounds[i * 2 + 1] - bounds[i * 2]);

        # Create sample points in middle 75% of grid
        pointSrc = vtk.vtkPointSource()
        pointSrc.SetCenter(center)
        pointSrc.SetRadius(np.min((dims * 0.75) / 2))
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
        for i in xrange(toPoints.GetNumberOfPoints()):
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


