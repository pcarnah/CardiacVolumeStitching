import gc
import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
import vtk.util.numpy_support as VN
from six.moves import range
import SimpleITK as sitk
import sitkUtils
import SimpleElastix
from timeit import default_timer as timer


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
        parameterSets = [['Par0001Rigid.txt'], ['Par0001NonRigid.txt'], ['Par0001Initialize.txt']]

        moduleDir = os.path.dirname(os.path.abspath(__file__))
        registrationResourcesDir = os.path.abspath(os.path.join(moduleDir, 'Resources', 'RegistrationParameters'))

        self._parameterMaps = [SimpleElastix.ReadParameterFile(os.path.join(registrationResourcesDir, parSet[0])) for
                               parSet in parameterSets]

        self._initializationPresets = [['Transesophageal', 'TEE'], ['Transgastric', 'TG'], ['Transthoracic', 'TTE']]

        # Create initialization transforms as identity
        self._initializationPresetTransforms = [None] * 3

        # TEE Transform
        tr = vtk.vtkTransform()
        tr.Identity()
        tr.RotateWXYZ(98.4210581181494, 0.35740674433659325, -0.8628562094610169, -0.35740674433659314)

        tr = sitk.Similarity3DTransform()
        tr.SetRotation([-0.35740674433659325, 0.8628562094610169, -0.35740674433659314], np.deg2rad(98.4210581181494))
        # Rotates the en face view to the correct orientation with respect to transgastric views
        self._initializationPresetTransforms[InitializationPresets_TEE] = tr.GetInverse()

        # TG Transform
        # Uses identity
        tr = vtk.vtkTransform()
        tr.Identity()

        tr = sitk.Euler3DTransform()
        self._initializationPresetTransforms[InitializationPresets_TG] = tr

        # TODO TTE Transform
        tr = vtk.vtkTransform()
        tr.Identity()

        tr = sitk.Euler3DTransform()
        self._initializationPresetTransforms[InitializationPresets_TTE] = tr

    def getInitializationPresets(self):
        return self._initializationPresets

    def stitchSequences(self, sequences, presets, outputSequence, masterSequence):
        """
        Stitch a group of sequences together. First applies rigid registration using semi-simultaneous method,
        then for each frame in sequence applies deformable registration. Aligned volumes are re-sampled and merged
        to generate output sequence.
        :param sequences: List of Sequence nodes
        :param presets: TODO remove
        :param outputSequence: Sequence node to store final result
        :param masterSequence: Sequence node to act as master for frame rate
        :return: None
        """

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


        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            numberOfDataNodes = masterSequence.GetNumberOfDataNodes()

            seqBrowser.SetSelectedItemNumber(0)
            slicer.modules.sequencebrowser.logic().UpdateProxyNodesFromSequences(seqBrowser)

            sitkIms = [sitkUtils.PullVolumeFromSlicer(n) for n in proxyNodes]

            # Apply positional initialization and align centers for maximum overlap to begin automatic registration
            initTrs = [sitk.Transform(self._initializationPresetTransforms[presets[i]]) for i, _ in enumerate(sitkIms)]
            initTrs = self.initialAlignment(sitkIms, initialTrs=initTrs)


            # Perform semi-simultaneous rigid registration step
            rigidTrs = self.semiSimultaneousRegister(sitkIms, initialTrs=initTrs, numCycles=5)

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

                sitkIms = [sitkUtils.PullVolumeFromSlicer(n) for n in proxyNodes]

                finalTrs = self.semiSimultaneousRegister(sitkIms, initialTrs=rigidTrs,
                                                            numCycles=1, parMap=self._parameterMaps[1])

                outputImage = self.mergeVolumesSITK(sitkIms, finalTrs)
                sitkUtils.PushVolumeToSlicer(outputImage, outputVolume)

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

            # Create sequence browser for output if it does not already exist
            if not slicer.modules.sequencebrowser.logic().GetFirstBrowserNodeForSequenceNode(outputSequence):
                outputSequenceBrowser = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSequenceBrowserNode",
                                                                           outputSequence.GetName() + ' browser')
                slicer.modules.sequencebrowser.logic().AddSynchronizedNode(outputSequence, None, outputSequenceBrowser)

        return

    def semiSimultaneousRegister(self, volumeList, fixedVolume=0, initialTrs=None, numCycles=10, parMap=None):
        """
        Runs semi-simultaneous registration algorithm describe by Wachinger, C. et al., (2007). This method works by
        registering each moving volume to all other volumes in turn, for a set number of cycles.
        :param volumeList: List of volumes to group-register
        :param fixedVolume: Index of fixed volume in list
        :param initalsTrs: List of initial transforms
        :param numCycles: Number of cycles to run algorithm
        :param parMap: ParameterMap
        :return: List of transforms
        """

        if initialTrs:
            volumeList = [self.getResampledImage(im, tr) for im, tr in zip(volumeList, initialTrs)]

            trs = initialTrs[0:fixedVolume] + initialTrs[fixedVolume + 1:]

            # Make copies of transforms to avoid side effects
            trs = [sitk.Transform(tr) for tr in trs]
            fixedTr = sitk.Transform(initialTrs[fixedVolume])
        else:
            trs = [sitk.Transform() for _ in range(len(volumeList) - 1)]
            fixedTr = sitk.Transform()

        volumes = volumeList[0:fixedVolume] + volumeList[fixedVolume + 1:]
        fixed = volumeList[fixedVolume]
        masks = [self.generateMask(im) for im in volumes]
        fixedMask = self.generateMask(fixed)

        if not parMap:
            parMap = self._parameterMaps[0]

        parMap = SimpleElastix.ParameterMap(parMap)

        # Set fixed volume weight higher
        w = 3.0
        parMap["Metric0Weight"] = [str(w / (len(volumes) + w - 1))]
        for i in range(1, len(volumes)):
            parMap["Metric{}Weight".format(i)] = [str(1.0 / (len(volumes) + w - 1))]

        # Semi simultaneous registration algorithm (Wachinger, C. et al., 2007)
        for cycle in range(numCycles):
            print('Starting cycle {}'.format(cycle+1))
            slicer.app.processEvents(qt.QEventLoop.ExcludeUserInputEvents)

            for i, movingVolume in enumerate(volumes):
                fixedVolumes = [fixed] + volumes[0:i] + volumes[i+1:]
                fixedMasks = [fixedMask] + masks[0:i] + masks[i+1:]

                movingMask = masks[i]

                _, tr = self.registerVolumes(fixedVolumes, movingVolume, fixedMasks, movingMask, parMap)

                print('Registered volume {} in cycle {}'.format(i+1, cycle+1))

                volumes[i] = self.getResampledImage(movingVolume, tr)
                masks[i] = self.generateMask(volumes[i])

                trs[i].AddTransform(tr)

                slicer.app.processEvents(qt.QEventLoop.ExcludeUserInputEvents)

        finalTrs = trs[0:fixedVolume] + [fixedTr] + trs[fixedVolume:]
        for tr in finalTrs:
            tr.FlattenTransform()
            tr.MakeUnique()

        return finalTrs



    def registerVolumes(self, fixedVolumeList, movingVolume, fixedMaskList=(), movingMask=None,
                        parMap=None):
        """
        Runs registration from single moving image to 1 or more fixed images. Takes images, masks, and parameter map.
        :param fixedVolumeList: List of SimpleITK images for fixed images
        :param movingVolume: Single SimpleITK image for moving image
        :param fixedMaskList: List of SimpleITK images for fixed image masks
        :param movingMask: Single SimpleITK image for moving image mask
        :param parMap: ParameterMap specifying registration parameters
        :return: Re-sampled image, Transform
        """

        if not fixedVolumeList or not movingVolume:
            logging.debug('simultaneuousRegister: missing parameter')
            return

        if not isinstance(fixedVolumeList[0], sitk.Image) or not isinstance(movingVolume, sitk.Image):
            logging.debug('Expected volumes of type SimpleITK Image')
            return

        if not parMap:
            parMap = self._parameterMaps[0]

        fixedCount = len(fixedVolumeList)

        # Use copy constructor to get local copy of map
        parMap = SimpleElastix.ParameterMap(parMap)

        # Set parameters to match number of fixed images
        parMap['Registration'] = ['MultiMetricMultiResolutionRegistration']
        parMap['FixedImagePyramid'] = parMap['FixedImagePyramid'] * fixedCount
        parMap['Metric'] = parMap['Metric'] * fixedCount
        parMap['ImageSampler'] = parMap['ImageSampler'] * fixedCount

        # Create elastix filter and set parameters
        selx = SimpleElastix.ElastixImageFilter()
        selx.SetParameterMap(parMap)

        for im in fixedVolumeList:
            selx.AddFixedImage(im)

        for im in fixedMaskList:
            selx.AddFixedMask(im)

        selx.SetMovingImage(movingVolume)
        if movingMask:
            selx.SetMovingMask(movingMask)

        selx.LogToConsoleOn()
        selx.LogToFileOff()
        selx.SetOutputDirectory('')

        resultImage = selx.Execute()

        conv = SimpleElastix.TransformConverter()
        conv.SetParameterMap(selx.GetTransformParameterMap()[0])
        tr = conv.Execute()

        return resultImage, tr

    def initialAlignment(self, volumeList, fixedVolume=0, initialTrs=None):

        volumes = volumeList[0:fixedVolume] + volumeList[fixedVolume + 1:]
        fixed = volumeList[fixedVolume]
        masks = [self.generateMask(im, 0) for im in volumes]
        fixedMask = self.generateMask(fixed, 0)

        if initialTrs:
            trs = initialTrs[0:fixedVolume] + initialTrs[fixedVolume + 1:]

            # Make copies of transforms to avoid side effects
            trs = [sitk.Transform(tr) for tr in trs]
            fixedTr = sitk.Transform(initialTrs[fixedVolume])
        else:
            trs = [sitk.Euler3DTransform() for _ in volumes]
            fixedTr = sitk.Euler3DTransform()

        for i, im in enumerate(volumes):
            # Perform initial alignment using moments
            trs[i] = sitk.CenteredTransformInitializer(fixed,
                                                       im,
                                                       trs[i],
                                                       sitk.CenteredTransformInitializerFilter.GEOMETRY)

        finalTrs = trs[0:fixedVolume] + [fixedTr] + trs[fixedVolume:]
        for tr in finalTrs:
            tr.FlattenTransform()
            tr.MakeUnique()

        return finalTrs

    def mergeVolumesSITK(self, volumeList, transforms):

        # Get bounds of final volume
        bounds = np.zeros((len(volumeList), 6))
        for i, im in enumerate(volumeList):
            bounds[i, :] = self.getBounds(im, transforms[i], True)

        min = bounds.min(0)
        max = bounds.max(0)

        volumeBounds_ROI = np.array([min[0], max[1], min[2], max[3], min[4], max[5]])

        outputSpacing = 0.4

        roiDim = np.zeros(3)

        for i in range(3):
            roiDim[i] = (volumeBounds_ROI[i * 2 + 1] - volumeBounds_ROI[i * 2]) / outputSpacing

        roiDim = np.ceil(roiDim).astype('uint16')

        # Crete reference image for resampling
        refImage = sitk.Image(*roiDim.tolist(), sitk.sitkFloat32)
        refImage.SetOrigin([volumeBounds_ROI[1], volumeBounds_ROI[3], volumeBounds_ROI[4]])
        refImage.SetSpacing([outputSpacing] * 3)
        refImage.SetDirection([-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0])

        resampledIms = [None] * len(volumeList)
        structureIms = [None] * len(volumeList)

        alphaSumIm = sitk.Image(refImage)

        # for each volume, resample into output space
        for i, im in enumerate(volumeList):
            start = timer()
            resIm = self.getResampledImage(sitk.Cast(im, sitk.sitkFloat32), transforms[i], refImage, sitk.sitkLinear)

            end = timer()
            print("resample: {}".format(end - start))

            start = timer()
            mask = self.generateMask(resIm, 2)
            # mask = self.getResampledImage(mask, transforms[i], refImage, sitk.sitkNearestNeighbor)
            mask = sitk.Cast(mask, sitk.sitkFloat32)
            end = timer()
            print("generate mask: {}".format(end - start))

            # Get probe position for volume in IJK space
            probeOriginLPS = [(bounds[i, j * 2] + bounds[i, j * 2 + 1]) / 2 for j in range(3)]
            probeOriginLPS[2] = bounds[i, 4]

            probeOriginIJK = resIm.TransformPhysicalPointToIndex(probeOriginLPS)

            distanceSource = sitk.Image(refImage)
            distanceSource.MakeUnique()
            distanceSource[probeOriginIJK] = 1

            distance = sitk.SignedMaurerDistanceMap(sitk.Cast(distanceSource, sitk.sitkUInt8),
                                                    squaredDistance=False, useImageSpacing=True)
            distance = sitk.Cast(distance, sitk.sitkFloat32)

            minmax = sitk.MinimumMaximumImageFilter()
            minmax.Execute(distance)
            imageHighestIntensity = minmax.GetMaximum()
            imageLowestIntensity = minmax.GetMinimum()

            distance = ((imageHighestIntensity - distance) * (4.0 / imageHighestIntensity)) ** 4

            # Blur and threshold to identify structures
            structureIm = sitk.SmoothingRecursiveGaussian(resIm, 3)
            structureIm = sitk.BinaryThreshold(structureIm, 110, imageHighestIntensity, 20, 0)
            structureIm = sitk.Cast(structureIm, sitk.sitkFloat32)

            structureIm = (distance + structureIm) * mask

            resampledIms[i] = resIm
            structureIms[i] = structureIm

            refImage = refImage + (resIm * structureIm)
            alphaSumIm = alphaSumIm + structureIm

        refImage = refImage / alphaSumIm

        refImage = sitk.Cast(sitk.RescaleIntensity(refImage, 0, 255), sitk.sitkUInt8)

        return refImage

    def generateMask(self, image, erosionRadius=3):
        """
        Generates a mask using 1-max threshold to identify only areas with information from ultrasound volumes
        :param image: SimpleITK image
        :return: binary mask as SimpleITK image
        """
        minmax = sitk.MinimumMaximumImageFilter()
        minmax.Execute(image)
        imageHighestIntensity = minmax.GetMaximum()
        imageLowestIntensity = minmax.GetMinimum()

        thresh = sitk.BinaryThresholdImageFilter()
        thresh.SetLowerThreshold(1)
        thresh.SetUpperThreshold(imageHighestIntensity)
        mask = thresh.Execute(image)

        if erosionRadius:
            dist = sitk.SignedMaurerDistanceMap(mask, True, False, True)
            mask = sitk.BinaryThreshold(dist, erosionRadius, 1e6)

        return mask

    def getBounds(self, sitkImage, transform=None, inv=False, masked=False):
        """
        Gets the bounds of an image in LPS space (SimpleITK space). Takes optional transform.
        :param sitkImage: The input image
        :param tr: Optional transform
        :param inv: Invert transform?
        :param masked: Masked image to include only area with data?
        :return: bounds in form [minX, maxX, minY, maxY, minZ, maxZ]
        """

        if transform:
            tr = transform
        else:
            tr = sitk.Transform()

        if inv:
            invTr = None
            try:
                invTr = tr.GetInverse()
            except RuntimeError:
                pass

            if not invTr:
                # Create displacement field from transform

                # Use sparse displacement field as we only need approximate bounds
                outDirection = sitkImage.GetDirection()

                dir3 = np.array([outDirection[0], outDirection[4], outDirection[8]])
                outSpacing = [2, 2, 2]
                sizePhys = np.array(sitkImage.GetSize()) * sitkImage.GetSpacing()

                outOrigin = [outDirection[i] - 0.5 * sizePhys[i] * dir3[i] for i in range(3)]
                outSize = np.ceil((sizePhys / np.array(outSpacing)) * 2).astype('uint32').tolist()

                disFieldFilter = sitk.TransformToDisplacementFieldFilter()
                disFieldFilter.SetOutputDirection(outDirection)
                disFieldFilter.SetOutputOrigin(outOrigin)
                disFieldFilter.SetOutputSpacing(outSpacing)
                disFieldFilter.SetSize(outSize)
                dspf = disFieldFilter.Execute(tr)

                invTr = sitk.DisplacementFieldTransform(sitk.InvertDisplacementField(dspf))

            oldTr = tr
            tr = invTr

        if masked:
            mask = self.generateMask(sitkImage, 0)

            stats = sitk.LabelStatisticsImageFilter()
            stats.Execute(sitkImage, mask)

            boundsIJK = stats.GetBoundingBox(1)
            corners = [
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[0], boundsIJK[2], boundsIJK[4]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[0], boundsIJK[2], boundsIJK[5]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[0], boundsIJK[3], boundsIJK[4]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[0], boundsIJK[3], boundsIJK[5]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[1], boundsIJK[2], boundsIJK[4]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[1], boundsIJK[2], boundsIJK[5]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[1], boundsIJK[3], boundsIJK[4]])),
                tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([boundsIJK[1], boundsIJK[3], boundsIJK[5]]))
            ]
        else:
            dims = (np.array(sitkImage.GetSize()) - 1).tolist()

            # New bounds can be computed from transformed corners
            corners = [tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([0, 0, 0])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([0, 0, dims[2]])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([0, dims[1], dims[2]])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([dims[0], 0, dims[2]])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([dims[0], dims[1], dims[2]])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([0, dims[1], 0])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([dims[0], 0, 0])),
                       tr.TransformPoint(sitkImage.TransformIndexToPhysicalPoint([dims[0], dims[1], 0]))]

        min = np.min(corners, 0)
        max = np.max(corners, 0)

        bounds = [min[0], max[0],
                  min[1], max[1],
                  min[2], max[2]]

        return bounds

    def processLog(self, text):
        slicer.app.processEvents()

    def getResampledImage(self, sitkImage, tr=None, refImage=None, interpolator=sitk.sitkLinear):
        """
        Resamples an image with a given transform to its new location. Takes option spacing and interpolator parameters.
        Transform should be defined as transformation from final->initial.
        :param sitkImage: The image to resample
        :param tr: The transform to use for resampling
        :param spacing: Desired output spacing. Needs to be same length as image dimensions.
        :param interpolator: Resample interpolator to use
        :return: New SimpleITK Image resampled
        """

        if not tr:
            tr = sitk.Transform()

        resample = sitk.ResampleImageFilter()

        if not refImage:
            bounds = self.getBounds(sitkImage, tr, True)

            outOrigin = [bounds[1], bounds[3], bounds[4]]
            outDirections = [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0]
            outSpacing = sitkImage.GetSpacing()
            outSize = [0, 0, 0]
            for i in range(3):
                outSize[i] = int(np.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / outSpacing[i]))

            resample.SetOutputOrigin(outOrigin)
            resample.SetOutputDirection(outDirections)
            resample.SetOutputSpacing(outSpacing)
            resample.SetSize(outSize)
        else:
            resample.SetReferenceImage(refImage)

        resample.SetTransform(tr)
        resample.SetInterpolator(interpolator)
        resample.SetOutputPixelType(sitkImage.GetPixelID())

        return resample.Execute(sitkImage)






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
