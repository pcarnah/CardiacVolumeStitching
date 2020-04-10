/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// CVSAlgorithms Logic includes
#include "vtkSlicerCVSAlgorithmsLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLSequenceNode.h>
#include <vtkMRMLSequenceBrowserNode.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkSlicerSequenceBrowserLogic.h>
 
// ITK incudes
#include <itkImage.h>
#include <itkTransform.h>
#include <itkImageToVTKImageFilter.h>
#include <itkVTKImageToImageFilter.h>
#include <itkCommand.h>
#include <itkMultiResolutionGaussianSmoothingPyramidImageFilter.h>
#include <itkImageDuplicator.h>

// Elastix includes
#include "itkSemiSimultaneousImageRegistrationMethod.h"
#include "itkMultiMetricMultiResolutionImageRegistrationMethod.h"
#include "itkAdaptiveStochasticGradientDescentOptimizer.h"
#include "itkAdvancedImageToImageMetric.h"
#include "itkAdvancedNormalizedCorrelationImageToImageMetric.h"
#include "itkAdvancedEuler3DTransform.h"
#include "itkAdvancedLinearInterpolateImageFunction.h"
#include "itkImageRandomSampler.h"
#include "itkImageRandomCoordinateSampler.h"

// STD includes
#include <cassert>
#include <array>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerCVSAlgorithmsLogic);

//----------------------------------------------------------------------------
vtkSlicerCVSAlgorithmsLogic::vtkSlicerCVSAlgorithmsLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerCVSAlgorithmsLogic::~vtkSlicerCVSAlgorithmsLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::testImage(vtkMRMLNode* node)
{
	vtkSmartPointer<vtkMRMLScalarVolumeNode> volumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(node);

	constexpr unsigned int Dimension = 3;
	using PixelType = unsigned char;
	using ImageType = itk::Image< PixelType, Dimension >;

	using FilterType = itk::VTKImageToImageFilter< ImageType >;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( volumeNode->GetImageData());

	try
	{
		filter->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		vtkErrorMacro("Error: " << error << std::endl);
		return;
	}

	ImageType::ConstPointer myitkImage = filter->GetOutput();
	
	//Registration setup
	//TODO Make factory for this? Possibly python wrap for easier testing
	//Have this class wrap the registration, pass on Set, Get image methods and construct all relevent things
	
	using DuplicatorType = itk::ImageDuplicator<ImageType>;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(myitkImage);
	duplicator->Modified();

	//Declare types
	//using RegistrationType = itk::SemiSimultaneousImageRegistrationMethod<ImageType, ImageType>;
	using RegistrationType = itk::SemiSimultaneousImageRegistrationMethod<ImageType, ImageType>;
	using OptimizerType = itk::AdaptiveStochasticGradientDescentOptimizer;
	using BaseMetricType = itk::AdvancedImageToImageMetric<ImageType, ImageType>;
	using MetricType = itk::AdvancedNormalizedCorrelationImageToImageMetric<ImageType, ImageType>;
	using CombinationMetricType = itk::CombinationImageToImageMetric<ImageType, ImageType>;
	using ImagePyramidType = itk::MultiResolutionGaussianSmoothingPyramidImageFilter< ImageType, ImageType >;
	using InterpolatorType = itk::AdvancedLinearInterpolateImageFunction<ImageType, MetricType::CoordinateRepresentationType>;
	using ImageSamplerType = itk::ImageRandomSampler<ImageType>;
	using TransformType = itk::AdvancedEuler3DTransform<RegistrationType::CombinationMetricType::TransformType::ScalarType>;

	//Create registration method
	RegistrationType::Pointer reg = RegistrationType::New();
	reg->SetNumberOfFixedImages(1);
	reg->SetNumberOfFixedImagePyramids(1);
	reg->SetNumberOfFixedImageRegions(1);
	reg->SetNumberOfMovingImages(1);
	reg->SetNumberOfMovingImagePyramids(1);
	reg->SetNumberOfLevels(4);
	reg->SetGlobalIterations(1);
	reg->SetNumberOfInterpolators(1);

	//Add fixed image
	reg->SetFixedImage(myitkImage);
	reg->SetFixedImagePyramid(ImagePyramidType::New());
	reg->SetFixedImageRegion(myitkImage->GetLargestPossibleRegion());
	reg->SetTransform(TransformType::New());

	OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetNumberOfIterations(50);

	itk::Array<double> scales(6);
	scales[0] = 5000;
	scales[1] = 5000;
	scales[2] = 5000;
	scales[3] = 1;
	scales[4] = 1;
	scales[5] = 1;

	optimizer->SetScales(scales);
	optimizer->SetLearningRate(20);

	CombinationMetricType::Pointer combinationMetric = CombinationMetricType::New();
	combinationMetric->SetNumberOfMetrics(1);
	combinationMetric->SetUseAllMetrics();
	combinationMetric->SetMetricWeight(1, 0);

	
	for (int i = 0; i < 1; i++) {
		MetricType::Pointer m = MetricType::New();
		m->SetImageSampler(ImageSamplerType::New());
		m->SetRequiredRatioOfValidSamples(0.05);
		m->GetImageSampler()->SetNumberOfSamples(30000);

		combinationMetric->SetMetric(m, i);

		duplicator->Modified();
		duplicator->Update();
		double origin[] = { 5.0, -3.7, 2.0 };
		duplicator->GetOutput()->SetOrigin(origin);

		reg->SetMovingImage(duplicator->GetOutput(), i);
		reg->SetMovingImagePyramid(ImagePyramidType::New(), i);
		reg->SetInterpolator(InterpolatorType::New(), i);
	}

	reg->SetInitialTransformParameters(reg->GetTransform()->GetParameters());

	reg->SetMetric(combinationMetric);
	reg->SetOptimizer(optimizer);
	
	try {
		reg->Update();
	}
	catch (itk::ExceptionObject &err) {
		std::cerr << err.GetDescription() << std::endl;
	}

}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::testWrapping()
{
	vtkOutputWindowDisplayText("Hello");
}

//---------------------------------------------------------------------------
void vtkSlicerCVSAlgorithmsLogic::stitchSequences(std::vector<std::string> sequencesBase, vtkMRMLNode* outputSequenceBase, vtkMRMLNode* masterSequenceBase)
{
	assert(this->GetMRMLScene() != 0);

	vtkNew<vtkSlicerSequenceBrowserLogic> seqLogic;

	std::vector<vtkSmartPointer<vtkMRMLNode>> sequences;
	for (auto seq : sequencesBase)
	{
		vtkSmartPointer<vtkMRMLNode> n = this->GetMRMLScene()->GetNodeByID(seq);
		if (n)
		{
			sequences.push_back(n);
		}
	}

	vtkSmartPointer<vtkMRMLSequenceNode> outputSequence = vtkMRMLSequenceNode::SafeDownCast(outputSequenceBase);
	if (!outputSequence)
	{
		vtkErrorMacro("Invalid output sequence.");
		return;
	}

	vtkSmartPointer<vtkMRMLSequenceNode> masterSequence = vtkMRMLSequenceNode::SafeDownCast(masterSequenceBase);
	if (!masterSequence)
	{
		vtkErrorMacro("Invalid master sequence.");
		return;
	}

	// Clear output sequence
	outputSequence->RemoveAllDataNodes();
	outputSequence->SetIndexType(masterSequence->GetIndexType());
	outputSequence->SetIndexName(masterSequence->GetIndexName());
	outputSequence->SetIndexUnit(masterSequence->GetIndexUnit());

	// Set up browser
	vtkNew<vtkMRMLSequenceBrowserNode> seqBrowser;
	this->GetMRMLScene()->AddNode(seqBrowser);
	bool masterNodeFound = false;
	for (auto seq : sequences)
	{
		seqLogic->AddSynchronizedNode(seq, nullptr, seqBrowser);
		if (seq->GetID() == masterSequence->GetID())
		{
			masterNodeFound = true;
		}

		if (!masterNodeFound)
		{
			seqLogic->AddSynchronizedNode(masterSequence, nullptr, seqBrowser);
		}
		seqBrowser->SetAndObserveMasterSequenceNodeID(masterSequence->GetID());

		// Create temporary volume for output
		vtkNew<vtkMRMLScalarVolumeNode> outputVolume;
		this->GetMRMLScene()->AddNode(outputVolume);

		// Get proxy nodes and generate mask
		std::vector<vtkSmartPointer<vtkMRMLScalarVolumeNode>> proxyNodes;
		for (auto n : sequences)
		{
			vtkSmartPointer<vtkMRMLScalarVolumeNode> node = vtkMRMLScalarVolumeNode::SafeDownCast(n);
			if (n)
			{
				proxyNodes.push_back(node);
			}
		}

		const int n = masterSequence->GetNumberOfDataNodes();

		seqBrowser->SetSelectedItemNumber(0);
		seqLogic->UpdateProxyNodesFromSequences(seqBrowser);
	}
}