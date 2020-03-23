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
#include <itkImageToVTKImageFilter.h>
#include <itkVTKImageToImageFilter.h>

// STD includes
#include <cassert>
#include <array>

#include <vnl\algo\>

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
	vtkOutputWindowDisplayText(volumeNode->GetName());

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
	myitkImage->Print(std::cout);
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
	for(auto seq : sequences)
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