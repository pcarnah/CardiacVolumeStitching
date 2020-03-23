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

// .NAME vtkSlicerCVSAlgorithmsLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerCVSAlgorithmsLogic_h
#define __vtkSlicerCVSAlgorithmsLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// VTK includes
#include <vtkDenseArray.h>

// STD includes
#include <cstdlib>

#include "vtkSlicerCVSAlgorithmsModuleLogicExport.h"


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_CVSALGORITHMS_MODULE_LOGIC_EXPORT vtkSlicerCVSAlgorithmsLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerCVSAlgorithmsLogic *New();
  vtkTypeMacro(vtkSlicerCVSAlgorithmsLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  void testImage(vtkMRMLNode* node);
  void testWrapping();
  void stitchSequences(std::vector<std::string> sequences, vtkMRMLNode* outputSequence, vtkMRMLNode* masterSequence);

protected:
  vtkSlicerCVSAlgorithmsLogic();
  virtual ~vtkSlicerCVSAlgorithmsLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);

private:

  vtkSlicerCVSAlgorithmsLogic(const vtkSlicerCVSAlgorithmsLogic&); // Not implemented
  void operator=(const vtkSlicerCVSAlgorithmsLogic&); // Not implemented
};

#endif