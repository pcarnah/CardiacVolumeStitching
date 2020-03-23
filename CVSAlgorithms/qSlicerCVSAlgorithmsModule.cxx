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
#include <vtkSlicerCVSAlgorithmsLogic.h>

// CVSAlgorithms includes
#include "qSlicerCVSAlgorithmsModule.h"
#include "qSlicerCVSAlgorithmsModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerCVSAlgorithmsModulePrivate
{
public:
  qSlicerCVSAlgorithmsModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerCVSAlgorithmsModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerCVSAlgorithmsModulePrivate::qSlicerCVSAlgorithmsModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerCVSAlgorithmsModule methods

//-----------------------------------------------------------------------------
qSlicerCVSAlgorithmsModule::qSlicerCVSAlgorithmsModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerCVSAlgorithmsModulePrivate)
{
	this->setWidgetRepresentationCreationEnabled(false);
}

//-----------------------------------------------------------------------------
qSlicerCVSAlgorithmsModule::~qSlicerCVSAlgorithmsModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerCVSAlgorithmsModule::helpText() const
{
  return "This is a hidden loadable module that contains logic for cardiac volume stitching";
}

//-----------------------------------------------------------------------------
QString qSlicerCVSAlgorithmsModule::acknowledgementText() const
{
  return "This work developed by Patrick Carnahan (Robarts Research Institute)";
}

//-----------------------------------------------------------------------------
QStringList qSlicerCVSAlgorithmsModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Patrick Carnahan (Robarts Research Institute)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerCVSAlgorithmsModule::icon() const
{
  return QIcon(":/Icons/CVSAlgorithms.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerCVSAlgorithmsModule::categories() const
{
  return QStringList() << "Cardiac";
}

//-----------------------------------------------------------------------------
QStringList qSlicerCVSAlgorithmsModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerCVSAlgorithmsModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerCVSAlgorithmsModule
::createWidgetRepresentation()
{
  return new qSlicerCVSAlgorithmsModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerCVSAlgorithmsModule::createLogic()
{
  return vtkSlicerCVSAlgorithmsLogic::New();
}
