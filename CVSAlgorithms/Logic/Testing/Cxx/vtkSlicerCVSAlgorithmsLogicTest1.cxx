/*=auto=========================================================================
  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.
  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.
  Program:   3D Slicer
=========================================================================auto=*/

// Logic includes
#include "vtkSlicerCVSAlgorithmsLogic.h"

// MRML includes
#include "vtkMRMLScene.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLStorageNode.h"

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkImageReader.h>
#include <vtkNrrdReader.h>
#include <vtkImageData.h>

// STL includes
#include <fstream>
#include <limits>

// OS includes
#include <math.h>



//-----------------------------------------------------------------------------
int vtkSlicerCVSAlgorithmsLogicTest1(int argc, char* argv [])
{

  if (argc < 2)
  {
    std::cerr << "Missing volume file name." << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkMRMLScene> scene;

  vtkNew<vtkSlicerCVSAlgorithmsLogic> moduleLogic;
  moduleLogic->SetMRMLScene(scene);
  if (moduleLogic->GetMRMLScene() != scene)
  {
    std::cerr << "A MRML Scene must be set to go further." << std::endl;
    return EXIT_FAILURE;
  }
	
  vtkNew<vtkMRMLScalarVolumeNode> volumeNode;
  vtkNew<vtkImageData> imageData;
  imageData->SetDimensions(40, 40, 40);
  imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  unsigned char* ptr = reinterpret_cast<unsigned char*>(
	  imageData->GetScalarPointer(0, 0, 0));
  for (int z = 0; z < 40; ++z)
  {
	  for (int y = 0; y < 40; ++y)
	  {
		  for (int x = 0; x < 40; ++x)
		  {
			  double normalizedIntensity = (static_cast<double>(x + (y * 40) + (z * 40 * 40)) / static_cast<double>(40 * 40 * 40 - 1));
			  *(ptr++) = 255 - static_cast<unsigned char>(255. * normalizedIntensity);
		  }
	  }
  }

  volumeNode->SetAndObserveImageData(imageData.GetPointer());
  scene->AddNode(volumeNode.GetPointer());


  moduleLogic->testImage(volumeNode);

  return EXIT_SUCCESS;
}