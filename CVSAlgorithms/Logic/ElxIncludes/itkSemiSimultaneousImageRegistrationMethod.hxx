#include "itkSemiSimultaneousImageRegistrationMethod.h"

#include <itkTimeProbe.h>
#include <itkResampleImageFilter.h>

/** Macro that implements the set methods. */
#define itkImplementationSetMacro( _name, _type ) \
  template< typename TFixedImage, typename TMovingImage > \
  void \
  SemiSimultaneousImageRegistrationMethod < TFixedImage, TMovingImage > \
  ::Set##_name( _type _arg, unsigned int pos ) \
  { \
    if( pos == 0 ) \
    { \
      this->Superclass::Set##_name( _arg ); \
    } \
    if( pos >= this->GetNumberOf##_name##s() ) \
    { \
      this->SetNumberOf##_name##s( pos + 1 ); \
    } \
    if( this->m_##_name##s[ pos ] != _arg ) \
    { \
      this->m_##_name##s[ pos ] = _arg; \
      this->Modified(); \
    } \
  } // comment to allow ; after calling macro

/** Macro that implements the get methods. */
#define itkImplementationGetMacro( _name, _type1, _type2 ) \
  template< typename TFixedImage, typename TMovingImage > \
  _type1 typename \
  SemiSimultaneousImageRegistrationMethod < TFixedImage, TMovingImage > \
  ::_type2 \
  SemiSimultaneousImageRegistrationMethod < TFixedImage, TMovingImage > \
  ::Get##_name( unsigned int pos ) const \
  { \
    if( pos >= this->GetNumberOf##_name##s() ) \
    { \
      return 0; \
    } \
    else \
    { \
      return this->m_##_name##s[ pos ].GetPointer(); \
    } \
  } // comment to allow ; after calling macro

namespace itk {


	itkImplementationSetMacro(Transform, TransformType *);

	itkImplementationGetMacro(Transform, , TransformType *);

	/**
	 * ****************** Constructor ******************
	 */

	template <typename TFixedImage, typename TMovingImage>
	SemiSimultaneousImageRegistrationMethod<TFixedImage, TMovingImage>::SemiSimultaneousImageRegistrationMethod()
		: MultiMetricMultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage>()
	{
		this->m_GlobalIterations = 1;

		this->m_AfterEachIterationCommand = AfterEachIterationCommandType::New();
		this->m_AfterEachIterationCommand->SetCallbackFunction( this, &Self::AfterEachIteration );
	} // end Constructor

	/**
	 * ****************** Initialize *******************************
	 */

	template <typename TFixedImage, typename TMovingImage>
	void SemiSimultaneousImageRegistrationMethod <TFixedImage, TMovingImage>::Initialize(void)
	{
		this->CheckOnInitialize();

		//using ResampleImageType = itk::ResampleImageFilter<TMovingImage, TMovingImage>;
		//ResampleImageType::Pointer resample = ResampleImageType::New();
		//resample->SetInterpolator(this->GetInterpolator());
		//resample->SetSize(this->GetFixedImageRegion().GetSize());
		//resample->SetInput(this->GetFixedImagePyramid()->GetOutput(this->GetCurrentLevel()));
		//
		//LinearTransformPointer inverse = LinearTransformType::New();
		//LinearTransformPointer transform = dynamic_cast<LinearTransformType *>(this->GetTransform(this->GetCurrentVolume()));

		//if (transform.IsNotNull()) {
		//	bool success = transform->GetInverse(inverse.GetPointer());
		//	std::cout << "Invertible: " << success << std::endl;

		//	resample->SetTransform(inverse);
		//	resample->Update();
		//}

		
		/** Setup the metric. */
		this->GetCombinationMetric()->SetTransform(this->GetModifiableTransform(this->GetCurrentVolume()));

		this->GetCombinationMetric()->SetFixedImage(
			this->GetFixedImagePyramid()->GetOutput(this->GetCurrentLevel()));
		for (unsigned int i = 0; i < this->GetNumberOfFixedImagePyramids(); ++i)
		{
			this->GetCombinationMetric()->SetFixedImage(
				this->GetFixedImagePyramid(i)->GetOutput(this->GetCurrentLevel()), i);
		}

		this->GetCombinationMetric()->SetMovingImage(
			this->GetMovingImagePyramid(this->GetCurrentVolume())->GetOutput(this->GetCurrentLevel()));
		for (unsigned int i = 0; i < this->GetNumberOfMovingImagePyramids(); ++i)
		{
			if (i != this->GetCurrentVolume())
			{
				this->GetCombinationMetric()->SetMovingImage(
					this->GetMovingImagePyramid(i)->GetOutput(this->GetCurrentLevel()), i);
			}
		}

		this->GetCombinationMetric()->SetInterpolator(this->GetInterpolator());
		for (unsigned int i = 0; i < this->GetNumberOfInterpolators(); ++i)
		{
			this->GetCombinationMetric()->SetInterpolator(this->GetInterpolator(i), i);
		}

		this->GetCombinationMetric()->SetFixedImageRegion(
			this->m_FixedImageRegionPyramids[0][this->GetCurrentLevel()]);
		for (unsigned int i = 0; i < this->m_FixedImageRegionPyramids.size(); ++i)
		{
			this->GetCombinationMetric()->SetFixedImageRegion(
				this->m_FixedImageRegionPyramids[i][this->GetCurrentLevel()], i);
		}

		//this->GetMetric()->Initialize();
		this->GetCombinationMetric()->Initialize();

		/** Setup the optimizer. */
		this->GetModifiableOptimizer()->SetCostFunction(this->GetModifiableMetric());
		this->GetModifiableOptimizer()->SetInitialPosition(
			this->GetTransform(this->GetCurrentVolume())->GetParameters());

		cout << this->GetInitialTransformParametersOfNextLevel() << endl;
		cout << this->GetTransform()->GetParameters() << endl << endl;

		/** Connect the transform to the Decorator. */
		TransformOutputType *transformOutput = static_cast<TransformOutputType *>(this->ProcessObject::GetOutput(0));

		transformOutput->Set(this->GetTransform());

	} // end Initialize()

	/**
	 * ********************* GenerateData ***********************
	 */

	template <typename TFixedImage, typename TMovingImage>
	void SemiSimultaneousImageRegistrationMethod<TFixedImage, TMovingImage>::GenerateData(void)
	{
		itk::TimeProbe timer;
		timer.Start();

		this->m_Stop = false;

		unsigned long afterIterationID = this->GetModifiableOptimizer()->AddObserver(itk::IterationEvent(), this->m_AfterEachIterationCommand);

		/** Check the transform and set the initial parameters. */
		if (this->GetTransform() == 0)
		{
			itkExceptionMacro(<< "Transform is not present");
		}

		this->SetInitialTransformParametersOfNextLevel(
			this->GetInitialTransformParameters());

		if (this->GetInitialTransformParametersOfNextLevel().Size() != this->GetTransform()->GetNumberOfParameters())
		{
			itkExceptionMacro(<< "Size mismatch between initial parameter and transform");
		}

		/** Prepare the fixed and moving pyramids. */
		this->PrepareAllPyramids();

		/** Loop over the resolution levels. */
		for (unsigned int currentLevel = 0; currentLevel < this->GetNumberOfLevels();
			currentLevel++)
		{
			this->SetCurrentLevel(currentLevel);

			// Invoke an iteration event.
			// This allows a UI to reset any of the components between
			// resolution level.
			this->InvokeEvent(IterationEvent());

			// Check if there has been a stop request
			if (this->m_Stop)
			{
				break;
			}

			for (unsigned int cycle = 0; cycle < this->m_GlobalIterations; cycle++)
			{
				//Set current cycle
				this->SetCurrentCycle(cycle);

				for (unsigned int i = 0; i < this->GetNumberOfMovingImages(); i++)
				{
					//Set current moving image
					this->SetCurrentVolume(i);

					try
					{
						// initialize the interconnects between components
						this->Initialize();
					}
					catch (ExceptionObject &err)
					{
						this->m_LastTransformParameters = ParametersType(1);
						this->m_LastTransformParameters.Fill(0.0f);

						// pass exception to caller
						throw err;
					}

					try
					{
						// do the optimization
						this->GetModifiableOptimizer()->StartOptimization();
					}
					catch (ExceptionObject &err)
					{
						// An error has occurred in the optimization.
						// Update the parameters
						this->m_LastTransformParameters = this->GetOptimizer()->GetCurrentPosition();

						// Pass exception to caller
						throw err;
					}

					// get the results
					this->m_LastTransformParameters = this->GetOptimizer()->GetCurrentPosition();
					this->GetModifiableTransform()->SetParameters(this->m_LastTransformParameters);

					// setup the initial parameters for next level
					if (this->GetCurrentLevel() < this->GetNumberOfLevels() - 1)
					{
						this->SetInitialTransformParametersOfNextLevel(
							this->m_LastTransformParameters);
					}
				} // end for loop over moving images

			} // end for loop over cycles

		} // end for loop over res levels

		this->GetModifiableOptimizer()->RemoveObserver(afterIterationID);

		timer.Stop();
		std::cout << "Registration took: " << timer.GetMean() << std::endl;

	} // end GenerateData()

	/**
	 * ****************** CheckPyramids ******************
	 */

	template< typename TFixedImage, typename TMovingImage >
	void SemiSimultaneousImageRegistrationMethod< TFixedImage, TMovingImage >
		::CheckPyramids(void)
	{
		/** Check if at least one of the following are provided. */
		if (this->GetFixedImage() == 0)
		{
			itkExceptionMacro(<< "FixedImage is not present");
		}
		if (this->GetMovingImage() == 0)
		{
			itkExceptionMacro(<< "MovingImage is not present");
		}
		if (this->GetFixedImagePyramid() == 0)
		{
			itkExceptionMacro(<< "Fixed image pyramid is not present");
		}
		if (this->GetMovingImagePyramid() == 0)
		{
			itkExceptionMacro(<< "Moving image pyramid is not present");
		}

		/** Check if the number if fixed/moving pyramids >= nr of fixed/moving images,
		 * and whether the number of fixed image regions == the number of fixed images.
		 */
		if (this->GetNumberOfFixedImagePyramids() < this->GetNumberOfFixedImages())
		{
			itkExceptionMacro(<< "The number of fixed image pyramids should be >= "
				<< "the number of fixed images");
		}
		if (this->GetNumberOfMovingImagePyramids() < this->GetNumberOfMovingImages())
		{
			itkExceptionMacro(<< "The number of moving image pyramids should be >= "
				<< "the number of moving images");
		}
		if (this->GetNumberOfFixedImageRegions() != this->GetNumberOfFixedImages())
		{
			itkExceptionMacro(<< "The number of fixed image regions should equal "
				<< "the number of fixed images");
		}

	} // end CheckPyramids()


	/**
	 * ****************** CheckOnInitialize ******************
	 */

	template< typename TFixedImage, typename TMovingImage >
	void
		SemiSimultaneousImageRegistrationMethod< TFixedImage, TMovingImage >
		::CheckOnInitialize(void)
	{
		/** Check if at least one of the following is present. */
		if (this->GetMetric() == 0)
		{
			itkExceptionMacro(<< "Metric is not present");
		}
		if (this->GetOptimizer() == 0)
		{
			itkExceptionMacro(<< "Optimizer is not present");
		}
		if (this->GetTransform() == 0)
		{
			itkExceptionMacro(<< "Transform is not present");
		}
		if (this->GetInterpolator() == 0)
		{
			itkExceptionMacro(<< "Interpolator is not present");
		}
		if (this->GetModifiableTransform()->GetTransformCategory() != TransformCategoryType::Linear)
		{
			itkExceptionMacro(<< "Transform must be linear");
		}

		/** nrofmetrics >= nrofinterpolators >= nrofpyramids >= nofimages */
		unsigned int nrOfMetrics = this->GetCombinationMetric()->GetNumberOfMetrics();
		if (this->GetNumberOfInterpolators() > nrOfMetrics)
		{
			itkExceptionMacro(<< "NumberOfInterpolators can not exceed the "
				<< "NumberOfMetrics in the CombinationMetric!");
		}
		if (this->GetNumberOfFixedImagePyramids() > nrOfMetrics)
		{
			itkExceptionMacro(<< "NumberOfFixedImagePyramids can not exceed the "
				<< "NumberOfMetrics in the CombinationMetric!");
		}
		if (this->GetNumberOfMovingImagePyramids() > nrOfMetrics)
		{
			itkExceptionMacro(<< "NumberOfMovingImagePyramids can not exceed the "
				<< "NumberOfMetrics in the CombinationMetric!");
		}
		if (this->GetNumberOfMovingImagePyramids() >
			this->GetNumberOfInterpolators())
		{
			itkExceptionMacro(<< "NumberOfMovingImagePyramids can not exceed the "
				<< "NumberOfInterpolators!");
		}

		/** For all components: ==nrofmetrics of ==1. */
		if ((this->GetNumberOfInterpolators() != 1)
			&& (this->GetNumberOfInterpolators() != nrOfMetrics))
		{
			itkExceptionMacro(<< "The NumberOfInterpolators should equal 1 "
				<< "or equal the NumberOfMetrics");
		}
		if ((this->GetNumberOfFixedImagePyramids() != 1)
			&& (this->GetNumberOfFixedImagePyramids() != nrOfMetrics))
		{
			itkExceptionMacro(<< "The NumberOfFixedImagePyramids should equal 1 "
				<< "or equal the NumberOfMetrics");
		}
		if ((this->GetNumberOfMovingImagePyramids() != 1)
			&& (this->GetNumberOfMovingImagePyramids() != nrOfMetrics))
		{
			itkExceptionMacro(<< "The NumberOfMovingImagePyramids should equal 1 "
				<< "or equal the NumberOfMetrics");
		}
		if ((this->GetNumberOfFixedImages() != 1)
			&& (this->GetNumberOfFixedImages() != nrOfMetrics))
		{
			itkExceptionMacro(<< "The NumberOfFixedImages should equal 1 "
				<< "or equal the NumberOfMetrics");
		}
		if ((this->GetNumberOfMovingImages() != 1)
			&& (this->GetNumberOfMovingImages() != nrOfMetrics))
		{
			itkExceptionMacro(<< "The NumberOfMovingImages should equal 1 "
				<< "or equal the NumberOfMetrics");
		}

	} // end CheckOnInitialize()

	template< typename TFixedImage, typename TMovingImage >
	void
		SemiSimultaneousImageRegistrationMethod< TFixedImage, TMovingImage >
		::AfterEachIteration(void) {
		
		for (unsigned int i = 0; i < this->GetCombinationMetric()->GetNumberOfMetrics(); i++) {
			MetricPointer metric = dynamic_cast<MetricType *>(this->GetCombinationMetric()->GetMetric(i));
			if (metric.IsNotNull() && metric->GetUseImageSampler()) {
				metric->GetImageSampler()->SelectNewSamplesOnUpdate();
			}
		}
	}
}