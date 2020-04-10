#ifndef __itkSemiSimultaneousImageRegistrationMethod_h
#define __itkSemiSimultaneousImageRegistrationMethod_h

#include "itkMultiMetricMultiResolutionImageRegistrationMethod.h"

#include "itkCommand.h"
#include "itkAdvancedMatrixOffsetTransformBase.h"
#include <itkTransformBase.h>

/** defines a method that calls the same method
 * with an extra 0 argument */
#define elxOverrideSimpleSetMacro( _name, _type ) \
  void Set##_name( _type _arg ) override \
  { \
    this->Set##_name( _arg, 0 ); \
  }

 /** defines for example: SetNumberOfInterpolators() */
#define itkSetNumberOfMacro( _name ) \
  virtual void SetNumberOf##_name##s( unsigned int _arg ) \
  { \
    if( this->m_##_name##s.size() != _arg ) \
    { \
      this->m_##_name##s.resize( _arg ); \
      this->Modified(); \
    } \
  }

/** defines for example: GetNumberOfInterpolators() */
#define itkGetNumberOfMacro( _name ) \
  virtual unsigned int GetNumberOf##_name##s( void ) const \
  { \
    return this->m_##_name##s.size(); \
  }

namespace itk
{

	/** \class SemiSimultaneousImageRegistrationMethod
	 * \brief Class for multi-resolution semi-simulataneous image registration method
	 *
	 * This class is an extension of the itk class
	 * MultiMetricMultiResolutionImageRegistrationMethod. It allows the use
	 * of multiple metrics, which are summed, multiple images,
	 * multiple interpolators, and/or multiple image pyramids.
	 *
	 * Make sure the following is true:\n
	 *   nrofmetrics >= nrofinterpolators >= nrofmovingpyramids >= nrofmovingimages\n
	 *   nrofmetrics >= nroffixedpyramids >= nroffixedimages\n
	 *   nroffixedregions == nroffixedimages\n
	 *
	 *   nrofinterpolators == nrofmetrics OR nrofinterpolators == 1\n
	 *   nroffixedimages == nrofmetrics OR nroffixedimages == 1\n
	 *   etc...
	 *
	 * You may also set an interpolator/fixedimage/etc to NULL, if you
	 * happen to know that the corresponding metric is not an
	 * ImageToImageMetric, but a regularizer for example (which does
	 * not need an image.
	 *
	 *
	 * \sa ImageRegistrationMethod
	 * \sa MultiResolutionImageRegistrationMethod
	 * \ingroup RegistrationFilters
	 */

	template< typename TFixedImage, typename TMovingImage >
	class SemiSimultaneousImageRegistrationMethod :
		public MultiMetricMultiResolutionImageRegistrationMethod< TFixedImage, TMovingImage >
	{
	public:
		/** Standard class typedefs. */
		typedef SemiSimultaneousImageRegistrationMethod Self;
		typedef MultiMetricMultiResolutionImageRegistrationMethod<
			TFixedImage, TMovingImage >                               Superclass;
		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(MultiMetricMultiResolutionImageRegistrationMethod,
			MultiResolutionImageRegistrationMethod2);

		/**  Superclass types */
		typedef typename Superclass::FixedImageType          FixedImageType;
		typedef typename Superclass::FixedImageConstPointer  FixedImageConstPointer;
		typedef typename Superclass::FixedImageRegionType    FixedImageRegionType;
		typedef typename Superclass::MovingImageType         MovingImageType;
		typedef typename Superclass::MovingImageConstPointer MovingImageConstPointer;

		typedef typename Superclass::MetricType               MetricType;
		typedef typename Superclass::MetricPointer            MetricPointer;
		typedef typename Superclass::TransformType            TransformType;
		typedef typename Superclass::TransformPointer         TransformPointer;
		typedef typename Superclass::InterpolatorType         InterpolatorType;
		typedef typename Superclass::InterpolatorPointer      InterpolatorPointer;
		typedef typename Superclass::OptimizerType            OptimizerType;
		typedef typename OptimizerType::Pointer               OptimizerPointer;
		typedef typename Superclass::FixedImagePyramidType    FixedImagePyramidType;
		typedef typename Superclass::FixedImagePyramidPointer FixedImagePyramidPointer;
		typedef typename Superclass::MovingImagePyramidType   MovingImagePyramidType;
		typedef typename
			Superclass::MovingImagePyramidPointer MovingImagePyramidPointer;

		typedef typename Superclass::TransformOutputType    TransformOutputType;
		typedef typename Superclass::TransformOutputPointer TransformOutputPointer;
		typedef typename
			Superclass::TransformOutputConstPointer TransformOutputConstPointer;

		typedef typename Superclass::ParametersType    ParametersType;
		typedef typename Superclass::DataObjectPointer DataObjectPointer;

		typedef itk::SimpleMemberCommand< Self >                  AfterEachIterationCommandType;
		typedef typename AfterEachIterationCommandType::Pointer   AfterEachIterationCommandPointer;

		typedef itk::AdvancedMatrixOffsetTransformBase<double, 3U, 3U>		LinearTransformType;		
		typedef typename LinearTransformType::Pointer						LinearTransformPointer;

		typedef itk::TransformBase::TransformCategoryType			TransformCategoryType;

		/** Set/Get the Transform. */
		itkSetNumberOfMacro(Transform);
		itkGetNumberOfMacro(Transform);

		virtual void SetTransform(TransformType * _arg, unsigned int pos);

		virtual TransformType * GetTransform(unsigned int pos) const;
		virtual TransformType * GetModifiableTransform(unsigned int pos)
		{
			if (pos >= this->GetNumberOfTransforms()) {
				return 0;
			}
			else {
				return this->m_Transforms[pos].GetPointer();
			}
		}

		TransformType * GetTransform(void) override
		{
			return this->GetTransform(0);
		}
		TransformType * GetModifiableTransform(void) override
		{
			return this->GetModifiableTransform(0);
		}
		elxOverrideSimpleSetMacro(Transform, TransformType *);


		itkSetMacro(GlobalIterations, unsigned int);
		itkGetMacro(GlobalIterations, unsigned int);

		itkSetMacro(CurrentCycle, unsigned int);
		itkGetMacro(CurrentCycle, unsigned int);

		itkSetMacro(CurrentVolume, unsigned int);
		itkGetMacro(CurrentVolume, unsigned int);

	protected:

		SemiSimultaneousImageRegistrationMethod();
		~SemiSimultaneousImageRegistrationMethod() override {}

		/** Initialize by setting the interconnects between the components.
		 * This method is executed at every level of the pyramid with the
		 * values corresponding to this resolution.
		 */
		void Initialize(void) override;

		/** Method invoked by the pipeline in order to trigger the computation of
		   * the registration.
		   */
		void GenerateData(void) override;

		/** Function called by PrepareAllPyramids, which checks if the user input
		* regarding the image pyramids is ok.
		*/
		virtual void CheckPyramids(void);

		/** Function called by Initialize, which checks if the user input
		* is ok. Called by Initialize().
		*/
		virtual void CheckOnInitialize(void);

		void AfterEachIteration(void);


		std::vector< TransformPointer >       m_Transforms;

		unsigned int m_GlobalIterations;
		unsigned int m_CurrentCycle;
		unsigned int m_CurrentVolume;

		/** CallBack commands. */
  		AfterEachIterationCommandPointer   m_AfterEachIterationCommand;
	};

}

#undef itkSetNumberOfMacro
#undef itkGetNumberOfMacro


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSemiSimultaneousImageRegistrationMethod.hxx"
#endif

#endif