#ifndef __metrics_cpp
#define __metrics_cpp

#include "Metrics.h"

Metrics::Metrics(void)
{
	val = 0.0;

	transform = TransformBType::New();
	interpolator = InterpolatorType::New();

	// Dice Overlap index.
	th1 = ThresholdType::New();
	th2 = ThresholdType::New();
	overlap = OverlapMeasuresType::New();

	// Sum of Squares difference.
	SQmetric = MSQType::New();
	// Mututal Information
	MImetric = MIType::New();
	Mattesmetric =MattesMIType::New();

	// Normalized Cross-Correlation;
	NCCmetric = NCCType::New();

	// Gradient Correlation;
	radius.Fill(1);
		
	sobelX.SetDirection( 0 );
	sobelX.CreateToRadius( radius );

	sobelY.SetDirection(1);
	sobelY.CreateToRadius( radius );

	filterAX = NeighborhoodOperatorType::New();
	filterAY = NeighborhoodOperatorType::New();

	filterBX = NeighborhoodOperatorType::New();
	filterBY = NeighborhoodOperatorType::New();

	NCCmetricX = NCCType::New();
	NCCmetricY = NCCType::New();

	// Pattern intensity
	patternInt = PatternIntensityType::New();

	// gradient difference
	gradDiff = GradientDifferenceType::New();

	// Local correlation
//	localCorr = LocalCorrelationType::New();
}


Metrics::~Metrics(void)
{
}


/* Feature-base registration */
float Metrics::compute_Dice(ImageType::Pointer image1, 
				  ImageType::Pointer image2,
				  float threshold)
{
	//float dice = 0.0;
	
		th1->SetInput( image1 );
		th1->SetLowerThreshold( threshold ) ;
		th1->SetOutsideValue( 1 );
		th1->SetInsideValue( 0 );
	
		th2->SetInput( image2 );
		th2->SetLowerThreshold( threshold ) ;
		th2->SetOutsideValue( 1 );
		th2->SetInsideValue( 0 );

		overlap->SetSourceImage( th1->GetOutput() );
		overlap->SetTargetImage( th2->GetOutput() );
		overlap->Update();
	
	// dice = overlap->GetDiceCoefficient();
	val = overlap->GetDiceCoefficient();
	std::cout << "Dice Coefficient:  " << val << std::endl;
	// std::cout << "This function is empty, so far...!" << std::endl;
	// return dice;
	return val;
}

/* Intensity-based registration*/
float Metrics::compute_msq(ImageType::Pointer image1,
								  ImageType::Pointer image2)
{
	TransformBType::ParametersType params(transform->GetNumberOfParameters());
		params.Fill(0);
	//float val = 0.0;
		interpolator->SetInputImage( image1 );

		SQmetric->SetFixedImage( image1 );
		SQmetric->SetMovingImage( image2 );
		SQmetric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
		SQmetric->SetTransform( transform );
		SQmetric->SetInterpolator( interpolator );

		SQmetric->Initialize();
	
		std::cout << "Sum of Squares Differences: " << SQmetric->GetValue(params) << std::endl;
			
	val = SQmetric->GetValue(params);
	return val;
}

float Metrics::compute_mi( ImageType::Pointer image1,
								 ImageType::Pointer image2)
{
	TransformBType::ParametersType params(transform->GetNumberOfParameters());
		params.Fill(0);
	// float val = 0.0;
		interpolator->SetInputImage( image1 );
	
	// viola Wells Mutual information
	//	MImetric->SetFixedImage( image1 );
	//	MImetric->SetMovingImage( image2 );
	//	MImetric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
	//	MImetric->SetTransform( transform );
	//	MImetric->SetInterpolator( interpolator );
	//
	//	MImetric->Initialize();
	//
	// val = MImetric->GetValue(params);
	//	std::cout << "Mutual Information : " << val << std::endl;
	
	// Mattes Mutual information
		Mattesmetric->SetFixedImage( image1 );
		Mattesmetric->SetMovingImage( image2 );
		Mattesmetric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
		Mattesmetric->SetTransform( transform );
		Mattesmetric->SetInterpolator( interpolator );

		Mattesmetric->Initialize();
	
		std::cout << "Mattes Mutual Information : " << Mattesmetric->GetValue(params) << std::endl;
	
	val = Mattesmetric->GetValue(params);
	return val;
}
	
float Metrics::compute_ncc( ImageType::Pointer image1,
								  ImageType::Pointer image2)
{
	TransformBType::ParametersType params(transform->GetNumberOfParameters());
		params.Fill(0);
		// float val = 0.0;
		interpolator->SetInputImage( image1 );
	
		NCCmetric->SetFixedImage( image1 );
		NCCmetric->SetMovingImage( image2 );

		NCCmetric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
		NCCmetric->SetTransform( transform );
		
		NCCmetric->SetInterpolator( interpolator );

		NCCmetric->Initialize();
	
	std::cout <<  "Normalized cross-Correlation : " << NCCmetric->GetValue(params) << std::endl;
	val = NCCmetric->GetValue(params);

	return val;
}

/* Intensity+position registration*/

float Metrics::compute_gc( ImageType::Pointer image1,
								ImageType::Pointer image2)
{
	TransformBType::ParametersType params( transform->GetNumberOfParameters() );
		params.Fill( 0 );
	
	// X Gradient !
	filterAX->SetOperator( sobelX );
	filterAX->SetInput( image1 );
	filterAX->Update();

	filterBX->SetOperator( sobelX );
	filterBX->SetInput( image2 );
	filterBX->Update();

		NCCmetricX->SetFixedImage( filterAX->GetOutput() );
		NCCmetricX->SetMovingImage( filterBX->GetOutput() );
		NCCmetricX->SetFixedImageRegion( filterAX->GetOutput()->GetLargestPossibleRegion() );
		NCCmetricX->SetTransform( transform );
		NCCmetricX->SetInterpolator( interpolator );
		NCCmetricX->Initialize();

	double nccX = NCCmetricX->GetValue( params );
	std::cout << "Correlacion de Gradiante en X: " << nccX << std::endl;

	// Y Gradient !
	filterAY->SetOperator( sobelY );
	filterAY->SetInput( image1 );
	filterAY->Update();

	filterBY->SetOperator( sobelY );
	filterBY->SetInput( image2 );
	filterBY->Update();
			
		NCCmetricY->SetFixedImage( filterAY->GetOutput() );
		NCCmetricY->SetMovingImage( filterBY->GetOutput() );
		NCCmetricY->SetFixedImageRegion( filterAY->GetOutput()->GetLargestPossibleRegion() );
		NCCmetricY->SetTransform( transform );
		NCCmetricY->SetInterpolator( interpolator );
		NCCmetricY->Initialize();

	double nccY = NCCmetricY->GetValue( params );
	std::cout << "Correlacion de Gradiante en Y: " << nccY << std::endl;

	// GC 

	val = (nccX+nccY)/2;

	std::cout << "Gradient Correlation: " << val << std::endl;

	return val;
}

float Metrics::compute_pi( ImageType::Pointer image1, 
								 ImageType::Pointer image2 )
{
	TransformBType::ParametersType params(transform->GetNumberOfParameters());
		params.Fill(0);
		// float val = 0.0;
		interpolator->SetInputImage( image1 );

		patternInt->SetFixedImage( image1 );
		patternInt->SetMovingImage( image2 );
		patternInt->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
		patternInt->SetTransform( transform );
		patternInt->SetInterpolator( interpolator );
		
		patternInt->SetLambda( 3 ); // Lamba regula el radio de captura de la metrica
		patternInt->SetDelta( 10 ); // Delta es usado como el diferencial en el calculo de las defivada por medio de diferencias finitas

		patternInt->Initialize();

	std::cout << "Pattern intensity : " << patternInt->GetValue(params) << std::endl;
			
	val = patternInt->GetValue(params);
	return val;
}

float Metrics::compute_gd( ImageType::Pointer image1,
								ImageType::Pointer image2 )
{	
	TransformBType::ParametersType params(transform->GetNumberOfParameters());
		params.Fill(0);
		// float val = 0.0;
		interpolator->SetInputImage( image1 );

		gradDiff->SetFixedImage( image1 );
		gradDiff->SetMovingImage( image2 );
		gradDiff->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
		gradDiff->SetTransform( transform );
		gradDiff->SetInterpolator( interpolator );
		
		gradDiff->SetDerivativeDelta( 0.5 );

		gradDiff->Initialize();

	std::cout << "Gradient Difference : " << gradDiff->GetValue(params) << std::endl;
			
	val = gradDiff->GetValue(params);

	return val;
}

//float Metrics::compute_localCorrelation( ImageType::Pointer image1,
//												ImageType::Pointer image2)
//{
//	TransformBType::ParametersType params(transform->GetNumberOfParameters());
//		params.Fill(0);
//		// float val = 0.0;
//		interpolator->SetInputImage( image1 );
//		
//		localCorr->SetFixedImage( image1 );
//		localCorr->SetMovingImage( image2 );
////		localCorr->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
//		
//		localCorr->SetTransform( transform );
////		localCorr->SetInterpolator( interpolator );
//		
//			radius.Fill(5);
//		localCorr->SetRadius(radius);
//		// localCorr->SetDerivativeDelta( 0.5 );
//
//		localCorr->Initialize();
//
//	std::cout << "Local Correlation : " << localCorr->GetValue() << std::endl;
//			
//	val = localCorr->GetValue();
//
//	return val;
//}


#endif