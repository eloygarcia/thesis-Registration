#ifndef __metrics_h
#define __metrics_h

#pragma once

#include "itkImage.h"
typedef itk::Image<unsigned short, 2 > ImageType;

#include "itkSize.h"

// FEATURE-BASED
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkImageFileWriter.h"

typedef itk::Image<unsigned short,2> LabelBType;
typedef itk::BinaryThresholdImageFilter<ImageType,LabelBType> ThresholdType;
typedef itk::BinaryImageToLabelMapFilter<LabelBType> ImageToLabelType;
typedef itk::LabelOverlapMeasuresImageFilter<LabelBType> OverlapMeasuresType;
typedef itk::ImageFileWriter<ImageType> Writer2DLType;

// INTENSITY-BASED 
#include "itkLinearInterpolateImageFunction.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::TranslationTransform<double,2> TransformBType;

//MEAN SQUARE METRIC
#include "itkMeanSquaresImageToImageMetric.h"
typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MSQType;
// MUTUAL INFORMATION 
#include "itkMutualInformationImageToImageMetric.h"
typedef itk::MutualInformationImageToImageMetric<ImageType, ImageType> MIType;

#include "itkMattesMutualInformationImageToImageMetric.h"
typedef itk::MattesMutualInformationImageToImageMetric<ImageType,ImageType >    MattesMIType;

// NORMALIZED CROSS-CORRELATION
#include "itkNormalizedCorrelationImageToImageMetric.h"
typedef itk::NormalizedCorrelationImageToImageMetric<ImageType,ImageType> NCCType;

// INTENSITY + POSITION 
#include "itkSobelOperator.h"
#include "itkNeighborhoodOperatorImageFilter.h"
typedef itk::SobelOperator<unsigned short,2> SobelOperatorType;
typedef itk::NeighborhoodOperatorImageFilter<ImageType, ImageType> NeighborhoodOperatorType;

// Pattern intensity
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
typedef itk::MeanReciprocalSquareDifferenceImageToImageMetric< ImageType, ImageType > PatternIntensityType;

// Gradient Difference
#include "itkGradientDifferenceImageToImageMetric.h"
typedef itk::GradientDifferenceImageToImageMetric< ImageType, ImageType> GradientDifferenceType;

// Local Correlation
// #include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
// typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4< ImageType, ImageType> LocalCorrelationType;

class Metrics
{
public:
	Metrics(void);
	~Metrics(void);

private: 
	float val;

	TransformBType::Pointer transform;
	InterpolatorType::Pointer interpolator;

	// Dice overlap index.
	ThresholdType::Pointer th1;
	ThresholdType::Pointer th2;
	OverlapMeasuresType::Pointer overlap;

	// Sum of squares Difference.
	MSQType::Pointer SQmetric;
	// Mutual Information.
	MIType::Pointer MImetric;
	MattesMIType::Pointer Mattesmetric;

	// Normalized Cross-Correlation.
	NCCType::Pointer NCCmetric;

	// Gradient Correlation.
	itk::Size<2> radius;
	SobelOperatorType sobelX;
	SobelOperatorType sobelY;

	NeighborhoodOperatorType::Pointer filterAX;
	NeighborhoodOperatorType::Pointer filterAY;

	NeighborhoodOperatorType::Pointer filterBX;
	NeighborhoodOperatorType::Pointer filterBY;

	NCCType::Pointer NCCmetricX;
	NCCType::Pointer NCCmetricY;

	// Pattern Intensity.
	PatternIntensityType::Pointer patternInt;

	// Gradient difference.
	GradientDifferenceType::Pointer gradDiff;

	// Local Correlation.
	// LocalCorrelationType::Pointer localCorr;



public:
	float compute_Dice( ImageType::Pointer image1, ImageType::Pointer image2, float threshold);
	
	float compute_msq( ImageType::Pointer image1, ImageType::Pointer image2);
	float compute_mi(ImageType::Pointer image1, ImageType::Pointer image2);
	float compute_ncc(ImageType::Pointer image1, ImageType::Pointer image2);
	
	float compute_gc(ImageType::Pointer image1, ImageType::Pointer image2);
	float compute_pi(ImageType::Pointer image1, ImageType::Pointer image2);
	float compute_gd(ImageType::Pointer image1, ImageType::Pointer image2);

//	float compute_localCorrelation( ImageType::Pointer image1, ImageType::Pointer image2 );

};

#include "Metrics.cpp"
#endif