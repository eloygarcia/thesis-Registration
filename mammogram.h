#pragma once

#ifndef __mammogram_h
#define __mammogram_h

#define _USE_MATH_DEFINES
#include <math.h>

// includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkFixedArray.h"
#include "itkFlipImageFilter.h"
#include "itkImageRegionIterator.h"

#include "itkOtsuThresholdImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

#include "itkBinaryContourImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkMetaDataDictionary.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkMaskImageFilter.h"

#include "itkImageLinearConstIteratorWithIndex.h"

//typedef

typedef itk::Image<unsigned short,2> ImageMammoType;

typedef itk::ImageFileReader<ImageMammoType> ReaderMammoType;
typedef itk::ImageFileWriter<ImageMammoType> WriterMammoType;

typedef itk::FlipImageFilter< ImageMammoType > FlipImageFilterType;

typedef itk::GDCMImageIO  IOType;
typedef itk::MetaDataDictionary MetaDataType;

typedef itk::OtsuThresholdImageFilter<ImageMammoType, ImageMammoType> OtsuMammoType;
typedef itk::BinaryContourImageFilter<ImageMammoType, ImageMammoType> ContourMammoType;

typedef itk::ShapeLabelObject< unsigned short, 2 > ShapeLabelObjectType;
typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;
typedef itk::LabelImageToShapeLabelMapFilter< ImageMammoType, LabelMapType> I2LType;

typedef itk::ResampleImageFilter<ImageMammoType, ImageMammoType> ResampleFilterType;
typedef itk::IdentityTransform<double,2> TransformType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageMammoType,double> NNType;

typedef itk::MaskImageFilter< ImageMammoType, ImageMammoType, ImageMammoType> MaskFilterType;

typedef itk::ImageLinearConstIteratorWithIndex< ImageMammoType > IteratorType;

class mammogram
{
public:
	mammogram(void);
	~mammogram(void);

private:
	
	// filename
	std::string m_filename;

	// mammogram
	ImageMammoType::Pointer m_image;
	bool has_mammogram;
	
	// image Io
	bool has_dictionary;
	IOType::Pointer m_imageIO;
	MetaDataType metadata;
	
	// segmentation
	bool has_segmentation;
	ImageMammoType::Pointer m_segmentation;

	// density map.
	bool has_densityMap;
	ImageMammoType::Pointer m_densityMap;
	ImageMammoType::Pointer im_dense_tissue;
	bool has_im_dense_tissue;

	// pectoral muscle
	ImageMammoType::Pointer m_pectoralMuscleMask;
	
	void computeDenseTissue( bool verbose );
	void computeAnglePectoralMuscle(ImageMammoType::Pointer image);
	double m_pectoralMuscleAngle;

public:
	void read_mammogram(std::string inputfilname );
	void read_mammogram_with_metadata(std::string inputfilname );
	void read_densityMap(std::string inputfilname);
	void read_pectoralmuscleMask(std::string inputfilename);

	void write_mammogram( std::string outputfilename );
	void write_densitymap( std::string outputfilename);
	
	ImageMammoType::Pointer get_image();
	ImageMammoType::Pointer get_otsu_segmentation();

	ImageMammoType::Pointer get_DensityMap();
	ImageMammoType::Pointer get_dense_tissue(bool verbose);
	ImageMammoType::Pointer get_otsu_densityMap();

	void otsu_segmentation();
	std::vector<float> getCentroid(bool verbose);
	void getContour();

	int getDistanceSourceToDetector(bool verbose);
	float getFocalSpot(bool verbose);
	int getBodyPartThickness(bool verbose);
	int getCompressionForce(bool verbose);
	std::string getViewPossition(bool verbose);
	int getPositionerPrimaryAngle(bool verbose);
	std::string getAnodeTargetMaterial(bool verbose);
	std::string getFilterMaterial(bool verbose);
	int getKVP(bool verbose);
	std::string getImageLaterality(bool verbose);

	double getPectoralAngle(){ return m_pectoralMuscleAngle;}; // En radianes!

};


#include "mammogram.cpp"

#endif