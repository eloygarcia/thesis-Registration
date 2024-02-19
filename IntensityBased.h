#pragma once

#ifndef __intensityBased_h
#define __intensityBased_h

//_USE_MATH_DEFINES// for C++
#include <cmath>

#include "auxiliarDefinitions.h"
#include "auxiliarRegistrationParameters.h"
// #include "auxiliarFunctions_Transformations.h"
#include "MechanicalProperties.h"

// #include "NiftySimEjecutable.h"
#ifdef _GPU_
	#include "CudaProjection.h"
#endif

#include "RegularGrid.h"
#include "Transformations.h"
#include "Metrics.h"

//
#include "itkImage.h"
typedef itk::Image<unsigned short, 2> Image2D;
typedef itk::Image<float, 3> Image3D;

// VTK
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"

#include "itkImportImageFilter.h"
typedef itk::ImportImageFilter<unsigned short, 2> ImportFilterType;

class IntensityBased
{
public:
	IntensityBased();
	~IntensityBased();
//
private:
	// float* im;
	registrationParameters * registration;

	float * m_source;
	Image2D::Pointer m_targetImage;
	Image3D::Pointer m_sourceImage;

	Image2D::Pointer m_outputImage;

	// ======= De: Initialize() =====
	std::vector<double> i_points; // initial points 
	std::vector<double> f_points; // final points

	std::vector<double> translation;

	float RotX ;
	float RotY ;
	float RotZ ;

	float angle_rad;

	std::vector<float> com ;
	std::vector<float> cen ;

	bool m_isinitialized;

	// ======= De: Update()  ===========
	std::vector<int> elem ;
	std::vector<double> T_points; // Temporal points -- puntos intermedios!

	modelParameters m_parameters; // parametros bases 
	modelParameters m_temp_parameters; // parametros temporales
	modelParameters ff_parameters; // paremetros finales

//	NewProjection * projection;

	Image2D::SizeType s2d;
	Image2D::PointType o2d;
	Image2D::SpacingType sp2d;
	
	int * size;
	float * origen;
	float * spacing;

	Image2D::Pointer impoint;
	std::string temp_pos;
	std::string outname;

	// ======= De: Set Metric ========
	int m_switch;
	// ======= De: registration metrics ===========
	float val;
	float value;

	// ======= De:  ============
//	Image3D::Pointer m_compressedImage;

	std::string m_outputDir;
	float * m_transformation;

	// ================ Compute Regular Grid =====================
	// int* correspondingElement;
//	float * final_image;

	float * m_origen2d;
	float * m_spacing2d;
	int * m_size2d;
	unsigned short * m_imageMammoPointer;

	float * m_origen3d;
	float * m_spacing3d;
	int * m_size3d;
	float * m_imageMRIPointer;
//	std::string m_side;

	// function getCompressedImage
	float * bounding_box;
//	float * image_input;

	int n;
	int m;

	std::vector<int> mm_size;
	std::vector<double> mm_origen; 
	std::vector<double> mm_spacing;

	Image2D::IndexType start;
	Image2D::SizeType size_im;
	Image2D::PointType origen_im;
	Image2D::SpacingType spacing_im;
	Image2D::DirectionType direction2d;
	ImportFilterType::Pointer importfilter;
	Metrics met;

	float nuPoisson ;
	float * new_transformation ;

#ifdef _GPU_
	// CudaProjection * CDProjection;
	// CudaProjection * Cp; // = new CudaProjection();
#endif

	Transformations * transformation; // = new Transformations();

	//function get_boundingbox
	double x_max, x_min, y_max, y_min, z_max, z_min;

	// WRITER
	Writer2DLType::Pointer writer2del; // = Writer2DLType::New();

//
private:
	void Initialize();

	void boundingBox(float* points, int numberOfPoints, float* bounding_box);

public:
	void SetTargetImage(Image2D::Pointer image);

	void SetSourceImage(Image3D::Pointer image);

	void SetModelParameters(modelParameters myParameters) { 
		m_parameters = myParameters; 
		m_temp_parameters = myParameters;
	};
	void SetMetric( int metric );

	void SetOutputDirectory( std::string output_dir) { m_outputDir = output_dir; };
	void SetTransformation(float * transformation){
			new_transformation[0] = transformation[0];
			new_transformation[1] = transformation[1];
			new_transformation[2] = transformation[2];
			new_transformation[3] = transformation[3];
			new_transformation[4] = transformation[4];
			new_transformation[5] = shearModulus( transformation[5], transformation[7] );
			new_transformation[6] = bulkModulus( transformation[5], transformation[7] );
			new_transformation[7] = shearModulus( transformation[5]*transformation[6], transformation[7]);
			new_transformation[8] = bulkModulus( transformation[5]*transformation[6], transformation[7] );
			new_transformation[9] = transformation[8];
		
		m_transformation = new_transformation;
	};

	float Update();

	Image2D::Pointer GetOutputImage(){ return this->m_outputImage; };

	friend class HillClimbingOpt;
};


#include "IntensityBased.cpp"

#endif