#pragma once

//#include "CudaProjection.cu"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>

#include <vector>
#include <string.h>

#include "auxiliarRegistrationParameters.h"

class CudaProjection
{
public:
	CudaProjection();
	~CudaProjection();

// VARIABLES !
private:
// host 

	registrationParameters * m_parameters;

	int * m_mri_size;
	float * m_mri_spacing;
	float * m_mri_origen;
	float * m_mri_imagepointer;

	// === Mesh and Grid ! ===
	float * m_i_points;
	float * m_f_points;
	int * m_elements;

	float * m_grid_origen;
	float * m_grid_spacing;
	int * m_grid_size;

	int * m_flags;
	int * m_cumsum;
	int * m_correspondingElements;
	// =======================

	int * m_simulada_size;
	float * m_simulada_spacing;
	float * m_simulada_origen;
	unsigned short * m_simulada_imagepointer;

	float * m_source;

// device
	int* dev_mri_size;
	float* dev_mri_spacing;
	float* dev_mri_origen;
int numberOfPixels_MRI; // esta es host !
	float* dev_mri_imagepointer;
	
	// === Mesh and Grid ! ===
int numberOfPoints;
	float * dev_i_points;
	float * dev_f_points;
int numberOfElements;
	int * dev_elements;

	float * dev_grid_origen;
	float * dev_grid_spacing;
	int * dev_grid_size;
int numberOfVoxelsGrid;

	int * dev_flags;
	int * dev_cumsum;
int maximumCorrespondingElements;
	int * dev_correspondingElements;
	// =======================

	int* dev_simulada_size;
	float* dev_simulada_spacing;
	float* dev_simulada_origen;
int numberOfPixels_Simulada; // esta también es host !
	//float* dev_simulada_imagepointer;
	unsigned short * dev_simulada_imagepointer;

	float* dev_source;

	cudaError_t cudaStatus;

	int bl; // Numero de bloques a lanzar, en conjuntos de 512 threads
// METODOS !

private:
	void Initialize();

public:
	/*
	//int action();
	//cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

	// metodos de construccion de la imagen:
	void Set3DImageSize(int *imageSize ){ m_mri_size = imageSize;
	//printf("Image Size : [ %d,%d,%d]\n",imageSize[0], imageSize[1], imageSize[2]);
	printf("3D Size : [ %d,%d,%d]\n",m_mri_size[0], m_mri_size[1], m_mri_size[2]);
		numberOfPixels_MRI = m_mri_size[0] * m_mri_size[1] * m_mri_size[2];
		//dev_3d_imagepointer = new float[ numberOfPixels_3D ];
	};

	void Set3DImageSpacing( float *imageSpacing ){ m_mri_spacing = imageSpacing;
	printf("3D Spacing : [ %f,%f,%f]\n",m_mri_spacing[0], m_mri_spacing[1], m_mri_spacing[2]);
	};

	void Set3DImageOrigin( float *imageOrigin ){ m_mri_origen = imageOrigin; 
	printf("3D Origen : [ %f,%f,%f]\n",m_mri_origen[0], m_mri_origen[1], m_mri_origen[2]);
	};
	void Set3DImagePointer( float *imagePointer ){ m_mri_imagepointer = imagePointer; };

	void Set2DImageSize( int *imageSize ){ m_simulada_size = imageSize;
	printf("2D Size : [ %d,%d,%d]\n",m_simulada_size[0], m_simulada_size[1], m_simulada_size[2]);
		numberOfPixels_Simulada = m_simulada_size[0] * m_simulada_size[1] * m_simulada_size[2];
		//dev_2d_imagepointer = new float[ numberOfPixels_2D ];
	};

	void Set2DImageSpacing( float *imageSpacing ){ m_simulada_spacing = imageSpacing; 
	printf("2D Spacing : [ %f,%f,%f]\n",m_simulada_spacing[0], m_simulada_spacing[1], m_simulada_spacing[2]);
	};

	void Set2DImageOrigin( float *imageOrigin ){ m_simulada_origen = imageOrigin; 
	printf("2D Origen : [ %f,%f,%f]\n",m_simulada_origen[0], m_simulada_origen[1], m_simulada_origen[2]);
	};

	void Set2DImagePointer( float *imagePointer ){m_simulada_imagepointer = imagePointer;};

	void SetSourceProjection( float * source ){ m_source = source; 
	printf("Source : [ %f,%f,%f]\n",m_source[0], m_source[1], m_source[2]);
	};
	*/

	void SetParameters( registrationParameters * parameters) { m_parameters = parameters; } ;

	// Métodos Get
	const int * Get3DImageSize(){return this->m_mri_size;};
	const float * Get3DImageSpacing(){return this->m_mri_spacing;};
	const float * Get3DImageOrigin(){return this->m_mri_origen;};
	const float * Get3DImagePointer(){return this->m_mri_imagepointer;};

	int * Get2DImageSize(){return this->m_simulada_size;};
	float * Get2DImageSpacing(){return this->m_simulada_spacing;};
	float * Get2DImageOrigin(){return this->m_simulada_origen;};
	unsigned short * Get2DImagePointer(){return this->m_simulada_imagepointer;};

	float * GetSourceProjection(){return this->m_source;};

	// Proyección, método update()
	void Update();

};
