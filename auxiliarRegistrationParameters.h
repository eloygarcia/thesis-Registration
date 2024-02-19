#ifndef __auxiliarRegistrationParameters_h
#define __auxiliarRegistrationParameters_h

#include <iostream>
#include <vector>

struct registrationParameters
{
	registrationParameters(){
	//	mri_origen = 0;
	//	mri_size = 0;
	//	mri_spacing = 0;
	//
	//	mri_imagePointer = 0;

	};

	~registrationParameters(){
	//	delete [] mri_origen;
	//	delete [] mri_size;
	//	delete [] mri_spacing;
	//
	//	delete [] mri_imagePointer;

	};

	// =============== MRI image ! ============
	float * mri_origen;
	int * mri_size;
	float * mri_spacing;

	float * mri_imagePointer;

	// =============== target image ! =============
	float * mamo_origen;
	int * mamo_size;
	float * mamo_spacing;

	unsigned short * mamo_imagePointer;

	// ================= Simulada ! =================

	unsigned short * simulada_imagePointer;

	// =============== Tetrahedral mesh ! ===========
	int * elements;
	int numberOfElements;
	
	int numberOfPoints;
	float * initial_points;
	float * final_points;

	// ================ Regular Grid ! ===============
	float * grid_origen;
	int * grid_size;
	float * grid_spacing;

	int * flags;
	int * cumsum;
	int * correspondingElements;

	// ============= Source ! ================
	float * source;
};


#endif