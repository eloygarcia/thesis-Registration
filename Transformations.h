#ifndef __Transformations_h
#define __Transformations_h

#pragma once

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <vector>

#include <iostream>
#include <math.h>
#include <cmath>

#include "auxiliarDefinitions.h"

#include "NiftySimEjecutable.h"

class Transformations
{
public:
	Transformations(void);
	~Transformations(void);

private:
	float * m_transformation;
	// m_transformation = [ transl_X, transl_Y, RotX, RotY, Thickness];

	modelParameters m_parameters;
	modelParameters temp_parameters;
	NiftySimEjecutable * m_niftysimulation;

	int m_angle;
	std::string m_side;
	std::vector<float> m_cen;
	std::string m_outputDir;
	float RotX;
	float RotY;
	float RotZ;

	std::vector<double> m_initial_points;
	std::vector<double> m_initial_boundingBox;
	std::vector<float> m_initial_centerOfMass;

	std::vector<double> m_med_points;
	std::vector<double> m_med_boundingBox;
	std::vector<float> m_med_centerOfMass;

	std::vector<double> m_final_points;
	std::vector<double> m_final_boundingBox;
	std::vector<float> m_final_centerOfMass;

	float grad2rad( float degree );
	void get_boundingBox(std::vector<double> points, std::vector<double> &boundingBox);
	void get_centerofmass(std::vector<double> points, std::vector<float> &centerOfMass);

	void set_translation(std::vector<double> i_points,  std::vector<double> translation, std::vector<double> &f_points);
	void set_rotation(std::vector<double> i_points, std::vector<double> &f_points, std::vector<float> center, float RotX, float RotY, float RotZ);
	void set_transform(std::vector<double> i_points, std::vector<double> &f_points, std::vector<double> translation, std::vector<float> center, float RotX, float RotY, float RotZ);

	void Initialize();
	void updateParameters(std::vector<double> T_points, std::vector<double> T_boundingBox);

public:
	void SetModelParameters( modelParameters parameters );
	void SetOutputDir( std::string outputDir ){ m_outputDir = outputDir; };

	void SetTransformation( float* transformation ){ m_transformation = transformation;};
	// void setPoints(std::vector<double> points){ m_initial_points = points; };

	std::vector<double> getInitialPoints(){ return m_initial_points; };
	std::vector<double> getInitialBoundingBox(){ return m_initial_boundingBox;};
	std::vector<float> getInitialCenterOfMass(){ return m_initial_centerOfMass;};

	std::vector<double> getFinalPoints(){ return m_final_points; };
	std::vector<double> getFinialBoundingBox(){ return m_final_boundingBox;};
	std::vector<float> getfinalCenterOfMass(){ return m_final_centerOfMass;};

	int Update();

};


#include "Transformations.cpp"
#endif