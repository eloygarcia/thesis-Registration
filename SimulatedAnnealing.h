#ifndef __SimulatedAnnealing_h
#define __SimulatedAnnealing_h

#include <iostream>

#include "MechanicalProperties.h"

#pragma once

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

#include <ctime>

#include "IntensityBased.h"

#include "itkImage.h"
typedef itk::Image<unsigned short,2> ImageSHC2d;
#include "itkImageFileWriter.h"
typedef itk::ImageFileWriter<ImageSHC2d> WriterSHC2d;

class SimulatedAnnealing
{
public:
	SimulatedAnnealing();
	~SimulatedAnnealing();
//
public:

private:
	std::string m_position;
	int n_iter;
	int m_numberofparameters;
	
	float m_epsilon;
	float* m_currentPoint;

	bool m_maximize;
	float bestScore;

	bool cont;
	int cont_member;
	
	IntensityBased * m_intensityBased;
	//typedef IntensityBased IBType;

	bool is_up;
	float m_aceptance;
	float m_Temperature;

	// Generador de aleatorios en una gaussiana
	float * m_searchSpace;
	std::default_random_engine generator;
	float m_sd;
	
	// float* temp_old_pos;
	float temp_new_pos;

	ofstream file;
	std::string m_outputdir;
	bool has_output;
	
	stringstream ss;
	std::string temp_string;

	WriterSHC2d::Pointer writer;

//
public:
	void setNumberOfParameters( int numberofparameters ){ m_numberofparameters = numberofparameters; };
	void setInitialPoint( float * iniPoint ) { m_currentPoint = iniPoint; };
	void defineSearchSpace( float * searchSpace ){ m_searchSpace = searchSpace; };
	
	void setTemperature( float T ){ m_Temperature = T; };

	void setMaximizeOn(){ m_maximize = true; };
	void setMinimizeOn(){ m_maximize = false; };

	void setEpsilon( float epsilon ) { m_epsilon = epsilon;  };
	void setMaxIter(int maxIter){ n_iter = maxIter; };

	void setPosition(std::string position ){ m_position = position; };

	void setOptimFunction( IntensityBased*  ib){ m_intensityBased = ib;};

	void setOutputDir( std::string outputdir){
		m_outputdir = outputdir;
		has_output= true;
	};

	int Update();

};


#include "SimulatedAnnealing.cpp"

#endif