#ifndef __Transformations_cpp
#define __Transformations_cpp

#include "Transformations.h"

Transformations::Transformations(void)
{
	m_angle = 0;

	RotX = 0.;
	RotY = 0.; 
	RotZ = 0.;

	m_niftysimulation = new NiftySimEjecutable();
}

Transformations::~Transformations(void)
{
	//delete [] m_transformation;
	m_cen.~vector();

	m_initial_points.~vector();
	m_initial_boundingBox.~vector();
	m_initial_centerOfMass.~vector();

	m_med_points.~vector();
	m_med_boundingBox.~vector();
	m_med_centerOfMass.~vector();

	m_final_points.~vector();
	m_final_boundingBox.~vector();
	m_final_centerOfMass.~vector();

	delete m_niftysimulation;
}

float Transformations::grad2rad( float degree)
{
	float aa = (degree*((float)M_PI)) /180;
	return aa;
}

// ========== Bounding Box & Center Of Mass ==================
void Transformations::get_boundingBox(std::vector<double> points, std::vector<double> &boundingBox)
{
	double x_max = -9999; double x_min = 9999;
	double y_max = -9999; double y_min = 9999;
	double z_max = -9999; double z_min = 9999;

	for(unsigned int i=0; i<(points.size()/3); i++)
	{
		if(points[3*i] < x_min) x_min= points[3*i];
		if(points[3*i] > x_max) x_max= points[3*i];

		if(points[3*i+1] < y_min) y_min= points[3*i+1];
		if(points[3*i+1] > y_max) y_max= points[3*i+1];

		if(points[3*i+2] < z_min) z_min= points[3*i+2];
		if(points[3*i+2] > z_max) z_max= points[3*i+2];
	}

	boundingBox = std::vector<double>();
		boundingBox.push_back(x_min); 
		boundingBox.push_back(x_max);
		boundingBox.push_back(y_min);
		boundingBox.push_back(y_max);
		boundingBox.push_back(z_min);
		boundingBox.push_back(z_max);

//	std::cout << "BoundingBox: ["<< boundingBox[0]<< ", "<< boundingBox[1] << ", "<< boundingBox[2] <<
//		 ", "<< boundingBox[3] << ", "<< boundingBox[4] << ", "<< boundingBox[5] << std::endl;

}

void Transformations::get_centerofmass(std::vector<double> points, std::vector<float> &centerOfMass)
{
	float sum_x = 0.0;
	float sum_y = 0.0;
	float sum_z = 0.0;
	int n_points = points.size()/3;

	for(int i=0; i<(n_points); i++)
	{
		sum_x += (float) points[3*i];
		sum_y += (float) points[3*i+1];
		sum_z += (float) points[3*i+2];
	}

	centerOfMass = std::vector<float>();
		centerOfMass.push_back( sum_x/n_points );
		centerOfMass.push_back( sum_y/n_points );
		centerOfMass.push_back( sum_z/n_points );

//	std::cout << std::endl;
//	std::cout << "Center of mass: " << centerOfMass[0] << ", " << 
//		centerOfMass[1] << ", " << centerOfMass[2] << "]" << std::endl;
}

// ============  Translation & Rotation =======================
void Transformations::set_translation(std::vector<double> i_points,  std::vector<double> translation, std::vector<double> &f_points)
{
		for(unsigned int i=0; i<(i_points.size()/3); i++)
		{
			f_points[3*i]   =  i_points[3*i]   +  translation[0];
			f_points[3*i+1] =  i_points[3*i+1] +  translation[1];
			f_points[3*i+2] =  i_points[3*i+2] +  translation[2];
		}
}

void Transformations::set_rotation(std::vector<double> i_points, std::vector<double> &f_points, std::vector<float> center, float RotX, float RotY, float RotZ)
{
	double * temp_point = new double[3]; // = {0.0, 0.0, 0.0};
	double * temp_rot = new double[3]; // = {0.0, 0.0, 0.0};
	std::vector<double> temp_vector;
	for(unsigned int i=0; i<(i_points.size() /3); i++)
	{
		temp_point[0] = i_points[3*i];
		temp_point[1] = i_points[3*i+1];
		temp_point[2] = i_points[3*i+2];

		// ================== Rotation ====================
		// --------- Tranlation: center of mass to 0,0,0 ----------
		temp_point[0] -= center[0];
		temp_point[1] -= center[1];
		temp_point[2] -= center[2];
		
		// ------------------- Rotation ------------------
			if(RotX != 0)
			{
				// RotX -- radians
				temp_rot[0] = temp_point[0];
				temp_rot[1] = temp_point[1] * cos(RotX) - temp_point[2] * sin(RotX);
				temp_rot[2] = temp_point[1] * sin(RotX) + temp_point[2] * cos(RotX);
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
			if(RotY != 0)
			{
				// RotY -- radians
				temp_rot[0] = temp_point[0] * cos(RotY) + temp_point[2] * sin(RotY);
				temp_rot[1] = temp_point[1];
				temp_rot[2] = - temp_point[0] * sin(RotY) + temp_point[2] * cos(RotY);
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
			if(RotZ != 0)
			{
				// RotZ -- radians
				temp_rot[0] = temp_point[0] * cos(RotZ) - temp_point[1] * sin(RotZ);
				temp_rot[1] = temp_point[0] * sin(RotZ) + temp_point[1] * cos(RotZ);
				temp_rot[2] = temp_point[2];
				// Actualizar;
				temp_point[0] = temp_rot[0];
				temp_point[1] = temp_rot[1];
				temp_point[2] = temp_rot[2];
			}
		// ------------ Back ---------------
		temp_point[0] += center[0];
		temp_point[1] += center[1];
		temp_point[2] += center[2];


		// ========== Compose vector ============
		temp_vector.push_back( temp_point[0] );
		temp_vector.push_back( temp_point[1] );
		temp_vector.push_back( temp_point[2] );
	}

	f_points.clear();
	f_points = temp_vector;

	// ==== delete ===
	delete [] temp_point;
	delete [] temp_rot; 
	temp_vector.~vector();

}

void Transformations::set_transform(std::vector<double> i_points, std::vector<double> &f_points, std::vector<double> translation, std::vector<float> center, float RotX, float RotY, float RotZ)
{
	// f_points = i_points;
/*
	ROTACION!
*/
	if( RotX!=0 || RotY!=0 ||  RotZ!=0) {
		double * temp_point = new double[3]; // = {0.0, 0.0, 0.0};
		double * temp_rot = new double[3]; // = {0.0, 0.0, 0.0};
		std::vector<double> temp_vector;

		for(unsigned int i=0; i<(f_points.size() /3); i++)
		{
			temp_point[0] = f_points[3*i];
			temp_point[1] = f_points[3*i+1];
			temp_point[2] = f_points[3*i+2];
	
			// ================== Rotation ====================
			// --------- Translation: center of mass to 0,0,0 ----------
			temp_point[0] = temp_point[0] - center[0];
			temp_point[1] = temp_point[1] - center[1];
			temp_point[2] = temp_point[2] - center[2];
			
			// ------------------- Rotation ------------------
				if(RotX != 0)
				{
					// RotX -- radians
					temp_rot[0] = temp_point[0];
					temp_rot[1] = temp_point[1] * cos(RotX) - temp_point[2] * sin(RotX);
					temp_rot[2] = temp_point[1] * sin(RotX) + temp_point[2] * cos(RotX);
					// Actualizar;
					temp_point[0] = temp_rot[0];
					temp_point[1] = temp_rot[1];
					temp_point[2] = temp_rot[2];
				}
				if(RotY != 0)
				{
					// RotY -- radians
					temp_rot[0] = temp_point[0] * cos(RotY) + temp_point[2] * sin(RotY);
					temp_rot[1] = temp_point[1];
					temp_rot[2] = - temp_point[0] * sin(RotY) + temp_point[2] * cos(RotY);
					// Actualizar;
					temp_point[0] = temp_rot[0];
					temp_point[1] = temp_rot[1];
					temp_point[2] = temp_rot[2];
				}
				if(RotZ != 0)
				{
					// RotZ -- radians
					temp_rot[0] = temp_point[0] * cos(RotZ) - temp_point[1] * sin(RotZ);
					temp_rot[1] = temp_point[0] * sin(RotZ) + temp_point[1] * cos(RotZ);
					temp_rot[2] = temp_point[2];
					// Actualizar;
					temp_point[0] = temp_rot[0];
					temp_point[1] = temp_rot[1];
					temp_point[2] = temp_rot[2];
				}
			// ------------ Back ---------------
			temp_point[0] = temp_point[0] + center[0];
			temp_point[1] = temp_point[1] + center[1];
			temp_point[2] = temp_point[2] + center[2];
	
			// ========== Compose vector ============
			temp_vector.push_back( temp_point[0] );
			temp_vector.push_back( temp_point[1] );
			temp_vector.push_back( temp_point[2] );
		}

	f_points.clear();
	f_points = temp_vector;

	// = ==  delete == =
	delete [] temp_point;
	delete [] temp_rot;
	temp_vector.~vector();

	}

/*
	TRANSLACION!!
*/
		for(unsigned int i=0; i<(f_points.size()/3); i++)
		{
			f_points[3*i]= f_points[3*i] + translation[0];
			f_points[3*i+1] = f_points[3*i+1] + translation[1];
			f_points[3*i+2] = f_points[3*i+2] + translation[2];
		}

}

// ================ Initializing transformations =============

void Transformations::SetModelParameters( modelParameters parameters ){ 
	m_parameters = parameters;

	// == initial points
	m_initial_points = m_parameters.nodes;
	get_centerofmass( m_initial_points, m_initial_centerOfMass );
	get_boundingBox( m_initial_points, m_initial_boundingBox);

	// == angle and side
	m_angle = m_parameters.angle;
	m_side = m_parameters.side;
	m_cen = m_parameters.cen;

	// == final points (initialize)
	m_med_points = std::vector<double>();
	for(int i=0; i< m_initial_points.size(); i++) m_med_points.push_back( 0 );
//	m_final_points = std::vector<double>();
//	for(int i=0; i< m_initial_points.size(); i++) m_final_points.push_back( 0 );
};

void Transformations::Initialize()
{
	// == Temporal points (initialize)
	std::vector<double> temp_points;
	for(int i=0; i< m_initial_points.size(); i++) temp_points.push_back( 0.0 );
	std::vector<double> temp_points2;
	for(int i=0; i< m_initial_points.size(); i++) temp_points2.push_back( 0.0 );

	// == Set Rotation
	RotX = M_PI_2 ; // M_PI_2 -- Rotaci'on para ponerlo en orden. ;
	RotY = 0.; // 
	// RotZ = 0.;

	if (m_side=="R ") {
		RotZ = M_PI_2; // Right M_PI_2;
	}else if (m_side=="L "){
		RotZ = -M_PI_2; // Left -M_PI_2;
	}
	set_rotation( m_initial_points, temp_points, m_initial_centerOfMass, RotX, RotY, RotZ);


	/*  CHECKEANDO LAS MLO LEFT !! */
	// == Set Rotation
	
	/**/


	// -- Rotación  Para acercarlo a la posición de mamografía
	if( m_side=="L "){
		RotX = (float)(m_parameters.angle)*M_PI/180;  // Left está bien.
	}else if(m_side=="R "){
		RotX = -(float)(m_parameters.angle)*M_PI/180;  // Right ??
	}
	RotY = 0.;
	if (m_side=="R ") {
		RotZ = -(fabs(m_parameters.pectoralAngle));
	}else if (m_side=="L "){
		RotZ = (fabs(m_parameters.pectoralAngle));
	}
	set_rotation( temp_points, temp_points2, m_initial_centerOfMass, RotX, RotY, RotZ);

	// == Set Translation
	std::vector<double> tra;
		tra.push_back((m_cen[0] - m_initial_centerOfMass[0])); 
		tra.push_back((m_cen[1] - m_initial_centerOfMass[1])); 
		tra.push_back((m_cen[2] - m_initial_centerOfMass[2]));
	set_translation( temp_points2, tra, m_med_points );

	get_boundingBox(m_med_points, m_med_boundingBox);
	get_centerofmass(m_med_points, m_med_centerOfMass);

	
	// == Delete ==
	tra.~vector();
	temp_points.~vector();
	temp_points2.~vector();
}

void Transformations::updateParameters(std::vector<double> T_points, std::vector<double> T_boundingBox)
{
	// == paramters
	temp_parameters.number_of_tissues = m_parameters.number_of_tissues;
	
	// Nodes
	temp_parameters.nodes = T_points;
	// char* nodes_Dof; --> 3 
	temp_parameters.boundingBox = T_boundingBox;
	
	// Elements
	temp_parameters.elements = m_parameters.elements;
	// char* type; --> T4
	temp_parameters.elementTissue = m_parameters.elementTissue; 
	
	// Subsets
	temp_parameters.subsets_on = m_parameters.subsets_on; 

	// System Params;
	// char* timeStep;
	// char* totalTime;
	// char* DampingCoeff;
	// char* Density;
	// char* DoDeformableCollision;
	
	// Set Material;
	// char* material; --> "NH"
	// int numParams; --> 2
	// std::vector<double> material_parameters; 
	temp_parameters.material_parameters[ 0 ] = m_transformation[ 5 ];
	temp_parameters.material_parameters[ 1 ] = m_transformation[ 6 ];
	temp_parameters.material_parameters[ 2 ] = m_transformation[ 7 ];
	temp_parameters.material_parameters[ 3 ] = m_transformation[ 8 ];

	// Constraint;
	temp_parameters.constraint_bC_list = m_parameters.constraint_bC_list;
	// char* constraint_type; --> "Fix"
	// char* constraint_axis; --> "0"

	// Compression;
	temp_parameters.position = m_parameters.position;
	// mri data !

	// mammogram data !
	temp_parameters.cen = m_parameters.cen;
	temp_parameters.angle = m_parameters.angle;
	temp_parameters.mm_origen = m_parameters.mm_origen;
	temp_parameters.mm_spacing = m_parameters.mm_spacing;
	temp_parameters.mm_size = m_parameters.mm_size;

	temp_parameters.side = m_parameters.side;

// ----->> // Update Thickness 
	temp_parameters.thickness = m_transformation[9];  // THICKNESS !!

	temp_parameters.pectoralAngle = m_parameters.pectoralAngle;
}

int Transformations::Update()
{
	Initialize();
	
	// == Temporal points (initialize)
	std::vector<double> temp_points;
	for(int i=0; i< m_initial_points.size(); i++) temp_points.push_back( 0 );
	std::vector<double> T_points;
	for(int i=0; i< m_initial_points.size(); i++) T_points.push_back( 0 );

	// == Set Rotation
	RotX = m_transformation[3];  // ROT_X
	RotY = m_transformation[4];  // ROT_Y
	RotZ = m_transformation[2];  // ROT_Z
	set_rotation(m_med_points, temp_points, m_med_centerOfMass ,RotX,  RotY, RotZ);

	// == Set Translation
	std::vector<double> tra_2;  // TRANSLATION
		tra_2.push_back( m_transformation[0] );
		tra_2.push_back( m_transformation[1] );
		tra_2.push_back( 0 );
	set_translation(temp_points, tra_2, T_points);

	std::vector<double> T_boundingBox;
	get_boundingBox( T_points, T_boundingBox );
	
	// == Update Parameters !!
	updateParameters( T_points, T_boundingBox);

//  =================================== NiftySim ejecution!! ====================================
	m_niftysimulation->SetModelParameters( temp_parameters );
	m_niftysimulation->SetOutputDir( m_outputDir );
	int sys_b = m_niftysimulation->Update();
		if(sys_b!=0) return 1;
	
	m_final_points = m_niftysimulation->GetFinalPoints();

	// == Destroying temporal vectors 
	temp_points.~vector();
	T_points.~vector();

	tra_2.~vector();
	
	return 0;
}

#endif

