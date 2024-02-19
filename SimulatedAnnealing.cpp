#ifndef __SimulatedAnnealing_cpp
#define __SimulatedAnnealing_cpp

#include "SimulatedAnnealing.h"

SimulatedAnnealing::SimulatedAnnealing()
{
	m_maximize = false;
	is_up = false;

	has_output=false;
	writer = WriterSHC2d::New();
}

SimulatedAnnealing::~SimulatedAnnealing()
{
}

//

int SimulatedAnnealing::Update()
{
	clock_t bg = clock();
	int max_j=50;

if( has_output ) { 
	//std::string output = m_outputdir + "\\" + m_position + "_optimizationValues.csv";
	std::string output = m_outputdir + "\\optimizationValues.csv";
	file.open( output.c_str() );
	file << "translationX,translationY,RotZ,RotX,RotY,E_fat,E_gland,thickness,metric\n";
}

		float eval = 0;// = IntensityBased::Update();
		eval = m_intensityBased->Update();
			if(eval==0.0) return 1;

if( has_output ) {
	temp_string.clear();
	//temp_string.c_str("");
	for(int r = 0; r<m_numberofparameters; r++) {
		ss.str("");
		ss << m_currentPoint[r];
		
		temp_string = temp_string + ss.str() + ",";
	}
	
	ss.str("");
	ss << eval;
	temp_string = temp_string + ss.str() + "\n";
	
	file << temp_string; // .c_str();
}

	std::string initialImage = m_outputdir+"\\"+m_position+"_initialImage.nrrd";
		writer->SetInput( m_intensityBased->GetOutputImage() );
		writer->SetFileName( initialImage );
		try{
			writer->Update();
		}catch( itk::ExceptionObject & excp ) {
			std::cout << " Simulated Annealing, Initial Image Writer exception! " << std::endl;
			std::cout << excp << std::endl;
			std::cout << std::endl;
			return 1;
		}
	


	bestScore = eval;
	float temp_bestScore = eval;
	float * temp_old_position = new float[m_numberofparameters];
		for(int m=0; m<(m_numberofparameters) ; m++) temp_old_position[m] = m_currentPoint[m];

	float * temp_new_position = new float[m_numberofparameters];
		for(int m=0; m<(m_numberofparameters) ; m++) temp_new_position[m] = 0.0;

	cont = true;
	cont_member = 0;
	
	float temp_eval = eval;

	if ( !m_maximize ) std::cout << "ESTAMOS EN MINIMIZE !! " << std::endl;
	int i=0;
	int j=0;

	while( cont && i<n_iter)
	{
		temp_eval = eval;
		cont_member = 0;
		int best=-1;

		is_up = false;
		j=0;
		while (j<max_j ) // j estaba en 20
		{
			// std::cout << " Actual point : [" ;
		 	// for(int r = 0; r<m_numberofparameters-1; r++) std::cout << temp_old_position[r] << ", ";
			// std::cout <<temp_old_position[m_numberofparameters-1] <<"] " << std::endl;
			
			std::cout << std::endl;
			std::cout << "Iteration : : " << i << std::endl;
			//std::cout << "CurrentPoint index : : " << m;
			std::cout << "Candidate = " << j << std::endl; // No evalua más de un candidato!
			std::cout << "Temp Best Score : " << temp_bestScore << std::endl;
			std::cout << std::endl;

			/* Elección aleatoria de la nueva posición */
			// std::default_random_engine generator;
//			std::cout << "Initial Point Value := " << m_currentPoint[m] << std::endl;
				

			for(int m=0; m<(m_numberofparameters) ; m++) {
				m_sd = abs( m_searchSpace[2*m+1] - m_searchSpace[2*m] )*0.025;
				// std::cout << std::endl;
				// std::cout << "Distribución normal en: (" << m_currentPoint[m] <<", " << m_sd <<") " << std::endl;
				std::normal_distribution<float> distribution(temp_old_position[m], m_sd);
				
				temp_new_pos = distribution( generator );

				if( temp_new_pos > m_searchSpace[2*m+1] ) temp_new_pos = m_searchSpace[2*m+1];
				else if( temp_new_pos < m_searchSpace[2*m] ) temp_new_pos = m_searchSpace[2*m];

				temp_new_position[m] = temp_new_pos;
			}
			
			// printing ...
			std::cout << "Evaluating : [" ;
			//for(int r = 0; r<m_numberofparameters-1; r++) 
				std::cout << temp_new_position[0] << ", "; // Translacion X
				std::cout << temp_new_position[1] << ", "; // Translacion Y
				std::cout << (temp_new_position[2] * 180) / M_PI<< ", "; // Rotacion Z
				std::cout << (temp_new_position[3] * 180) / M_PI<< ", "; // Rotacion X
				std::cout << (temp_new_position[4] * 180) / M_PI<< ", "; // Rotacion Y
				for(int r = 5; r<m_numberofparameters-1; r++) std::cout << temp_new_position[r] << ", ";
			std::cout <<temp_new_position[m_numberofparameters-1] <<"] " << std::endl;

			// evaluating ....
				m_intensityBased->SetTransformation( temp_new_position );
			temp_eval = m_intensityBased->Update();
				if(temp_eval==0.0){
					j++;
					break;
					//
					//m_intensityBased->SetTransformation( temp_old_position );
					//temp_eval = m_intensityBased->Update();

					//std::string finalImage = m_outputdir+"\\uncompletedImage.nrrd";
					//writer->SetInput( m_intensityBased->GetOutputImage() );
					//writer->SetFileName( finalImage );
					//try{
					//	writer->Update();
					//}catch( itk::ExceptionObject & excp ) {
					//	std::cout << " Simulated Annealing, Uncompleted imae Image Writer exception! " << std::endl;
					//	std::cout << excp << std::endl;
					//	std::cout << std::endl;
					//	return 1;
					//}
					//return 1;
				}


//			std::cout << std::endl;
			std::cout << "Best Score " << bestScore << " vs temp_bestScore " << temp_bestScore << " vs. temp_eval " << temp_eval << std::endl;
			std::cout << std::endl;
							
			// m_currentPoint[m] = temp_old_pos;
			
			// Es esto lo que hay que cambiar en la optimización
			/* ACEPTANCE */
			m_aceptance = float(float(rand()%100)/100);
				
			std::cout << "Probabilidad de aceptacion : " << m_aceptance << std::endl;
			std::cout << "Probabilidad resultado : " << float(1/(1+exp ((temp_eval-temp_bestScore)/m_Temperature))) << std::endl;
			std::cout << std::endl;

			if( m_maximize ) {		
				// if(( m_aceptance < float(1/(1+exp ((temp_eval-temp_bestScore)/m_Temperature))) )  || (temp_bestScore < temp_eval) || (temp_bestScore == temp_eval))	{
				if((temp_bestScore < temp_eval)||(temp_bestScore == temp_eval)){
					temp_bestScore = temp_eval; 
					for(int m=0; m<(m_numberofparameters) ; m++) temp_old_position[m] = temp_new_position[m];

if( has_output ) {
	temp_string.clear();
	//temp_string.c_str("");
	for(int r = 0; r<m_numberofparameters; r++) {
		ss.str("");
		ss << temp_new_position[r];
		
		temp_string = temp_string + ss.str() + ",";
	}
	
	ss.str("");
	ss << temp_eval;
	temp_string = temp_string + ss.str() + "\n";
	
	file << temp_string; // .c_str();
}

					break;
				} else {
					j++;
					std::cout << "Numero de intentos: " << j << std::endl;
				}						
			}  else {	
				//if(( m_aceptance < float( 1/(1+exp ((temp_bestScore - temp_eval)/m_Temperature))) ) || (temp_bestScore > temp_eval) || (temp_bestScore == temp_eval))	{
				if((temp_bestScore > temp_eval) || (temp_bestScore == temp_eval)){
					temp_bestScore = temp_eval; 
					for(int m=0; m<(m_numberofparameters) ; m++) temp_old_position[m] = temp_new_position[m];

if( has_output ) {
	temp_string.clear();
	//temp_string.c_str("");
	for(int r = 0; r<m_numberofparameters; r++) {
		ss.str("");
		ss << temp_new_position[r];
		
		temp_string = temp_string + ss.str() + ",";
	}
	
	ss.str("");
	ss << temp_eval;
	temp_string = temp_string + ss.str() + "\n";
	
	file << temp_string; // .c_str();
}

					break;
				} else {
					j++;
					std::cout << "Numero de intentos: " << j << std::endl;
				}
			}
		}

		// Criterio de parada :: si no mejora lo suficiente o se pasan todos los miembros...
		if( m_maximize ) {	
			if(temp_bestScore > bestScore) bestScore = temp_bestScore;
			else if( j==max_j ) cont = false;
		} 
		else { 
			if(temp_bestScore < bestScore) bestScore = temp_bestScore;
			else if( j==max_j) cont=false;
		}
		
		i++;
	}



if( has_output ) {
	temp_string.clear();
	//temp_string.c_str("");
	for(int r = 0; r<m_numberofparameters; r++) {
		ss.str("");
		ss <<temp_old_position[r];
		
		temp_string = temp_string + ss.str() + ",";
	}
	
	ss.str("");
	ss << bestScore;
	temp_string = temp_string + ss.str() + "\n";
	
	file << temp_string; // .c_str();

	 file.close();
}
	

	m_intensityBased->SetTransformation( temp_old_position );
		temp_eval = m_intensityBased->Update();
			if(temp_eval==0.0) 	return 1;

	
	//std::string finalImage = m_outputdir+"\\"+m_position+"_finalImage.nrrd";
	std::string finalImage = m_outputdir+"\\finalImage.nrrd";
		writer->SetInput( m_intensityBased->GetOutputImage() );
		writer->SetFileName( finalImage );
		try{
			writer->Update();
		}catch( itk::ExceptionObject & excp ) {
			std::cout << " Simulated Annealing, Final Image Writer exception! " << std::endl;
			std::cout << excp << std::endl;
			std::cout << std::endl;
			return 1;
		}

	std::cout << std::endl;
	std::cout << "Final Image Writed !"<< std::endl;
	std::cout << std::endl;

	clock_t en = clock();
	double elapse_time = double(en-bg)/CLOCKS_PER_SEC;
	std::cout << std::endl;
	std::cout << "Time : " << elapse_time << std::endl;
	std::cout << std::endl;

	return 0;
}

#endif