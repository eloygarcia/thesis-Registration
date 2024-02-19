#include "cuda_runtime.h"
#include "device_launch_parameters.h"


__device__ void computealpha( const int  ind, const float origen, const float voxelSize,const float p1, const float p2, float &alpha) {	// alpha = [ X_plane(i) - X1] / [X2-X1] ;
		alpha = ((ind*voxelSize + origen)-p1)/(p2-p1); };

__device__ void Determinante( float* v0, float* v1, float* v2, float* v3, float & det )
{
	//float det = 0.0;
	det  = (v0[2]*v1[1]*v2[0]*1 - v0[1]*v1[2]*v2[0]*1 - 
			v0[2]*v1[0]*v2[1]*1 + v0[0]*v1[2]*v2[1]*1 + 
			v0[1]*v1[0]*v2[2]*1 - v0[0]*v1[1]*v2[2]*1 - 
			v0[2]*v1[1]*1*v3[0] + v0[1]*v1[2]*1*v3[0] + 
			v0[2]*1*v2[1]*v3[0] - 1*v1[2]*v2[1]*v3[0] - 
			v0[1]*1*v2[2]*v3[0] + 1*v1[1]*v2[2]*v3[0] + 
			v0[2]*v1[0]*1*v3[1] - v0[0]*v1[2]*1*v3[1] - 
			v0[2]*1*v2[0]*v3[1] + 1*v1[2]*v2[0]*v3[1] + 
			v0[0]*1*v2[2]*v3[1] - 1*v1[0]*v2[2]*v3[1] - 
			v0[1]*v1[0]*1*v3[2] + v0[0]*v1[1]*1*v3[2] + 
			v0[1]*1*v2[0]*v3[2] - 1*v1[1]*v2[0]*v3[2] - 
			v0[0]*1*v2[1]*v3[2] + 1*v1[0]*v2[1]*v3[2])/6;

	// return det;
}

__device__ void triLinearInterpolator(const float* imagepointer, const int* size, const float* pixel, float & value) // es un pixel continuo, con lo cual, vale
{
//	float value = 0.0;
	int x0, y0, z0, x1, y1, z1;
		x0 = (int) (floor(pixel[0])); x1=x0+1; // x1 = (int) (ceilf(pixel[0]));
		y0 = (int) (floor(pixel[1])); y1=y0+1; // y1 = (int) (ceilf(pixel[1]));
		z0 = (int) (floor(pixel[2])); z1=z0+1; // z1 = (int) (ceilf(pixel[2]));

	float xd = (pixel[0]-x0)/(x1-x0);
	float yd = (pixel[1]-y0)/(y1-y0);
	float zd = (pixel[2]-z0)/(z1-z0);

	int idx_1 = 0;		int idx_2 = 0;

		idx_1 = x0 + (size[0] * y0) + (size[0]*size[1]*z0);
		idx_2 = x1 + (size[0] * y0) + (size[0]*size[1]*z0);
	float c00 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y0) + (size[0]*size[1]*z1);
		idx_2 = x1 + (size[0] * y0) + (size[0]*size[1]*z1);
	float c01 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y1) + (size[0]*size[1]*z0);
		idx_2 = x1 + (size[0] * y1) + (size[0]*size[1]*z0);
	float c10 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );
		idx_1 = x0 + (size[0] * y1) + (size[0]*size[1]*z1);
		idx_2 = x1 + (size[0] * y1) + (size[0]*size[1]*z1);
	float c11 = ( imagepointer[idx_1] * (1-xd) ) + ( imagepointer[ idx_2 ] * (xd) );

		float c0 = c00 * (1-yd) + c10 * yd;
		float c1 = c01 * (1-yd) + c11 * yd;

		float c = c0 * (1-zd) + c1 * zd;

		value = c;
}

__device__ void computeBaricentricCoordinates(float* posicionPunto, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* baricentricCoordinates, bool & is_inside)
{
	baricentricCoordinates[0] = 0.0;
	baricentricCoordinates[1] = 0.0;
	baricentricCoordinates[2] = 0.0;
	baricentricCoordinates[3] = 0.0;

	//float V = 0;
	float V = 0;
	Determinante( vertex_0, vertex_1, vertex_2, vertex_3, V );

	//float v1 = 0; 
	float v1 = 0;
	Determinante( posicionPunto, vertex_1, vertex_2, vertex_3, v1);
	//float v2 = 0;
	float v2 = 0;
	Determinante( vertex_0, posicionPunto, vertex_2, vertex_3, v2);
	//float v3 = 0;
	float v3 = 0;
	Determinante( vertex_0, vertex_1, posicionPunto, vertex_3, v3);
	//float v4 = 0;
	float v4 = 0;
	Determinante( vertex_0, vertex_1, vertex_2 ,posicionPunto, v4);
		
	//bool is_inside = false;
	is_inside = false;

	if( ((v1/V)>= -.01 && (v1/V)<= 1.01) && ((v2/V)>= -0.01 && (v2/V)<= 1.01) && ((v3/V)>= -0.01 && (v3/V)<=1.01) && ((v4/V)>= -0.01 && (v4/V)<= 1.01) )
	//if( ((v1/V)>=0 && (v1/V)<= 1) && ((v2/V)>=0 && (v2/V)<= 1) && ((v3/V)>=0 && (v3/V)<=1) && ((v4/V)>=0 && (v4/V)<=1) )
	{
		/*
		std::cout << std::endl;
		std::cout << " punto : [" << posicionPunto[0] << ", " << posicionPunto[1] << ", " << posicionPunto[2] << "] " << std::endl;
		std::cout << "vertex_0 : [" << vertex_0[0] << ", " << vertex_0[1] << ", " << vertex_0[2] << "] " << std::endl;
		std::cout << "vertex_1 : [" << vertex_1[0] << ", " << vertex_1[1] << ", " << vertex_1[2] << "] " << std::endl;
		std::cout << "vertex_2 : [" << vertex_2[0] << ", " << vertex_2[1] << ", " << vertex_2[2] << "] " << std::endl;
		std::cout << "vertex_3 : [" << vertex_3[0] << ", " << vertex_3[1] << ", " << vertex_3[2] << "] " << std::endl;
		std::cout << std::endl;

		std::cout << "Baricentric Coordinates : [" << v1/V <<", " << v2/V << ", " << v3/V << ", " << v4/V << "] " << std::endl;
		*/
		is_inside = true;
		/**/
		baricentricCoordinates[0] = v1/V; 
		baricentricCoordinates[1] = v2/V;
		baricentricCoordinates[2] = v3/V;
		baricentricCoordinates[3] = v4/V;

	}

	// std::cout << "Baricentric Coordinates : [" << v1/V <<", " << v2/V << ", " << v3/V << ", " << v4/V << "] " << std::endl;

	// Sleep(2000);

	// return is_inside;
}

__device__ void computeCartessianCoordinates(float* baricentricCoordinates, float* vertex_0, float* vertex_1, float* vertex_2, float* vertex_3, float* posicionPunto)
{
	posicionPunto[0] = baricentricCoordinates[0]*vertex_0[0] + baricentricCoordinates[1]*vertex_1[0] + baricentricCoordinates[2]*vertex_2[0] + baricentricCoordinates[3]*vertex_3[0];
	posicionPunto[1] = baricentricCoordinates[0]*vertex_0[1] + baricentricCoordinates[1]*vertex_1[1] + baricentricCoordinates[2]*vertex_2[1] + baricentricCoordinates[3]*vertex_3[1];
	posicionPunto[2] = baricentricCoordinates[0]*vertex_0[2] + baricentricCoordinates[1]*vertex_1[2] + baricentricCoordinates[2]*vertex_2[2] + baricentricCoordinates[3]*vertex_3[2];
}

/*
__global__ void kernel_projection(const int* dev_3d_size, const float* dev_3d_spacing, const float* dev_3d_origen, const float* dev_3d_imagepointer,
								  const int* dev_2d_size, const float* dev_2d_spacing, const float* dev_2d_origen, float* dev_2d_imagepointer,
								  const float* source)
{
	int i = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	int numberOfPixels2d = dev_2d_size[0] * dev_2d_size[1]; // *size2d[2];

	if( i<numberOfPixels2d){

	// Fila y columna de la imagen !
	int row = (int) floorf( i / dev_2d_size[0] ); 
	int col = (int) ( i - (row*dev_2d_size[0]) ); 
	

	// Posición del pixel del detector !!
	float x2 = dev_2d_origen[0] + (col * dev_2d_spacing[0] );
	float y2 = dev_2d_origen[1] + (row * dev_2d_spacing[1] );

	// Vector de dirección !!
	float vect[3];
		vect[0] = x2 - source[0];
		vect[1] = y2 - source[1];
		vect[2] = dev_2d_origen[2] - source[2];  // Revisar este punto.

	// Distancia de la fuente al detector !!
	float xa = pow(vect[0],2);
	float ya = pow(vect[1],2);
	float za = pow(vect[2],2);

	float dist12 = sqrt( xa + ya + za );

	// Cálculo del alpha en Z !!
	float temp = 0.0f;
	computealpha( 0, dev_3d_origen[2], dev_3d_spacing[2], source[2], dev_2d_origen[2], temp );

	// Resolver la ecuación de la recta !!
	float temp_dist[3] = {0.0f, 0.0f, 0.0f};
	float point[3] = {0.0f, 0.0f, 0.0f};
	float step = 0.0005f;
	float t=temp;

	float pixel[3] = {0.0f, 0.0f, 0.0f};

	float value = 0.0f;
	float length = 0.0f;

	// Longitud del paso !!
	temp_dist[0] = pow(step*vect[0],2);
	temp_dist[1] = pow(step*vect[1],2);
	temp_dist[2] = pow(step*vect[2],2);
	float l_step = sqrt(temp_dist[0] + temp_dist[1] + temp_dist[2]);

	while( t<1 )
	{
		// Calculo del siguiente punto en la recta
		point[0] = source[0] + t * vect[0];
		point[1] = source[1] + t * vect[1];
		point[2] = source[2] + t * vect[2];
		// Posición del voxel parcial !!
		pixel[0] = (point[0] - dev_3d_origen[0]) / dev_3d_spacing[0];
		pixel[1] = (point[1] - dev_3d_origen[1]) / dev_3d_spacing[1];
		pixel[2] = (point[2] - dev_3d_origen[2]) / dev_3d_spacing[2];
		// Interpolación trilienal
		if((pixel[0]<0 || pixel[0]>dev_3d_size[0]-1) || (pixel[1]<0 || pixel[1]>dev_3d_size[1]-1) || (pixel[2]<0 || pixel[2]>dev_3d_size[2]-1))  value = 0;
		else triLinearInterpolator( dev_3d_imagepointer, dev_3d_size, pixel, value);
		// Longitud acumulada !!
		length += (1000 * value * l_step);
		// Nuevo punto!! 
		t+=step;
	}

	if(length>0 & length < 65535 ) dev_2d_imagepointer[ i ] = length;
	else dev_2d_imagepointer[ i ] = 0.0f;
	
	}
}
*/
__global__ void kernel_projection(const int* dev_3d_size, const float* dev_3d_spacing, const float* dev_3d_origen, const float* dev_3d_imagepointer,
								  const float* dev_i_points, const float* dev_f_points, const int* dev_elements,
								  const float* dev_grid_origen, const float* dev_grid_spacing, const int* dev_grid_size,
								  const int* dev_flags, const int* dev_cumsum, const int* dev_correspondingElements,
								  const int* dev_2d_size, const float* dev_2d_spacing, const float* dev_2d_origen, unsigned short* dev_2d_imagepointer, //float* dev_2d_imagepointer,
								  const float* dev_source)
{
	int i = (blockDim.x * blockIdx.x) + threadIdx.x;
	
	int numberOfPixels2d = dev_2d_size[0] * dev_2d_size[1]; // *size2d[2];

//	if( i<numberOfPixels2d){
		// Fila y columna de la imagen !
		int row = (int) floorf( i / dev_2d_size[0] ); 
		int col = (int) ( i - (row*dev_2d_size[0]) ); 

		// Posición del pixel del detector !!
		float x2 = dev_2d_origen[0] + (col * dev_2d_spacing[0] );
		float y2 = dev_2d_origen[1] + (row * dev_2d_spacing[1] );

		// Vector de dirección !!
		float vect[3];
			vect[0] = x2 - dev_source[0];
			vect[1] = y2 - dev_source[1];
			vect[2] = dev_2d_origen[2] - dev_source[2];  // Revisar este punto.

		// Distancia de la fuente al detector !!
		float xa = pow(vect[0],2);
		float ya = pow(vect[1],2);
		float za = pow(vect[2],2);

		float dist12 = sqrt( xa + ya + za );

		// Cálculo del alpha en Z !!
		float temp = 0.0f;
		// En este caso el alpha se calcula en relación a la grid !!
		computealpha( 0, dev_grid_origen[2], dev_grid_spacing[2], dev_source[2], dev_2d_origen[2], temp ); 

		// Resolver la ecuación de la recta !!
		float step = 0.0005f;
		float t=temp;
		// on the grid !
		float temp_dist[3] = {0.0f, 0.0f, 0.0f};  // temporal distance
		float point[3] = {0.0f, 0.0f, 0.0f};  // physical point on the grid
		
		int voxel_grid[3] = {0,0,0};
		float baricentricCoordinates[4] = {0.0f, 0.0f, 0.0f, 0.0f};

		int index = 0;

		int index_start = 0;
		int number_of_elements_here = 0;
		// int count = 0;

		int element_index_number = 0;
		int element_number = 0;

		int numberofVoxelsGrid = dev_grid_size[0] * dev_grid_size[1] * dev_grid_size[2];

		float vertex_0[3] = {0.0f,0.0f,0.0f};
		float vertex_1[3] = {0.0f,0.0f,0.0f};
		float vertex_2[3] = {0.0f,0.0f,0.0f};
		float vertex_3[3] = {0.0f,0.0f,0.0f};

		bool is_inside = false;

		//float old_vertex_0[3] = {0.0f,0.0f,0.0f};
		//float old_vertex_1[3] = {0.0f,0.0f,0.0f};
		//float old_vertex_2[3] = {0.0f,0.0f,0.0f};
		//float old_vertex_3[3] = {0.0f,0.0f,0.0f};

		// on the mri
		float position_mri[3] = {0.0f, 0.0f, 0.0f};
		float pixel_mri[3] = {0.0f, 0.0f, 0.0f}; // el pixel es el voxel de a grid !! // No lo necesito?'
		
		// Solving the step !
		float value = 0.0f;
		float length = 0.0f;

		// Longitud del paso !!
		temp_dist[0] = pow(step*vect[0],2);
		temp_dist[1] = pow(step*vect[1],2);
		temp_dist[2] = pow(step*vect[2],2);
		float l_step = sqrt(temp_dist[0] + temp_dist[1] + temp_dist[2]);

		int max_num_of_elem_here = 0;

		while( t<1 )
		{
			// Calculo del siguiente punto en la recta
			point[0] = dev_source[0] + t * vect[0];
			point[1] = dev_source[1] + t * vect[1];
			point[2] = dev_source[2] + t * vect[2];
			
			// Voxel on the grid //
			voxel_grid[0] = (int)(floorf((point[0] - dev_grid_origen[0])/dev_grid_spacing[0]));
			voxel_grid[1] = (int)(floorf((point[1] - dev_grid_origen[1])/dev_grid_spacing[1]));
			voxel_grid[2] = (int)(floorf((point[2] - dev_grid_origen[2])/dev_grid_spacing[2]));

			if( ((voxel_grid[0]>0) && (voxel_grid[0]<dev_grid_size[0]) ) &&
				((voxel_grid[1]>0) && (voxel_grid[1]<dev_grid_size[1]) ) &&
				((voxel_grid[2]>0) && (voxel_grid[2]<dev_grid_size[2]) ) )
			{
				index =(int)( (voxel_grid[2]*(dev_grid_size[0]*dev_grid_size[1])) + (voxel_grid[1]*dev_grid_size[0]) + voxel_grid[0]);
				
				index_start = (int)dev_cumsum[ index ];
				number_of_elements_here = (int) dev_flags[ index ];
				
				if(number_of_elements_here!=0){
					for( int j=0; j<number_of_elements_here; j++){
						element_index_number = index_start + j;
						if( element_index_number < dev_cumsum[ numberofVoxelsGrid -1] ){
							element_number =  dev_correspondingElements[ element_index_number ];

							// float temp_vertex_0[3] = {0.0,0.0,0.0};
							vertex_0[0] = dev_f_points[ 3*dev_elements[ 4*element_number ]  ];
							vertex_0[1] = dev_f_points[ 3*dev_elements[ 4*element_number ] +1 ];
							vertex_0[2] = dev_f_points[ 3*dev_elements[ 4*element_number ] +2 ];
					
							// float temp_vertex_1[3] = {0.0,0.0,0.0};
							vertex_1[0] = dev_f_points[ 3*dev_elements[ 4*element_number +1 ] ];
							vertex_1[1] = dev_f_points[ 3*dev_elements[ 4*element_number +1 ] +1 ];
							vertex_1[2] = dev_f_points[ 3*dev_elements[ 4*element_number +1 ] +2 ];
				
							// float temp_vertex_2[3] = {0.0,0.0,0.0};
							vertex_2[0] = dev_f_points[ 3*dev_elements[ 4*element_number +2 ] ];
							vertex_2[1] = dev_f_points[ 3*dev_elements[ 4*element_number +2 ] +1 ];
							vertex_2[2] = dev_f_points[ 3*dev_elements[ 4*element_number +2 ] +2 ];
				
							// float temp_vertex_3[3] = {0.0,0.0,0.0};
							vertex_3[0] = dev_f_points[ 3*dev_elements[ 4*element_number +3 ] ];
							vertex_3[1] = dev_f_points[ 3*dev_elements[ 4*element_number +3 ] +1 ];
							vertex_3[2] = dev_f_points[ 3*dev_elements[ 4*element_number +3 ] +2 ];

							computeBaricentricCoordinates(point, vertex_0, vertex_1, vertex_2, vertex_3, baricentricCoordinates, is_inside);
				
							if( is_inside ){
								vertex_0[0] = dev_i_points[ 3*dev_elements[ 4*element_number ]  ];
								vertex_0[1] = dev_i_points[ 3*dev_elements[ 4*element_number ] +1 ];
								vertex_0[2] = dev_i_points[ 3*dev_elements[ 4*element_number ] +2 ];
					
								// float temp_vertex_1[3] = {0.0,0.0,0.0};
								vertex_1[0] = dev_i_points[ 3*dev_elements[ 4*element_number +1 ] ];
								vertex_1[1] = dev_i_points[ 3*dev_elements[ 4*element_number +1 ] +1 ];
								vertex_1[2] = dev_i_points[ 3*dev_elements[ 4*element_number +1 ] +2 ];
					
								// float temp_vertex_2[3] = {0.0,0.0,0.0};
								vertex_2[0] = dev_i_points[ 3*dev_elements[ 4*element_number +2 ] ];
								vertex_2[1] = dev_i_points[ 3*dev_elements[ 4*element_number +2 ] +1 ];
								vertex_2[2] = dev_i_points[ 3*dev_elements[ 4*element_number +2 ] +2 ];
				
								// float temp_vertex_3[3] = {0.0,0.0,0.0};
								vertex_3[0] = dev_i_points[ 3*dev_elements[ 4*element_number +3 ] ];
								vertex_3[1] = dev_i_points[ 3*dev_elements[ 4*element_number +3 ] +1 ];
								vertex_3[2] = dev_i_points[ 3*dev_elements[ 4*element_number +3 ] +2 ];
		
								computeCartessianCoordinates( baricentricCoordinates, vertex_0, vertex_1, vertex_2, vertex_3, position_mri);

								// Posición del voxel parcial !! // esto tiene que ser on the grid
								pixel_mri[0] = (float)((position_mri[0] - dev_3d_origen[0]) / dev_3d_spacing[0]);
								pixel_mri[1] = (float)((position_mri[1] - dev_3d_origen[1]) / dev_3d_spacing[1]);
								pixel_mri[2] = (float)((position_mri[2] - dev_3d_origen[2]) / dev_3d_spacing[2]);
									
								// Interpolación trilineal
								if(((pixel_mri[0]>0.0f) && (pixel_mri[0]<(float)(dev_3d_size[0]))) &&
								   ((pixel_mri[1]>0.0f) && (pixel_mri[1]<(float)(dev_3d_size[1]))) &&
								   ((pixel_mri[2]>0.0f) && (pixel_mri[2]<(float)(dev_3d_size[2])))) 
									{ 
									//	value = 1.0f;
										triLinearInterpolator( dev_3d_imagepointer, dev_3d_size, pixel_mri, value);
										//max_num_of_elem_here = 0.0f;
								} else {
									value = 0.0f;
									//max_num_of_elem_here = 1.0f;
								}
	
								// Longitud acumulada !!
								length += (1000 * value * l_step);
								break;
							}
						}
					}
				}
			}
			// Nuevo Punto !!
			t+=step;
		}

		if(length>0 & length < 65535 ) dev_2d_imagepointer[ i ] = (unsigned short) length;
		else dev_2d_imagepointer[ i ] = 0;
		//else dev_2d_imagepointer[ i ] = 0.0f;

		//dev_2d_imagepointer[ i ] = length;
		
}


__global__ void fill_dos( float * imagepointer)
{
	int i = blockDim.x * blockIdx.x +threadIdx.x;
	imagepointer[i]=2.0f;
}