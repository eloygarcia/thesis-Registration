#include "CudaProjection.h"
#include "CudaLibrary_kernels.cu"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Constructor y Destructor

CudaProjection::CudaProjection()
{
	// mri
	m_mri_size = new int[3];
	m_mri_spacing = new float[3];
	m_mri_origen = new float[3];
//	m_mri_imagepointer = 0;

	// mesh and grid !
//	m_i_points = 0;
//	m_f_points = 0;
//	m_elements = 0;

	m_grid_origen = new float[3];
	m_grid_spacing = new float[3];
	m_grid_size = new int[3];
	
//	m_flags = 0;
//	m_cumsum = 0;
//	m_correspondingElements = 0;

	// mamo simulada
	m_simulada_size = new int[3];
	m_simulada_spacing = new float[3];
	m_simulada_origen = new float[3];
	//m_2d_imagepointer = = new float[ numberOfPixels_2D ];

	m_source = new float[3];

	// Variables en:  device
//	dev_mri_size = new int[3];
//	dev_mri_spacing = new float[3];
//	dev_mri_origen = new float[3];
//numberOfPixels_3D = 0 ;
//	dev_3d_imagepointer = new float[ numberOfPixels_3D ];

//	dev_simulada_size = new int[3];
//	dev_simulada_spacing = new float[3];
//	dev_simulada_origen = new float[3]; 
//numberOfPixels_2D = ( m_2d_size[0] * m_2d_size[1] );
//	dev_2d_imagepointer = new float[ numberOfPixels_2D ];

//	dev_source = new float[3];

	cudaError_t cudaStatus;
}

CudaProjection::~CudaProjection()
{
	// delete [] m_parameters;
	// host
	delete[] m_mri_size;
	delete[] m_mri_spacing;
	delete[] m_mri_origen;
	delete[] m_mri_imagepointer;

	delete [] m_i_points;
	delete [] m_f_points;
	delete [] m_elements;

	delete [] m_grid_origen;
	delete [] m_grid_size;
	delete [] m_grid_spacing;

	delete [] m_flags;
	delete [] m_cumsum;
	delete [] m_correspondingElements;

	delete[] m_simulada_size;
	delete[] m_simulada_spacing;
	delete[] m_simulada_origen;
//	delete[] m_simulada_imagepointer;

	delete[] m_source;

	// device
//	delete[] dev_mri_size;
//	delete[] dev_mri_spacing;
//	delete[] dev_mri_origen;
//	delete[] dev_mri_imagepointer;
//
//	delete [] dev_i_points;
//	delete [] dev_f_points;
//	delete [] dev_elements;
//
//	delete [] dev_grid_origen;
//	delete [] dev_grid_size;
//	delete [] dev_grid_spacing;
//
//	delete [] dev_flags;
//	delete [] dev_cumsum;
//	delete [] dev_correspondingElements;
//
//	delete[] dev_simulada_size;
//	delete[] dev_simulada_spacing;
//	delete[] dev_simulada_origen;
//	delete[] dev_simulada_imagepointer;
//
//	delete[] dev_source;
}

void CudaProjection::Initialize()
{
	m_mri_size[0] = m_parameters->mri_size[0];
	m_mri_size[1] = m_parameters->mri_size[1];
	m_mri_size[2] = m_parameters->mri_size[2];
numberOfPixels_MRI = m_mri_size[0] * m_mri_size[1] * m_mri_size[2];

	m_mri_spacing[0] = m_parameters->mri_spacing[0];
	m_mri_spacing[1] = m_parameters->mri_spacing[1];
	m_mri_spacing[2] = m_parameters->mri_spacing[2];

	m_mri_origen[0] = m_parameters->mri_origen[0];
	m_mri_origen[1] = m_parameters->mri_origen[1];
	m_mri_origen[2] = m_parameters->mri_origen[2];

	m_mri_imagepointer = new float[ numberOfPixels_MRI ];
		for(int i=0; i<numberOfPixels_MRI; i++) m_mri_imagepointer[i] = m_parameters->mri_imagePointer[i];

numberOfPoints = m_parameters->numberOfPoints;
	m_i_points = new float[ 3*numberOfPoints ];
		for( int i=0; i<3*numberOfPoints; i++) m_i_points[i] = m_parameters->initial_points[i];
	m_f_points = new float[ 3*numberOfPoints ];
		for(int i=0; i<3*numberOfPoints; i++) m_f_points[i] = m_parameters->final_points[i];
numberOfElements = m_parameters->numberOfElements;
	m_elements = new int[4*numberOfElements];	
	for(int i=0; i<4*numberOfElements; i++) m_elements[i] = m_parameters->elements[i];

	m_grid_origen[0] = m_parameters->grid_origen[0];
	m_grid_origen[1] = m_parameters->grid_origen[1];
	m_grid_origen[2] = m_parameters->grid_origen[2];

	m_grid_spacing[0] = m_parameters->grid_spacing[0];
	m_grid_spacing[1] = m_parameters->grid_spacing[1];
	m_grid_spacing[2] = m_parameters->grid_spacing[2];

	m_grid_size[0] = m_parameters->grid_size[0];
	m_grid_size[1] = m_parameters->grid_size[1];
	m_grid_size[2] = m_parameters->grid_size[2];
numberOfVoxelsGrid = m_grid_size[0] * m_grid_size[1] * m_grid_size[2] ; 

	m_flags = new int[ numberOfVoxelsGrid ];
		for(int i=0; i<numberOfVoxelsGrid; i++) m_flags[i] = m_parameters->flags[i];
	m_cumsum = new int[ numberOfVoxelsGrid ];
		for(int i=0; i<numberOfVoxelsGrid; i++) m_cumsum[i] = m_parameters->cumsum[i];
maximumCorrespondingElements = m_cumsum[ numberOfVoxelsGrid -1];
	m_correspondingElements = new int[ maximumCorrespondingElements];
	for (int i=0; i< maximumCorrespondingElements; i++)	m_correspondingElements[i] = m_parameters->correspondingElements[i];

	m_simulada_size[0] = m_parameters->mamo_size[0];
	m_simulada_size[1] = m_parameters->mamo_size[1];
	m_simulada_size[2] = m_parameters->mamo_size[2];
numberOfPixels_Simulada = m_simulada_size[0] * m_simulada_size[1];

	m_simulada_origen[0] = m_parameters->mamo_origen[0];
	m_simulada_origen[1] = m_parameters->mamo_origen[1];
	m_simulada_origen[2] = m_parameters->mamo_origen[2];

	m_simulada_spacing[0] = m_parameters->mamo_spacing[0];
	m_simulada_spacing[1] = m_parameters->mamo_spacing[1];
	m_simulada_spacing[2] = m_parameters->mamo_spacing[2];
	
	//m_simulada_imagepointer = new float[numberOfPixels_Simulada];
	m_simulada_imagepointer = new unsigned short[numberOfPixels_Simulada];

	m_source[0] = m_parameters->source[0];
	m_source[1] = m_parameters->source[1];
	m_source[2] = m_parameters->source[2];
}

// Metodos
void CudaProjection::Update()
{
	/*
	printf("\n");
	printf("Entra en cuda\n");
	printf("\n");

// Timer !
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start,0);
	*/

	Initialize();

// Allocacion de la memoria GPU !
	// MRI !
    cudaStatus = cudaMalloc((void**)&dev_mri_size, 3*sizeof(int));
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 3d Size!\n");   // MRI Size
	cudaStatus = cudaMalloc((void**)&dev_mri_spacing, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 3d Spacing!\n"); // MRI Spacing
	cudaStatus = cudaMalloc((void**)&dev_mri_origen, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 3d Origen!\n");  // MRI Origen

	cudaStatus = cudaMalloc((void**)&dev_mri_imagepointer, numberOfPixels_MRI*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 3d image pointer!\n");  // MRI imagen !!

	// Mesh & Grid !!
	cudaStatus = cudaMalloc((void**)&dev_i_points, 3*numberOfPoints*sizeof(float)); 
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc i_point !\n"); // i_points
	cudaStatus = cudaMalloc((void**)&dev_f_points, 3*numberOfPoints*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc f_point !\n"); // f_points

	cudaStatus = cudaMalloc((void**)&dev_elements, 4*numberOfElements*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc elements !\n"); // elements;

	cudaStatus = cudaMalloc((void**)&dev_grid_origen, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Grid Origen !\n"); // grid_origen
	cudaStatus = cudaMalloc((void**)&dev_grid_spacing, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Grid spacing!\n"); // grid_Spacing
	cudaStatus = cudaMalloc((void**)&dev_grid_size, 3*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Grid Size !\n"); // grid_Size

	cudaStatus = cudaMalloc((void**)&dev_flags, numberOfVoxelsGrid*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Flags !\n"); // flags
	cudaStatus = cudaMalloc((void**)&dev_cumsum, numberOfVoxelsGrid*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc CumSum !\n"); // cumsum
	cudaStatus = cudaMalloc((void**)&dev_correspondingElements, maximumCorrespondingElements*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc corresponding element !\n"); // corresponding Elements !

	// Imagen Simulada !
	cudaStatus = cudaMalloc((void**)&dev_simulada_size, 3*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 2d Size!\n");  // Simulada Size
	cudaStatus = cudaMalloc((void**)&dev_simulada_spacing, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 2d Spacing!\n");  // Simulada Spacing
	cudaStatus = cudaMalloc((void**)&dev_simulada_origen, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 2d Origen!\n");  // simulada Origen

	cudaStatus = cudaMalloc((void**)&dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(unsigned short));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 2d image pointer!\n");  // Simulada Imagen !!
	cudaStatus = cudaMemset((void*)dev_simulada_imagepointer, 0, numberOfPixels_Simulada*sizeof(unsigned short));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc 2d image pointer!\n");  // Inicialización de la imagen simulada a Zeros !!

	// Source !
	cudaStatus = cudaMalloc((void**)&dev_source, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc Source!\n");  // source !"


// Copia a memoria device
	// mri
	cudaStatus = cudaMemcpy(dev_mri_size, (const int*) m_mri_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	// cudaStatus = cudaMemcpy(dev_mri_size, (const int*) m_mri_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Size Host2Dev!\n");
	cudaStatus = cudaMemcpy(dev_mri_spacing, (const float*) m_mri_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
	// cudaStatus = cudaMemcpy(dev_mri_spacing, this->m_mri_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Spacing Host2Dev!\n");
	cudaStatus = cudaMemcpy(dev_mri_origen, (const float*) m_mri_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
	// cudaStatus = cudaMemcpy(dev_mri_origen, this->m_mri_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Origen Host2Dev!\n");
	cudaStatus = cudaMemcpy(dev_mri_imagepointer, (const float*) m_mri_imagepointer, numberOfPixels_MRI*sizeof(float), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_mri_imagepointer, this->m_mri_imagepointer, numberOfPixels_MRI*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d image pointer Host2Dev!\n");

	// Mesh & Grid !!
	cudaStatus = cudaMemcpy(dev_i_points, (const float*) m_i_points, 3*numberOfPoints*sizeof(float), cudaMemcpyHostToDevice); 
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy i_point Host2Dev !\n"); // i_points
	cudaStatus = cudaMemcpy(dev_f_points, (const float*) m_f_points,  3*numberOfPoints*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy f_point Host2Dev !\n"); // f_points

	cudaStatus = cudaMemcpy(dev_elements, (const int*) m_elements, 4*numberOfElements*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy elements Host2Dev !\n"); // elements;

	cudaStatus = cudaMemcpy(dev_grid_origen, (const float*) m_grid_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy Grid Origen !\n"); // grid_origen
	cudaStatus = cudaMemcpy(dev_grid_spacing, (const float*) m_grid_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy Grid spacing!\n"); // grid_Spacing
	cudaStatus = cudaMemcpy(dev_grid_size, (const int*) m_grid_size, 3*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy Grid Size !\n"); // grid_Size

	cudaStatus = cudaMemcpy(dev_flags, (const int*) m_flags, numberOfVoxelsGrid*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcpy Flags !\n"); // flags
	cudaStatus = cudaMemcpy(dev_cumsum, (const int*) m_cumsum, numberOfVoxelsGrid*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcpy CumSum !\n"); // cumsum
	cudaStatus = cudaMemcpy(dev_correspondingElements, (const int*) m_correspondingElements, maximumCorrespondingElements*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcpy corresponding element !\n"); // corresponding Elements !

	// mamo simulada
	cudaStatus = cudaMemcpy(dev_simulada_size, (const int*) m_simulada_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_simulada_size, this->m_simulada_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d Size Host2Dev!\n");
	cudaStatus = cudaMemcpy(dev_simulada_spacing, (const float*) m_simulada_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_simulada_spacing, this->m_simulada_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d Spacing Host2Dev!\n");
	cudaStatus = cudaMemcpy(dev_simulada_origen, (const float*) m_simulada_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_simulada_origen, this->m_simulada_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d Origen Host2Dev!\n");
//	cudaStatus = cudaMemcpy(dev_2d_imagepointer, this->m_2d_imagepointer, numberOfPixels_2D*sizeof(float), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d image pointer Host2Dev!\n");

	cudaStatus = cudaMemcpy(dev_source, (const float*) m_source, 3*sizeof(float), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_source, this->m_source, 3*sizeof(float), cudaMemcpyHostToDevice);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Source Host2Dev!\n");

/*
		// Checkin memory !!
		printf("\n");
		printf("3d size: [ %i,%i,%i]\n", m_3d_size[0],m_3d_size[1],m_3d_size[2]);
		printf("3D Spacing : [ %f,%f,%f]\n",m_3d_spacing[0], m_3d_spacing[1], m_3d_spacing[2]);
		printf("3D Origen : [ %f,%f,%f]\n",m_3d_origen[0], m_3d_origen[1], m_3d_origen[2]);
		
		printf("2D Size : [ %d,%d,%d]\n",m_2d_size[0], m_2d_size[1], m_2d_size[2]);
		printf("2D Spacing : [ %f,%f,%f]\n",m_2d_spacing[0], m_2d_spacing[1], m_2d_spacing[2]);
		printf("2D Origen : [ %f,%f,%f]\n",m_2d_origen[0], m_2d_origen[1], m_2d_origen[2]);
		
		printf("Source : [ %f,%f,%f]\n",m_source[0], m_source[1], m_source[2]);
		printf("\n");
*/

		bl = (int)ceilf((float)(numberOfPixels_Simulada/128))+1;
// Kernel de proyección?
//	 fill_dos <<< bl,128 >>> (dev_simulada_imagepointer);

//	printf("Entra en el kernel\n" );

	//cudaStatus = cudaSetDevice(0);
	kernel_projection <<< bl,128 >>> (dev_mri_size, dev_mri_spacing, dev_mri_origen, dev_mri_imagepointer,
										dev_i_points, dev_f_points, dev_elements,
										dev_grid_origen, dev_grid_spacing, dev_grid_size,
										dev_flags, dev_cumsum, dev_correspondingElements,
									  dev_simulada_size, dev_simulada_spacing, dev_simulada_origen, dev_simulada_imagepointer,
									  dev_source);								   
//	printf("Sale del kernel\n" );							   
	cudaDeviceSynchronize();
// Copia a memoria host
/*	int temp_3dsize[3] = {0,0,0};
	cudaStatus = cudaMemcpy( temp_3dsize, (const int*) dev_3d_size, 3*sizeof(int),  cudaMemcpyDeviceToHost);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Size Dev2Host!\n");
		else m_3d_size = temp_3dsize;

	float temp_3dspacing[3] = {0.0, 0.0,0.0};
	cudaStatus = cudaMemcpy( temp_3dspacing, (const float*) dev_3d_spacing, 3*sizeof(float),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Spacing Dev2Host!\n");
		else m_3d_spacing = temp_3dspacing;

	float temp_3dorigen[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_3dorigen, (const float*) dev_3d_origen, 3*sizeof(float),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d Origen Dev2Host!\n");
		else m_3d_origen = temp_3dorigen;

	float* image3dpointer = new float[numberOfPixels_3D];
	cudaStatus = cudaMemcpy( image3dpointer, (const float*) dev_3d_imagepointer, numberOfPixels_3D*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 3d image pointer Dev2Host!\n");
		else m_3d_imagepointer = image3dpointer;

	int temp_2dsize[3] = {0,0,0};
	cudaStatus = cudaMemcpy( temp_2dsize, (const int*) dev_2d_size, 3*sizeof(int), cudaMemcpyDeviceToHost);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d Size Dev2Host!\n");
		else m_2d_size = temp_2dsize;

	float temp_2dspacing[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_2dspacing, (const float*) dev_2d_spacing, 3*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d Spacing Dev2Host!\n");
		else m_2d_spacing = temp_2dspacing;

	float temp_2dorigen[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_2dorigen, (const float*) dev_2d_origen, 3*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d Origen Dev2Host!\n");
		else m_2d_origen = temp_2dorigen;
*/
//	printf("va a inicializar con pixels \n" );
	//float * temp_image2dpointer;
	 //float * temp_image2dpointer = new float[numberOfPixels_Simulada]; // El m_... no convenció quizá porque no esta ba inicializado...

	//cudaStatus = cudaMemcpy(temp_image2dpointer, (const float*) dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(float),  cudaMemcpyDeviceToHost);
	cudaStatus = cudaMemcpy(m_simulada_imagepointer, (const unsigned short*)dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(unsigned short),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d image pointer Dev2Host!\n");
		else  m_parameters->simulada_imagePointer = m_simulada_imagepointer;
		

/*  AQUI VA LA RECUPERACION DE LA IMAGEN ORIGINAL !! RECUERDALO PORQUE ESTO HABRA QUE CAMBIARLO

	printf("va a inicializar con pixels \n" );
	//float * temp_image2dpointer;
	 float * temp_image2dpointer = new float[numberOfPixels_Simulada]; // El m_... no convenció quizá porque no esta ba inicializado...
	//cudaStatus = cudaMemcpy(temp_image2dpointer, (const float*) dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(float),  cudaMemcpyDeviceToHost);
	cudaStatus = cudaMemcpy(temp_image2dpointer, (const float*)dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(float),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy 2d image pointer Dev2Host!\n");
		else {
			m_simulada_imagepointer = temp_image2dpointer;
			m_parameters->simulada_imagePointer = temp_image2dpointer;
		}


*/


/*	float temp_source[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_source, (const float*) dev_source, 3*sizeof(float), cudaMemcpyDeviceToHost);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy Source Dev2Host!\n");
		else m_source = temp_source;
*/

/*
	printf("\n");
	printf("3d size: [ %i,%i,%i]\n", m_3d_size[0],m_3d_size[1],m_3d_size[2]);
	printf("3D Spacing : [ %f,%f,%f]\n",m_3d_spacing[0], m_3d_spacing[1], m_3d_spacing[2]);
	printf("3D Origen : [ %f,%f,%f]\n",m_3d_origen[0], m_3d_origen[1], m_3d_origen[2]);
	printf("3D Image Pointer : [ %f,%f,%f, ...]\n",m_3d_imagepointer[0], m_3d_imagepointer[1], m_3d_imagepointer[2]);


	printf("2D Size : [ %d,%d,%d]\n",m_2d_size[0], m_2d_size[1], m_2d_size[2]);
	printf("2D Spacing : [ %f,%f,%f]\n",m_2d_spacing[0], m_2d_spacing[1], m_2d_spacing[2]);
	printf("2D Origen : [ %f,%f,%f]\n",m_2d_origen[0], m_2d_origen[1], m_2d_origen[2]);
	printf("2D Image Pointer : [ %f,%f,%f, ...]\n",m_2d_imagepointer[0], m_2d_imagepointer[1], m_2d_imagepointer[2]);

	printf("Source : [ %f,%f,%f]\n",m_source[0], m_source[1], m_source[2]);
	printf("\n");
*/

// Liberando memoria !!
	cudaFree( (void*) dev_mri_size);				//cudaFree( temp_3dsize);
	cudaFree( (void*) dev_mri_spacing);			//cudaFree( temp_3dspacing);
	cudaFree( (void*) dev_mri_origen );			//cudaFree( temp_3dorigen);
	cudaFree( (void*) dev_mri_imagepointer );	//cudaFree( image3dpointer);

	cudaFree( (void*) dev_i_points );
	cudaFree( (void*) dev_f_points );
	cudaFree( (void*) dev_elements );

	cudaFree( (void*) dev_grid_origen );
	cudaFree( (void*) dev_grid_spacing );
	cudaFree( (void*) dev_grid_size );

	cudaFree( (void*) dev_flags );
	cudaFree( (void*) dev_cumsum );
	cudaFree( (void*) dev_correspondingElements );

	cudaFree( (void*) dev_simulada_size);				//cudaFree( temp_2dsize);
	cudaFree( (void*) dev_simulada_spacing);			//cudaFree( temp_2dspacing);
	cudaFree( (void*) dev_simulada_origen);			//cudaFree( temp_2dorigen);
	cudaFree( (void*) dev_simulada_imagepointer);		//cudaFree( (void*) temp_image2dpointer);

	//cudaFree( numberOfPixels_3D);		cudaFree( numberOfPixels_2D);

	cudaFree( (void*) dev_source);				//cudaFree( temp_source);
	cudaDeviceReset();

	cudaFree( kernel_projection );

// Time !!
/*	cudaEventRecord(stop,0);
	//cudaEventSynchronize( stop);
	cudaEventElapsedTime( &time, start, stop);
	printf( "Time: %f ms.\n", time);

	printf("\n");
	printf("Sale de cuda\n");
	printf("\n");
	*/
}

