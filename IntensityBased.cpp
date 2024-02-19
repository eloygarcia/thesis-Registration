#ifndef __intensityBased_cpp
#define __intensityBased_cpp

#include "IntensityBased.h"

IntensityBased::IntensityBased()
{
	registration = new registrationParameters();
	//transformation = new Transformations();

	m_source = new float[3];
		m_source [0] = 0.0;
		m_source [1] = 0.0;
		m_source [2] = 0.0;
	 
	m_switch = 0;

	val=0;
	value = 0;

	// ======= De: Update() ============
	size = new int[3];
		size[0]=1;
		size[1]=1;
		size[2]=1;
	origen = new float[3];
		origen[0] = 0.0;
		origen[1] = 0.0;
		origen[2] = 0.0;
	spacing = new float[3];
		spacing[0] = 1.0;
		spacing[1] = 1.0;
		spacing[2] = 1.0;

	// ============= De: Compute Regular Grid ========================
#ifdef _GPU_
	// CDProjection = new CudaProjection;
	// Cp = new CudaProjection();
#endif

	importfilter = ImportFilterType::New();

	bounding_box = new float[6];
	for(unsigned int i=0; i<6; i++){
		bounding_box[i]=0;
	}

	 importfilter = ImportFilterType::New();
	// Write
	writer2del = Writer2DLType::New();

	nuPoisson = 0.495;
	new_transformation = new float[10];
	
}

IntensityBased::~IntensityBased()
{
#ifdef _GPU_
	// delete [] CDProjection;
	// delete [] Cp;
#endif

	delete []  m_source;

	delete [] size;
	delete [] origen;
	delete [] spacing;

	delete [] bounding_box;

	// delete [] im; // std::cout << "IM !" << std::endl;

}

void IntensityBased::SetTargetImage( Image2D::Pointer image ){ 
	m_targetImage = image;
	
	Image2D::SizeType sa = image->GetLargestPossibleRegion().GetSize();
	m_size2d = new int[3];	m_size2d[0] = sa[0];	m_size2d[1] = sa[1];	m_size2d[2] = 1;
	mm_size.clear(); mm_size.push_back(sa[0]); mm_size.push_back(sa[1]);  mm_size.push_back(sa[2]);
	
	Image2D::PointType or = image->GetOrigin();
	m_origen2d = new float[3];	m_origen2d[0] = or[0];	m_origen2d[1] = or[1];	m_origen2d[2] = 0.0;
	mm_origen.clear(); mm_origen.push_back(or[0]); mm_origen.push_back(or[1]);  mm_origen.push_back(or[2]);

	Image2D::SpacingType sp = image->GetSpacing();
	m_spacing2d = new float[3];	m_spacing2d[0] = sp[0];	m_spacing2d[1] = sp[1];	m_spacing2d[2] = 1;
	mm_spacing.clear(); mm_spacing.push_back( sp[0] ); mm_spacing.push_back( sp[1] ); mm_spacing.push_back( sp[2] );
	
	m_imageMammoPointer = image->GetBufferPointer();
}

void IntensityBased::SetSourceImage( Image3D::Pointer image ) { 
	m_sourceImage = image; 
	Image3D::SizeType sa = image->GetLargestPossibleRegion().GetSize();
	m_size3d = new int[3];
		m_size3d[0] = sa[0];
		m_size3d[1] = sa[1];
		m_size3d[2] = sa[2];
	Image3D::PointType or = image->GetOrigin();
	m_origen3d = new float[3];
		m_origen3d[0] = or[0];
		m_origen3d[1] = or[1];
		m_origen3d[2] = or[2];
	Image3D::SpacingType sp = image->GetSpacing();
	m_spacing3d = new float[3];
		m_spacing3d[0] = sp[0];
		m_spacing3d[1] = sp[1];
		m_spacing3d[2] = sp[2];
	m_imageMRIPointer = image->GetBufferPointer();
}

/* ============= Registration: Similarity Metrics ======== */
void IntensityBased::SetMetric( int metric ){ 
	m_switch = metric; 
	switch ( m_switch )	{
		case 0: { // Mutual Information:
			std::cout << "Mutual infomation Metric" << std::endl;
			break;
			}
		case 1: { // Normalized Cross-Correlation:
			std::cout << "Normalized Cross Correlation Metric" << std::endl;
			break;
			}
		case 2: { // Sum of Square Differences:
			std::cout << "Sum of Square Differences Metric " << std::endl;
			break;
			}
		case 3: { // Dice Overlap Index:
			std::cout << "Dice Overlap Index" << std::endl;
			break;
			}
		case 4: { //Gradient Correlation metric:
			std::cout << "Gradient Correlation metric: " << std::endl;
			break;
			}
		case 5: {
			std::cout << "Pattern Intensity metric" << std::endl;
			break;
			}
		case 6: {
			std::cout << "Gradient Difference metric " << std::endl;
			break;
			}
		case 7: {
			std::cout << "Local Correlation metric " << std::endl;
			break;
			}
	}
}

/* =============== Bounding box =========================== */
void IntensityBased::boundingBox(float* points, int numberOfPoints, float* bounding_box)
{
	x_max = -9999; x_min = 9999;
	y_max = -9999; y_min = 9999;
	z_max = -9999; z_min = 9999;

	for(int i=0; i<numberOfPoints; i++)
	{
		if(points[3*i] < x_min) x_min= points[3*i];
		if(points[3*i] > x_max) x_max= points[3*i];

		if(points[3*i+1] < y_min) y_min= points[3*i+1];
		if(points[3*i+1] > y_max) y_max= points[3*i+1];

		if(points[3*i+2] < z_min) z_min= points[3*i+2];
		if(points[3*i+2] > z_max) z_max= points[3*i+2];
	}

	bounding_box[0]=x_min;
	bounding_box[1]=x_max;
	bounding_box[2]=y_min;
	bounding_box[3]=y_max;
	bounding_box[4]=z_min;
	bounding_box[5]=z_max;

	std::cout << "BoundingBox: ["<< bounding_box[0]<< ", "<< bounding_box[1] << ", "<< bounding_box[2] <<
		 ", "<< bounding_box[3] << ", "<< bounding_box[4] << ", "<< bounding_box[5] << std::endl;
}

float IntensityBased::Update()
{
// ======================== Transformación del modelo !! =====================
	m_temp_parameters.mm_origen = mm_origen;
	m_temp_parameters.mm_spacing = mm_spacing;
	m_temp_parameters.mm_size = mm_size;

	Transformations * transformation = new Transformations();
		transformation->SetModelParameters( m_temp_parameters );
		transformation->SetTransformation( m_transformation );
		transformation->SetOutputDir( m_outputDir );
		
	int sys = transformation->Update();
		if(sys == 1) return 0.0;

	// == points & elements !!
	i_points = transformation->getInitialPoints();
	f_points = transformation->getFinalPoints();

	elem = m_parameters.elements;

//  ===================================  Creación de la grid regular !! ============================
	float * i_a_points = new float[ i_points.size() ];
	for(unsigned int i=0; i<i_points.size(); i++) i_a_points[i] = i_points[i];

	float * f_a_points = new float[ f_points.size() ];
	for(unsigned int i=0; i<f_points.size(); i++) f_a_points[i] = f_points[i];

	float * bb = new float[6];
	boundingBox( f_a_points, (int)(f_points.size()/3), bb);

	int * elements = new int[ elem.size() ];
	for(unsigned int i=0; i<elem.size(); i++) elements[i] = elem[i];

	float spa[3] = {5,5,5};
	RegularGrid * rr  = new RegularGrid();
		rr->setSpacing( spa );
		rr->setPoints( (int)(f_points.size()/3), f_a_points );
		rr->setElements( (int)(elem.size()/4), elements );
		rr->Update();

	// Localización del souce-point!

	if (m_temp_parameters.side=="R ") {
		m_source[0] = ( m_origen2d[0] + (m_size2d[0] * m_spacing2d[0]) ); 
		m_source[1] = ( m_origen2d[1] + (m_size2d[1] * m_spacing2d[1])/2 );
		m_source[2] = ( bb[5]- 660 );  // cambiado, jodido
	}else if (m_temp_parameters.side=="L "){
		m_source[0] = ( m_origen2d[0] ); 
		m_source[1] = ( m_origen2d[1] + (m_size2d[1] * m_spacing2d[1])/2 );
		
		m_source[2] = ( bb[5]- 660 ); // cambiado,jodido
	}

	// Registration Parameters !
	// registrationParameters * registration = new registrationParameters();

		registration->mri_origen = m_origen3d; 
		registration->mri_spacing = m_spacing3d; 
		registration->mri_size = m_size3d; 
		registration->mri_imagePointer = m_imageMRIPointer; 

		registration->mamo_origen = m_origen2d; 
		registration->mamo_spacing = m_spacing2d; 
		registration->mamo_size = m_size2d; 
		registration->mamo_imagePointer = m_imageMammoPointer; 

		// registration->simulada_imagePointer ; // Se obtiene después de la proyección !!

		registration->numberOfPoints = (int)(i_points.size()/3);
		registration->initial_points = i_a_points;
		registration->final_points = f_a_points;
		
		registration->numberOfElements = (int)(elem.size()/4);
		registration->elements = elements;
		
		registration->grid_origen = rr->getOrigen();
		registration->grid_spacing = rr->getSpacing();
		registration->grid_size = rr->getSize();

		registration->flags = rr->getFlags();
		registration->cumsum = rr->getCumSum();
		registration->correspondingElements = rr->getElementList();

		registration->source = m_source; 

// Proyección !
#ifdef _GPU_
	CudaProjection * Cp = new CudaProjection();
		Cp->SetParameters( registration );
		Cp->Update();
#endif

/* Imagen ITK */
//	Image2D::IndexType start;
		start[0] = 0; start[1] = 0; //start[2] = 0;
//	Image2D::SizeType size_im;
		size_im[0] = m_size2d[0]; size_im[1] = m_size2d[1]; //size_im[2] = size[2];
	Image2D::RegionType region(start, size_im);
//	Image2D::PointType origen_im;
		origen_im[0] = m_origen2d[0]; origen_im[1] = m_origen2d[1]; //origen_im[2] = m_origen2d[2];
//	Image2D::SpacingType spacing_im;
		spacing_im[0] = m_spacing2d[0]; spacing_im[1] = m_spacing2d[1];// spacing_im[2] = spacing[2];
//	Image2D::DirectionType direction2d;
		direction2d[0][0] =1;  direction2d[0][1] = 0; //  direction[0][2] =0;
		direction2d[1][0] =0;  direction2d[1][1] = 1; //  direction[1][2] =0;
		// direction[2][0] =0;  direction[2][1] = 0;   direction[2][2] =1;
	int numOfPix = size_im[0] * size_im[1];

//	ImportFilterType::Pointer importfilter = ImportFilterType::New();
		importfilter->SetOrigin( origen_im );
		importfilter->SetRegion( region );
		importfilter->SetSpacing( spacing_im );
		importfilter->SetDirection( direction2d );
	const bool importImagefilterWillOwnTheBuffer = false;
		//importfilter->SetImportPointer( flags, numberOfPixels, importImagefilterWillOwnTheBuffer );
		//importfilter->SetImportPointer( final_image, numberOfPixels, importImagefilterWillOwnTheBuffer );
	//float* im = CDProjection->Get2DImagePointer();
	// im = CDProjection->Get2DImagePointer();
	
#ifdef _GPU_
	//float * im = Cp->Get2DImagePointer();
	unsigned short * im = Cp->Get2DImagePointer();
#else
	float* im = 0;
#endif
		// importfilter->SetImportPointer(  Cp->Get2DImagePointer(), numOfPix, importImagefilterWillOwnTheBuffer );
		importfilter->SetImportPointer( im, numOfPix, importImagefilterWillOwnTheBuffer );
		importfilter->Update();

// 	m_outputImage = projection->getImagen();
	temp_pos =  m_temp_parameters.position ;
	outname = m_outputDir + "\\" + temp_pos + "_intermedial.nrrd";

	this->m_outputImage = importfilter->GetOutput();

	//Writer2DLType::Pointer writer2del = Writer2DLType::New();
		writer2del->SetInput( importfilter->GetOutput() );
		writer2del->SetFileName( outname );
		writer2del->SetUseCompression(true);
		try{
			writer2del->Update();
			std::cout << "Intermedial image write!! " << std::endl;
		} catch( itk::ExceptionObject & excp ) {
			std::cout << "Simulated image exception !" << std::endl;
			std::cout << excp << std::endl;
		}
/**/
	

// Comparativa de imágenes según métrica:
	//float value = 0;
	switch ( m_switch )	{
		case 0: {// Mutual Information:
			//std::cout << "Mutual infomation Metric" << std::endl;
			value = met.compute_mi( m_targetImage, importfilter->GetOutput());
			break;
			}
		case 1: {
			//std::cout << "Normalized Cross Correlation Metric" << std::endl;
			value = met.compute_ncc( m_targetImage, importfilter->GetOutput() );
			break;
			}
		case 2: {
			//std::cout << "Sum of Square Differences Metric " << std::endl;
			value = met.compute_msq( m_targetImage, importfilter->GetOutput() );
			break;
			}
		case 3: {
			//std::cout << "Dice Overlap Index" << std::endl;
			value = met.compute_Dice( m_targetImage, importfilter->GetOutput(), 0.1 );
			break;
			}
		case 4: {
			//std::cout << "Gradient Correlation metric" << std::endl;
			value = met.compute_gc( m_targetImage, importfilter->GetOutput() );
			break;
			}
		case 5: {
			//std::cout << "Pattern Intensity metric" << std::endl;
			value = met.compute_pi( m_targetImage, importfilter->GetOutput() );
			break;
			}
		case 6: {
			//std::cout << "Gradient Difference metric " << std::endl;
			value = met.compute_gd( m_targetImage, importfilter->GetOutput() );
			break;
			}
		case 7: {
			//std::cout << "Local Correlation metric " << std::endl;
//			value = met.compute_localCorrelation( m_targetImage, importfilter->GetOutput() );
			break;
			}
	}
	
	delete [] i_a_points;  // std::cout << "I_A_POINTS !" << std::endl;
	delete [] f_a_points; // std::cout << "F_A_POINTS !" << std::endl;
	delete [] bb; // std::cout << "BB !" << std::endl;
	delete [] elements; // std::cout << "ELEMENTS !" << std::endl;
	
	delete [] im; // std::cout << "IM !" << std::endl;

	delete rr; // std::cout << "RR !" << std::endl;
#ifdef _GPU_
	delete Cp;  // std::cout << "CP !" << std::endl;
#endif
	delete transformation;

	return value;	
}

#endif