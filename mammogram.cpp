
#ifndef __mammogram_cpp
#define __mammogram_cpp

#include "mammogram.h"

mammogram::mammogram(void)
{
	this->m_image = ImageMammoType::New();
	this->has_mammogram = false;
	this->m_imageIO = IOType::New();
	
	this->has_segmentation = false;
	this->m_segmentation = ImageMammoType::New();

	this->has_densityMap = false;

	this->has_im_dense_tissue = false;

	m_pectoralMuscleAngle = 0.0;
}

mammogram::~mammogram(void)
{
}

// =========================================================

void mammogram::read_mammogram(std::string inputfilename)
{
	this->m_filename = inputfilename;
	ReaderMammoType::Pointer reader_mammo= ReaderMammoType::New();
		reader_mammo->SetFileName( inputfilename);
		try{
			reader_mammo->Update();
		}catch(itk::ExceptionObject & err){
			std::cout <<" " << std::endl;
			std::cout <<"Image File Reader Fails! " << std::endl;
			std::cout << err << std::endl;
			//return EXIT_FAILURE;
		}

	this->m_image = reader_mammo->GetOutput();
	this->has_mammogram  = true;
}

void mammogram::read_mammogram_with_metadata(std::string inputfilename)
{
	this->m_filename = inputfilename;
	ReaderMammoType::Pointer reader_mammo= ReaderMammoType::New();
		reader_mammo->SetFileName( inputfilename);
		reader_mammo->SetImageIO( this->m_imageIO);
		try{
			reader_mammo->Update();
		}catch(itk::ExceptionObject & err){
			std::cout <<" " << std::endl;
			std::cout <<"Image File Reader Fails! " << std::endl;
			std::cout << err << std::endl;
			//return EXIT_FAILURE;
		}

	this->m_image = reader_mammo->GetOutput();
	this->metadata = this->m_imageIO->GetMetaDataDictionary();
	this->has_mammogram = true;
	this->has_dictionary = true;
}

void mammogram::read_densityMap(std::string inputfilename)
{
	this->m_filename = inputfilename;
	ReaderMammoType::Pointer reader_dm= ReaderMammoType::New();
		reader_dm->SetFileName( inputfilename);
		try{
			reader_dm->Update();
		}catch(itk::ExceptionObject & err){
			std::cout <<" " << std::endl;
			std::cout <<"Image File Reader Fails! " << std::endl;
			std::cout << err << std::endl;
			//return EXIT_FAILURE;
		}

	if ( getImageLaterality(false) == "R ") 
	{
		std::cout << "\t Fliping Density Map..... " << std::endl;
		itk::FixedArray<bool, 2> flipAxes;
			flipAxes[0] = true;
			flipAxes[1] = false;
		FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
			flipFilter->SetInput( reader_dm->GetOutput() );
			flipFilter->SetFlipAxes(flipAxes);
			flipFilter->Update();

		this->m_densityMap = flipFilter->GetOutput();
		this->has_densityMap  = true;
	}else{
		this->m_densityMap = reader_dm->GetOutput();
		this->has_densityMap  = true;
	}

	if(this->has_mammogram)
	{
		this->m_densityMap->SetOrigin( this->m_image->GetOrigin());

		ImageMammoType::SpacingType dm_spa;
		ImageMammoType::SizeType dm_size = this->m_densityMap->GetLargestPossibleRegion().GetSize();
		
		ImageMammoType::SpacingType spacing = this->m_image->GetSpacing();
		ImageMammoType::SizeType size = this->m_image->GetLargestPossibleRegion().GetSize();

			dm_spa[0] = (spacing[0]*size[0])/dm_size[0];
			dm_spa[1] = (spacing[1]*size[1])/dm_size[1];
		this->m_densityMap->SetSpacing( dm_spa );
	}

	computeDenseTissue(false);
}

void mammogram::read_pectoralmuscleMask(std::string inputfilename)
{
	ReaderMammoType::Pointer reader_mask= ReaderMammoType::New();
		reader_mask->SetFileName( inputfilename);
		try{
			reader_mask->Update();
		}catch(itk::ExceptionObject & err){
			std::cout <<" " << std::endl;
			std::cout <<"Image File Reader Fails! " << std::endl;
			std::cout << err << std::endl;
			//return EXIT_FAILURE;
		}

	this->m_pectoralMuscleMask = reader_mask->GetOutput();
	
	this->m_pectoralMuscleMask->SetOrigin( this->m_densityMap->GetOrigin() );
	this->m_pectoralMuscleMask->SetSpacing( this->m_densityMap->GetSpacing() );
	this->m_pectoralMuscleMask->SetDirection( this->m_densityMap->GetDirection() );

	//Angulo del músculo pectoral!
	computeAnglePectoralMuscle(this->m_pectoralMuscleMask);

	if ( getImageLaterality(true) == "R ") 
	{
		std::cout << "\t Fliping Density Map..... " << std::endl;
		itk::FixedArray<bool, 2> flipAxes;
			flipAxes[0] = true;
			flipAxes[1] = false;
		FlipImageFilterType::Pointer flipmaskFilter = FlipImageFilterType::New();
			flipmaskFilter->SetInput( this->m_pectoralMuscleMask );
			flipmaskFilter->SetFlipAxes(flipAxes);
			flipmaskFilter->Update();

		this->m_pectoralMuscleMask = flipmaskFilter->GetOutput();
	} //else{
	//	this->m_pectoralMuscleMask = this->m_pectoralMuscleMask;
	// }

	this->m_pectoralMuscleMask->SetOrigin( this->m_densityMap->GetOrigin() );
	this->m_pectoralMuscleMask->SetSpacing( this->m_densityMap->GetSpacing() );
	this->m_pectoralMuscleMask->SetDirection( this->m_densityMap->GetDirection() );

	if (has_densityMap){
		//TransformType::Pointer transform2 = TransformType::New();
		//	transform2->SetIdentity();
		//NNType::Pointer nn2 = NNType::New();
		//ResampleFilterType::Pointer resample2 = ResampleFilterType::New();
		//	resample2->SetInput( this->m_pectoralMuscleMask );
		//	
		//	resample2->SetSize( this->m_densityMap->GetLargestPossibleRegion().GetSize() );
		//	resample2->SetOutputSpacing( this->m_densityMap->GetSpacing() );
		//	resample2->SetOutputOrigin( this->m_densityMap->GetOrigin() );
		//	resample2->SetOutputDirection( this->m_densityMap->GetDirection() );

		//	resample2->SetInterpolator( nn2 );

		//	resample2->SetTransform( transform2 );
		//	resample2->Update();

		MaskFilterType::Pointer maskFilter2 = MaskFilterType::New();
			maskFilter2->SetInput( this->m_densityMap );
			maskFilter2->SetMaskingValue( 0 );
			maskFilter2->SetMaskImage( this->m_pectoralMuscleMask );
			maskFilter2->Update();

		this->m_densityMap = maskFilter2->GetOutput();
	}

	if (has_im_dense_tissue){
		//TransformType::Pointer transform3 = TransformType::New();
		//	transform3->SetIdentity();
		//NNType::Pointer nn3 = NNType::New();
		//ResampleFilterType::Pointer resample3 = ResampleFilterType::New();
		//	resample3->SetInput( m_pectoralMuscleMask );
		//	
		//	resample3->SetSize( this->m_densityMap->GetLargestPossibleRegion().GetSize() );
		//	resample3->SetOutputSpacing( this->m_densityMap->GetSpacing() );
		//	resample3->SetOutputOrigin( this->m_densityMap->GetOrigin() );
		//	resample3->SetOutputDirection( this->m_densityMap->GetDirection() );
		//	
		//	resample3->SetInterpolator( nn3 );

		//	resample3->SetTransform( transform3 );
		//	resample3->Update();

		MaskFilterType::Pointer maskFilter3 = MaskFilterType::New();
			maskFilter3->SetInput( this->im_dense_tissue );
			maskFilter3->SetMaskingValue( 0 );
			maskFilter3->SetMaskImage( this->m_pectoralMuscleMask );
			maskFilter3->Update();

		this->im_dense_tissue = maskFilter3->GetOutput();
	}

	if (has_mammogram){
		TransformType::Pointer transform = TransformType::New();
			transform->SetIdentity();
		NNType::Pointer nn = NNType::New();
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
			resample->SetInput( this->m_pectoralMuscleMask );

			resample->SetSize( this->m_image->GetLargestPossibleRegion().GetSize() );
			resample->SetOutputSpacing( this->m_image->GetSpacing() );
			resample->SetOutputOrigin( this->m_image->GetOrigin() );
			resample->SetOutputDirection( this->m_image->GetDirection() );

			resample->SetInterpolator( nn );
			
			resample->SetTransform( transform );
			resample->Update();

		MaskFilterType::Pointer maskFilter = MaskFilterType::New();
			maskFilter->SetInput( this->m_image );
			maskFilter->SetMaskingValue( 0 );
			maskFilter->SetMaskImage( resample->GetOutput() );
			maskFilter->Update();

		this->m_image = maskFilter->GetOutput();
	}


}

// ==========================================================
void mammogram::write_mammogram( std::string outputfilename)
{
	WriterMammoType::Pointer writer = WriterMammoType::New();
		writer->SetInput( this->m_image );
		//writer->SetImageIO(this->m_imageIO );
		writer->SetMetaDataDictionary( this->m_imageIO->GetMetaDataDictionary() );
		writer->SetFileName(outputfilename);
	try{
		writer->Update();
	}catch(itk::ExceptionObject & excp ){
		std::cout << "Mammogram writer exception! " << std::endl;
		std::cout << excp << std::endl;
	}
}

void mammogram::write_densitymap( std::string outputfilename)
{
	WriterMammoType::Pointer writer = WriterMammoType::New();
		writer->SetInput( get_dense_tissue(false) );
		//writer->SetImageIO( this->m_imageIO );
		//writer->SetMetaDataDictionary( this->m_imageIO->GetMetaDataDictionary() );
		writer->SetFileName(outputfilename);
	try{
		writer->Update();
	}catch(itk::ExceptionObject & excp ){
		std::cout << "Density Map writer exception! " << std::endl;
		std::cout << excp << std::endl;
	}
}

// ==========================================================
ImageMammoType::Pointer mammogram::get_image()
{
	if(this->has_mammogram){
		return this->m_image;
	} else {
		std::cout << "class mammogram is empty" << std::endl;
		return 0;
	}
}

ImageMammoType::Pointer mammogram::get_DensityMap()
{
	if(this->has_densityMap) {
		return this->m_densityMap;
	} else {
		std::cout << "Class mammogram has no density map! " << std::endl;
		return 0;
	}
}

void mammogram::computeDenseTissue(bool verbose)
{
	if(verbose)	std::cout << "Computing Dense Tissue from Volpara(R) Density Map......." << std::endl;

	im_dense_tissue = ImageMammoType::New();
	im_dense_tissue = this->m_densityMap;
	
	if(this->has_densityMap & this->has_mammogram) {
		im_dense_tissue->SetOrigin( this->m_densityMap->GetOrigin() );
		im_dense_tissue->SetSpacing( this->m_densityMap->GetSpacing() );

		int H = getBodyPartThickness(false);
		// DENSE_TISSUE_HEIGHT_MM = (PIXEL_VALUE / 65535 -0.5) *2.0* H;
		itk::ImageRegionIterator<ImageMammoType> it= itk::ImageRegionIterator<ImageMammoType>(im_dense_tissue, im_dense_tissue->GetLargestPossibleRegion() );
		for(it.Begin(); !it.IsAtEnd(); it++)
		{
			//if(it.Get()!=32767)	std::cout << it.Get() << std::endl;
			float a =( (((float)it.Get()/65535) - 0.5) * 2.0 * H )*1000; // um of dense tissue.
			//if(it.Get()!=32767) std::cout << a << std::endl;
			if( a>0 ) {
				it.Set( ceil(a));
			} else {
				it.Set( 0 );
			}
		}

	} else {
		std::cout << " " << std::endl;
		std::cout << "There is no density map or mammogram avaliable. " << std::endl;
		std::cout << "Check both images & mammogram's dicomtags. " << std::endl;
	}

	this->has_im_dense_tissue = true;
}

ImageMammoType::Pointer mammogram::get_dense_tissue(bool verbose)
{ 
	if(!has_im_dense_tissue) computeDenseTissue( true );
	
	return this->im_dense_tissue;
}

void mammogram::computeAnglePectoralMuscle(ImageMammoType::Pointer image)
{
	IteratorType it( image, image->GetLargestPossibleRegion() );
	
	it.SetDirection(0); // X	
	it.GoToBegin();
	bool cont = true;
	itk::Index< 2 > x_b;

	while( cont && !it.IsAtEndOfLine() ){
		if( it.Get() != 0 ) {
			x_b=it.GetIndex();
			cont=false;
		}
		it++;
	}

	//std::cout << "Index En X " << x_b << std::endl;

	it.SetDirection(1); // Y
	it.GoToBegin();
	cont = true;
	itk::Index< 2 > y_a;
	
	while( cont && !it.IsAtEnd() ) {
		while( cont && !it.IsAtEndOfLine() )
		{
			if( it.Get() != 0 ) {
				y_a=it.GetIndex();
				cont=false;		
			}
			it++;
		}
		it.NextLine();
	}

	//std::cout << "Index En Y " << y_a << std::endl;
	
	double m[2] = {0.0,0.0};
		m[0] = x_b[0] - y_a[0];
		m[1] = x_b[1] - y_a[1];

	double m0[2] = {0.0,1.0};

	//std::cout << "Pectoral Muscle Direction: [" << m[0] << ", " << m[1] << "] " << std::endl;
	//std::cout << "Direction a tener en cuenta: [" << m0[0] << ", " << m0[1] << "] " << std::endl;

	//double angle = atan( abs(double(x_b[0]-y_a[0])/double(x_b[1]-y_a[1]) )) * 180 /M_PI ;
	this->m_pectoralMuscleAngle = atan( abs(double(x_b[0]-y_a[0])/double(x_b[1]-y_a[1]) )) +.2; 
	std::cout << "Angle : " <<atan(abs( double(x_b[0]-y_a[0])/double(x_b[1]-y_a[1]) )) * 180 / (double)M_PI<< std::endl;

}

// ===========================================================
void mammogram::otsu_segmentation()
{
	OtsuMammoType::Pointer otsu = OtsuMammoType::New();
		otsu->SetInput( this->m_image );
		//otsu->set
		otsu->SetOutsideValue(1);
		otsu->SetInsideValue(0);
		try {
			otsu->Update();
		} catch( itk::ExceptionObject & err) {
			std::cout << " " << std::endl;
			std::cout << "Otsu threshold Fails! " << std::endl;
			std::cout << err << std::endl;
			// return EXIT_FAILURE;
		}
	this->has_segmentation = true;
	this->m_segmentation = otsu->GetOutput();
}

ImageMammoType::Pointer mammogram::get_otsu_segmentation()
{
	if( this->has_segmentation ) {
		return this->m_segmentation;
	} else {
		try {
			otsu_segmentation();
			return this->m_segmentation;
		} catch ( int e )  {
			std::cout << "Unable Otsu Segmentation !! " << std::endl;
			std::cout << e << std::endl;
			return 0;
		}
	}

}

ImageMammoType::Pointer mammogram::get_otsu_densityMap()
{
	if( this->has_segmentation ) {
		TransformType::Pointer transform = TransformType::New();
			transform->SetIdentity();
		NNType::Pointer nn = NNType::New();
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
			resample->SetInput( this->m_segmentation );
			resample->SetSize( this->m_densityMap->GetLargestPossibleRegion().GetSize() );
			resample->SetOutputSpacing( this->m_densityMap->GetSpacing() );
			resample->SetInterpolator( nn );
			resample->SetTransform( transform );
			resample->Update();

		return resample->GetOutput();

	} else {
		try {
			otsu_segmentation();

			TransformType::Pointer transform = TransformType::New();
				transform->SetIdentity();
			NNType::Pointer nn = NNType::New();
			ResampleFilterType::Pointer resample = ResampleFilterType::New();
				resample->SetInput( this->m_segmentation );
				resample->SetSize( this->m_densityMap->GetLargestPossibleRegion().GetSize() );
				resample->SetOutputSpacing( this->m_densityMap->GetSpacing() );
				resample->SetInterpolator( nn );
				resample->SetTransform( transform );
				resample->Update();

			return resample->GetOutput();
		} catch ( int e )  {
			std::cout << "Unable Otsu Segmentation !! " << std::endl;
			std::cout << e << std::endl;
			return 0;
		}
	}

}

void mammogram::getContour()
{
	if(this->has_segmentation)
	{
		ContourMammoType::Pointer contour = ContourMammoType::New();
			contour->SetInput( this->m_segmentation );
			contour->SetBackgroundValue(0);
			contour->SetForegroundValue(1);
			contour->Update();
	}else if(this->has_mammogram){
		otsu_segmentation();
		ContourMammoType::Pointer contour = ContourMammoType::New();
			contour->SetInput( this->m_segmentation );
			contour->SetBackgroundValue(0);
			contour->SetForegroundValue(1);
			contour->Update();
	}else{ std::cout << "\t Mammogram has no segmentation" << std::endl; }
}

std::vector<float> mammogram::getCentroid(bool verbose)
{
	std::vector<float> centroid;
	I2LType::Pointer i2l = I2LType::New();
		i2l->SetInput( this->m_segmentation );
		i2l->SetComputePerimeter(true);
		i2l->Update();
		
	LabelMapType *labelMap = i2l->GetOutput();
	ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(0);
	
	if(verbose)	std::cout << "\t Centroid: "<< labelObject->GetCentroid() << std::endl;

	centroid.push_back(labelObject->GetCentroid()[0]);
	centroid.push_back(labelObject->GetCentroid()[1]);
	centroid.push_back( 0 );

	return centroid;
}

int mammogram::getDistanceSourceToDetector(bool verbose)
{
	if( this->has_dictionary )
	{
		std::string tagkey = "0018|1110";
		std::string labelId;

		bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
		if (found && verbose) std::cout << "\t Distance Source To Detector: " << labelId <<std::endl;

		return atoi( labelId.c_str() );
	} else { return 660;}
}

float mammogram::getFocalSpot(bool verbose)
{
	std::string tagkey = "0018|1190";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t Focal Spot: " << labelId <<std::endl;

	return atof(labelId.c_str());
}

int mammogram::getBodyPartThickness(bool verbose)
{	
	std::string tagkey = "0018|11a0";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t Body Part Thickness: " << labelId <<std::endl;

	return atof(labelId.c_str() );
}

int mammogram::getCompressionForce(bool verbose)
{	
	std::string tagkey = "0018|1190";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t Compression Force: " << labelId <<std::endl;

	return atoi(labelId.c_str());
}

std::string mammogram::getViewPossition(bool verbose)
{	
	std::string tagkey = "0018|5101";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t View Possiton: " << labelId <<std::endl;

	return labelId;
}

int mammogram::getPositionerPrimaryAngle(bool verbose)
{	
	std::string tagkey = "0018|1510";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t Positioner Primary Angle: " << labelId <<std::endl;

	return atof(labelId.c_str());
}

std::string mammogram::getAnodeTargetMaterial(bool verbose)
{	
	std::string tagkey = "0018|1191";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t Anode Target Material: " << labelId <<std::endl;

	return labelId;
}

std::string mammogram::getFilterMaterial(bool verbose)
{	
	std::string tagkey = "0018|7050";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t Filter Material: " << labelId <<std::endl;

	return labelId;
}

int mammogram::getKVP(bool verbose)
{	
	std::string tagkey = "0018|0060";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found && verbose) std::cout << "\t KVP: " << labelId <<std::endl;

	return atoi(labelId.c_str());
}

std::string mammogram::getImageLaterality(bool verbose)
{	
	std::string tagkey = "0020|0062";
	std::string labelId;

	bool found = this->m_imageIO->GetValueFromTag( tagkey, labelId );
	if (found){ if(verbose){ std::cout << "\t Image Laterality: " << labelId <<std::endl;}}
	else{ std::cout << "\t No Image Laterality. " << std::endl; }

	return labelId;
}

#endif