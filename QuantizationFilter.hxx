#ifndef __itkQuantizationFilter_hxx
#define __itkQuantizationFilter_hxx

#include "QuantizationFilter.h"

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <math.h>
#include <random>
namespace itk
{

template<class TImage>
QuantizationFilter<TImage>::QuantizationFilter()
	{
		this->m_NumberOfBits = 32;
		//this->m_Disable = false;
		//this->SetPitch(5);
	}

template<class TImage>
void QuantizationFilter<TImage>::BeforeThreadedGenerateData()
	{
		typedef itk::MinimumMaximumImageCalculator< TImage >  CalculatorType;
		typename CalculatorType::Pointer calculatorI = CalculatorType::New();
		calculatorI->SetImage(this->GetInput());
		calculatorI->Compute();
		m_Maximum = calculatorI->GetMaximum();
		m_Minimum = calculatorI->GetMinimum();
		typename TImage::ValueType DataRange = m_Maximum - m_Minimum; //Pixel type , typically float
		m_Pitch = DataRange / pow(2.0, m_NumberOfBits); // double

		// coutS
		//std::cout << "m_Pitch = " << m_Pitch << std::endl;
	}


template< class TImage>
void QuantizationFilter< TImage>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
{

	typename TImage::ConstPointer input = this->GetInput();
	typename TImage::Pointer output = this->GetOutput();

	itk::ImageRegionIterator<TImage> out(output, outputRegionForThread);
	itk::ImageRegionConstIterator<TImage> in(input, outputRegionForThread);

	
	long long int BinCounter = -1;
	//typename TImage::ValueType RelativePosition = -1.0; // (in.GetValue() - m_Minimum) / m_Pitch;
	double RelativePosition = 0.0;
	double NewPixelValue; // higher precision
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);


	while (!out.IsAtEnd())
	{
	
		RelativePosition = (in.Get() - m_Minimum) / m_Pitch;
		if (RelativePosition == std::floor(RelativePosition)) // use == for double, check if the difference is smaller than the min pos double
		{
			out.Set(in.Get());
			++in;
			++out;
			continue;
		}
		BinCounter = static_cast<long long int>(std::floor(RelativePosition));
		
		if (in.Get() == m_Maximum)
		{
			out.Set(m_Maximum);
			++in;
			++out;
			continue;
		}
		/*if (in.Get() == m_Minimum)
		{
			out.Set(m_Minimum);
			++in;
			++out;
			continue;
		}
		*/
		NewPixelValue = m_Minimum+BinCounter*m_Pitch;
		//out.Set(in.Get()); 
		/**/

		double MaxAllowed = pow(2.0, -1 * m_NumberOfBits);
		if (dist(mt) > (RelativePosition - BinCounter))
		{
			/*double diff = static_cast<typename TImage::ValueType>(NewPixelValue) - in.Get();
			if ( std::abs(diff) > MaxAllowed)
			{
				std::cout << "Error at Qfilter." << std::endl;
				std::cout << "MaxAllowed" << MaxAllowed << std::endl;
				std::cout << "diff." << diff << std::endl;
			}*/
			//if diff
			out.Set(static_cast<typename TImage::ValueType>(NewPixelValue));
		}
		else
			out.Set(static_cast<typename TImage::ValueType>(NewPixelValue + m_Pitch));
			
/*
		if (std::abs(in.Get() - out.Get()) > pow(2.0, -1 * m_NumberOfBits))
		{
			
			//std::cout<<
			std::cout << "Diff = "<<std::abs(in.Get() - out.Get()) << std::endl;
			std::cout << "Max allowed: " << pow(2.0, -1 * m_NumberOfBits) << std::endl;

			std::cout << "Error at Qfilter." << std::endl;
		}*/
		++in;
		++out;
	}
	//cout
	//std::cout << "BinCounter = " << BinCounter << std::endl;
	//DD(BinCounter)




}

}// end namespace


#endif
