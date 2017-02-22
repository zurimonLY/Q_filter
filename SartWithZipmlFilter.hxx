/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef rtkSartWithZipmlFilter_hxx
#define rtkSartWithZipmlFilter_hxx

#include "SartWithZipmlFilter.h"


#include <algorithm>
#include <itkTimeProbe.h>

namespace rtk
{
	template<class TInputImage, class TOutputImage>
	SartWithZipmlFilter<TInputImage, TOutputImage>
		::SartWithZipmlFilter()
	{
		this->SetNumberOfRequiredInputs(2);

		// Set default parameters
		m_EnforcePositivity = false;
		m_IsGated = false;
		m_Quantization = true;//temp
		m_MyVerbose = false;
		m_NumberOfIterations = 3;
		m_Lambda = 0.3;

		// Create each filter of the composite filter
		m_ExtractFilter = ExtractFilterType::New();
		m_ZeroMultiplyFilter = MultiplyFilterType::New();
		m_SubtractFilter = SubtractFilterType::New();
		m_AddFilter = AddFilterType::New();
		m_DisplacedDetectorFilter = DisplacedDetectorFilterType::New();
		m_MultiplyFilter = MultiplyFilterType::New();
		m_GatingWeightsFilter = GatingWeightsFilterType::New();
		m_ConstantVolumeSource = ConstantImageSourceType::New();

		// Create the filters required for correct weighting of the difference
		// projection
		m_ExtractFilterRayBox = ExtractFilterType::New();
		m_RayBoxFilter = RayBoxIntersectionFilterType::New();
		m_DivideFilter = DivideFilterType::New();
		m_ConstantProjectionStackSource = ConstantImageSourceType::New();

		// Create the filter that enforces positivity
		m_ThresholdFilter = ThresholdFilterType::New();

		m_MyWriter01 = ImageFileWriterType::New();
		m_MyWriter02 = ImageFileWriterType::New();
		//m_MyWriter->SetFileName("MyProbe.mhd");
		// Quantization filters
		m_QuantRBIFilter = QuantizationFilterProjectionType::New();				//[0]
		m_QuantProjExtractionFilter = QuantizationFilterProjectionType::New();	//[1]
		m_QuantFPInputFilter = QuantizationFilterVolumeType::New();				//[2]
		m_QuantFPResultFilter = QuantizationFilterProjectionType::New();		//[3]
		m_QuantDiffResultFilter = QuantizationFilterProjectionType::New();		//[4]
		m_QuantB4NormFilter = QuantizationFilterProjectionType::New();			//[5]
		m_QuantBPInput = QuantizationFilterProjectionType::New();				//[6]
		m_QuantAdderInput = QuantizationFilterVolumeType::New();				//[7]

		m_Subtractor = SubtractImageFilterType::New();
		m_MinMaxCalculator = MinMaxCalcType::New();
		/*
		m_QuantRBIFilter
		m_QuantProjExtractionFilter
		m_QuantFPInputFilter
		m_QuantFPResultFilter
		m_QuantDiffResultFilter
		m_QuantB4NormFilter
		m_QuantBPInput
		*/

		//Permanent internal connections
		m_ZeroMultiplyFilter->SetInput1(itk::NumericTraits<typename InputImageType::PixelType>::ZeroValue());
		m_ZeroMultiplyFilter->SetInput2(m_ExtractFilter->GetOutput());


		m_ExtractFilterRayBox->SetInput(m_ConstantProjectionStackSource->GetOutput());
		m_RayBoxFilter->SetInput(m_ExtractFilterRayBox->GetOutput());


		m_QuantProjExtractionFilter->SetInput(m_ExtractFilter->GetOutput());

		m_SubtractFilter->SetInput(0, m_QuantProjExtractionFilter->GetOutput());




		m_MultiplyFilter->SetInput1(itk::NumericTraits<typename InputImageType::PixelType>::ZeroValue());// will be replaced by lambda

		m_QuantDiffResultFilter->SetInput(m_SubtractFilter->GetOutput());

		m_MultiplyFilter->SetInput2(m_QuantDiffResultFilter->GetOutput());


		m_QuantRBIFilter->SetInput(m_RayBoxFilter->GetOutput());
		m_QuantB4NormFilter->SetInput(m_MultiplyFilter->GetOutput());


		m_DivideFilter->SetInput1(m_QuantB4NormFilter->GetOutput());
		m_DivideFilter->SetInput2(m_QuantRBIFilter->GetOutput());


		m_DisplacedDetectorFilter->SetInput(m_DivideFilter->GetOutput());



		// Default parameters
		m_ExtractFilter->SetDirectionCollapseToSubmatrix();
		m_ExtractFilterRayBox->SetDirectionCollapseToSubmatrix();
		m_IsGated = false;
		m_NumberOfProjectionsPerSubset = 1; //Default is the SART behavior
		m_DisplacedDetectorFilter->SetPadOnTruncatedSide(false);
		m_DisableDisplacedDetectorFilter = false;
	}

	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::SetForwardProjectionFilter(int _arg)
	{
		if (_arg != this->GetForwardProjectionFilter())
		{
			Superclass::SetForwardProjectionFilter(_arg);
			m_ForwardProjectionFilter = this->InstantiateForwardProjectionFilter(_arg);
		}
	}

	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::SetBackProjectionFilter(int _arg)
	{
		if (_arg != this->GetBackProjectionFilter())
		{
			Superclass::SetBackProjectionFilter(_arg);
			m_BackProjectionFilter = this->InstantiateBackProjectionFilter(_arg);
		}
	}

	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::SetGatingWeights(std::vector<float> weights)
	{
		m_GatingWeights = weights;
		m_IsGated = true;
	}

	/*
	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::SetNumberOfBits(std::vector<int> & NumberOfBits_arg)// pass a ref
	{
		m_NumberOfBits = NumberOfBits_arg;
		m_Quantization = true;
		
		m_QuantRBIFilter->SetNumberOfBits(m_NumberOfBits[0]);
		m_QuantProjExtractionFilter->SetNumberOfBits(m_NumberOfBits[1]);
		m_QuantFPInputFilter->SetNumberOfBits(m_NumberOfBits[2]);
		m_QuantFPResultFilter->SetNumberOfBits(m_NumberOfBits[3]);
		m_QuantDiffResultFilter->SetNumberOfBits(m_NumberOfBits[4]);
		m_QuantB4NormFilter->SetNumberOfBits(m_NumberOfBits[5]);
		m_QuantBPInput->SetNumberOfBits(m_NumberOfBits[6]);
		m_QuantAdderInput->SetNumberOfBits(m_NumberOfBits[7]);
		
		
	}*/


	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::GenerateInputRequestedRegion()
	{
		typename Superclass::InputImagePointer inputPtr =
			const_cast<TInputImage *>(this->GetInput());

		if (!inputPtr)
			return;

		if (m_EnforcePositivity)//please always enforce 
		{
			m_ThresholdFilter->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
			m_ThresholdFilter->GetOutput()->PropagateRequestedRegion();
		}
		else
		{
			m_QuantAdderInput->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
			m_QuantAdderInput->GetOutput()->PropagateRequestedRegion();
			//m_BackProjectionFilter->GetOutput()->SetRequestedRegion(this->GetOutput()->GetRequestedRegion());
			//m_BackProjectionFilter->GetOutput()->PropagateRequestedRegion();
		}
	}

	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::GenerateOutputInformation()
	{
		m_DisplacedDetectorFilter->SetDisable(m_DisableDisplacedDetectorFilter);

		// We only set the first sub-stack at that point, the rest will be
		// requested in the GenerateData function
		typename ExtractFilterType::InputImageRegionType projRegion;

		projRegion = this->GetInput(1)->GetLargestPossibleRegion();
		m_ExtractFilter->SetExtractionRegion(projRegion);
		m_ExtractFilterRayBox->SetExtractionRegion(projRegion);

		// Links with the forward and back projection filters should be set here
		// and not in the constructor, as these filters are set at runtime

		//m_QuantRBIFilter->SetInput(m_ExtractFilterRayBox->GetOutput());


		//
		if (m_MyVerbose)
		{
			std::cout << "Before m_ConstantVolumeSource->UpdateOutputInformation(); in GenerateOutputInformation" << std::endl; std::cin.get();
		}

		m_ConstantVolumeSource->SetInformationFromImage(const_cast<TInputImage *>(this->GetInput(0)));
		m_ConstantVolumeSource->SetConstant(0);
		m_ConstantVolumeSource->UpdateOutputInformation();

		m_QuantBPInput->SetInput(m_DisplacedDetectorFilter->GetOutput());


		m_BackProjectionFilter->SetInput(0, m_ConstantVolumeSource->GetOutput());
		m_BackProjectionFilter->SetInput(1, m_QuantBPInput->GetOutput());
		m_BackProjectionFilter->SetTranspose(false);

		m_QuantAdderInput->SetInput(m_BackProjectionFilter->GetOutput());
		m_AddFilter->SetInput1(m_QuantAdderInput->GetOutput());
		m_AddFilter->SetInput2(this->GetInput(0));

		m_QuantFPInputFilter->SetInput(this->GetInput(0));

		m_ForwardProjectionFilter->SetInput(0, m_ZeroMultiplyFilter->GetOutput());
		m_ForwardProjectionFilter->SetInput(1, m_QuantFPInputFilter->GetOutput());
		//m_ExtractFilter->SetInput(m_QuantFPInputFilter->GetOutput());
		m_ExtractFilter->SetInput(this->GetInput(1));

		m_QuantFPResultFilter->SetInput(m_ForwardProjectionFilter->GetOutput());
		m_SubtractFilter->SetInput(1, m_QuantFPResultFilter->GetOutput());

		// For the same reason, set geometry now
		// Check and set geometry
		if (this->GetGeometry().GetPointer() == ITK_NULLPTR)
		{
			itkGenericExceptionMacro(<< "The geometry of the reconstruction has not been set");
		}
		m_ForwardProjectionFilter->SetGeometry(this->m_Geometry);
		m_BackProjectionFilter->SetGeometry(this->m_Geometry.GetPointer());
		m_DisplacedDetectorFilter->SetGeometry(this->m_Geometry);

		if (m_IsGated) // For gated SART, insert a gating filter into the pipeline
		{
			m_GatingWeightsFilter->SetInput1(m_DivideFilter->GetOutput());
			m_GatingWeightsFilter->SetConstant2(1);
			m_DisplacedDetectorFilter->SetInput(m_GatingWeightsFilter->GetOutput());
		}
		if (m_MyVerbose)
		{
			std::cout << "Before m_ConstantProjectionStackSource->UpdateOutputInformation(); in GenerateOutputInformation" << std::endl; std::cin.get();
		}

		m_ConstantProjectionStackSource->SetInformationFromImage(const_cast<TInputImage *>(this->GetInput(1)));
		m_ConstantProjectionStackSource->SetConstant(0);
		m_ConstantProjectionStackSource->UpdateOutputInformation();


		// Create the m_RayBoxFiltersectionImageFilter
		m_RayBoxFilter->SetGeometry(this->GetGeometry().GetPointer());
		itk::Vector<double, 3> Corner1, Corner2;

		Corner1[0] = this->GetInput(0)->GetOrigin()[0];
		Corner1[1] = this->GetInput(0)->GetOrigin()[1];
		Corner1[2] = this->GetInput(0)->GetOrigin()[2];
		Corner2[0] = this->GetInput(0)->GetOrigin()[0] + this->GetInput(0)->GetLargestPossibleRegion().GetSize()[0] * this->GetInput(0)->GetSpacing()[0];
		Corner2[1] = this->GetInput(0)->GetOrigin()[1] + this->GetInput(0)->GetLargestPossibleRegion().GetSize()[1] * this->GetInput(0)->GetSpacing()[1];
		Corner2[2] = this->GetInput(0)->GetOrigin()[2] + this->GetInput(0)->GetLargestPossibleRegion().GetSize()[2] * this->GetInput(0)->GetSpacing()[2];

		m_RayBoxFilter->SetBoxMin(Corner1);
		m_RayBoxFilter->SetBoxMax(Corner2);

		if (m_EnforcePositivity)
		{
			m_ThresholdFilter->SetOutsideValue(0);
			m_ThresholdFilter->ThresholdBelow(0);
			m_ThresholdFilter->SetInput(m_AddFilter->GetOutput());

			// Update output information

			m_ThresholdFilter->UpdateOutputInformation();
			this->GetOutput()->SetOrigin(m_ThresholdFilter->GetOutput()->GetOrigin());
			this->GetOutput()->SetSpacing(m_ThresholdFilter->GetOutput()->GetSpacing());
			this->GetOutput()->SetDirection(m_ThresholdFilter->GetOutput()->GetDirection());
			this->GetOutput()->SetLargestPossibleRegion(m_ThresholdFilter->GetOutput()->GetLargestPossibleRegion());
		}
		else
		{
			// Update output information
			if (m_MyVerbose)
			{
				std::cout << "Before m_AddFilter->UpdateOutputInformation(); in GenerateOutputInformation" << std::endl; std::cin.get();


				DD(m_AddFilter->GetOutput()->GetLargestPossibleRegion())
			}

			m_AddFilter->UpdateOutputInformation();
			this->GetOutput()->SetOrigin(m_AddFilter->GetOutput()->GetOrigin());
			this->GetOutput()->SetSpacing(m_AddFilter->GetOutput()->GetSpacing());
			this->GetOutput()->SetDirection(m_AddFilter->GetOutput()->GetDirection());
			this->GetOutput()->SetLargestPossibleRegion(m_AddFilter->GetOutput()->GetLargestPossibleRegion());
			if (m_MyVerbose)
			{
				DD(m_AddFilter->GetOutput()->GetLargestPossibleRegion())
			}
		}

		// Set memory management flags
		m_ZeroMultiplyFilter->ReleaseDataFlagOn();
		m_ForwardProjectionFilter->ReleaseDataFlagOn();
		m_SubtractFilter->ReleaseDataFlagOn();
		m_MultiplyFilter->ReleaseDataFlagOn();
		m_RayBoxFilter->ReleaseDataFlagOn();
		m_DivideFilter->ReleaseDataFlagOn();
		m_DisplacedDetectorFilter->ReleaseDataFlagOn();

		if (m_EnforcePositivity)
			m_AddFilter->ReleaseDataFlagOn();

	}

	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::GenerateData()
	{
		m_QuantRBIFilter->SetNumberOfBits(m_NumberOfBitsQ00);
		m_QuantProjExtractionFilter->SetNumberOfBits(m_NumberOfBitsQ01);
		m_QuantFPInputFilter->SetNumberOfBits(m_NumberOfBitsQ02);
		m_QuantFPResultFilter->SetNumberOfBits(m_NumberOfBitsQ03);
		m_QuantDiffResultFilter->SetNumberOfBits(m_NumberOfBitsQ04);
		m_QuantB4NormFilter->SetNumberOfBits(m_NumberOfBitsQ05);
		m_QuantBPInput->SetNumberOfBits(m_NumberOfBitsQ06);
		m_QuantAdderInput->SetNumberOfBits(m_NumberOfBitsQ07);
		
		
		
		std::cout<<"The Quantization configuration: "<<std::endl;
		std::cout<<m_QuantRBIFilter->GetNumberOfBits()<<' ';
		std::cout<<m_QuantProjExtractionFilter->GetNumberOfBits()<<' ';
		std::cout<<m_QuantFPInputFilter->GetNumberOfBits()<<' ';
		std::cout<<m_QuantFPResultFilter->GetNumberOfBits()<<' ';
		std::cout<<m_QuantDiffResultFilter->GetNumberOfBits()<<' ';
		std::cout<<m_QuantB4NormFilter->GetNumberOfBits()<<' ';
		std::cout<<m_QuantBPInput->GetNumberOfBits()<<' ';
		std::cout<<m_QuantAdderInput->GetNumberOfBits()<<' ';
		std::cout<<std::endl;
		
		
		const unsigned int Dimension = this->InputImageDimension;

		// The backprojection works on one projection at a time
		typename ExtractFilterType::InputImageRegionType subsetRegion;
		subsetRegion = this->GetInput(1)->GetLargestPossibleRegion();
		unsigned int nProj = subsetRegion.GetSize(Dimension - 1);

		subsetRegion.SetSize(Dimension - 1, 1);// set particular dimension, here it means it takes only one.

		// Fill and shuffle randomly the projection order.
		// Should be tunable with other solutions.
		std::vector< unsigned int > projOrder(nProj);

		for (unsigned int i = 0; i < nProj; i++)
			projOrder[i] = i;
		//std::random_shuffle(projOrder.begin(), projOrder.end());

		m_MultiplyFilter->SetInput1((const float)m_Lambda / (double)m_NumberOfProjectionsPerSubset);

		// Create the zero projection stack used as input by RayBoxIntersectionFilter
		if (m_MyVerbose)
		{
			std::cout << "Before m_ConstantProjectionStackSource->Update(); in GenerateOutput" << std::endl; std::cin.get();
		}

		m_ConstantProjectionStackSource->Update();

		// Declare the image used in the main loop
		typename TInputImage::Pointer pimg;

		// For each iteration, go over each projection
		for (unsigned int iter = 0; iter < m_NumberOfIterations; iter++)
		{
			/*if (iter == 2)
			{
				m_MyWriter01->SetInput(m_RayBoxFilter->GetOutput());
				m_MyWriter01->SetFileName("Q01Input.mhd");
				m_MyWriter01->SetInput(m_QuantRBIFilter->GetOutput());
				m_MyWriter01->SetFileName("Q01_Output.mhd");
			}*/

			unsigned int projectionsProcessedInSubset = 0;

			for (unsigned int i = 0; i < nProj; i++)
			{
				if (i == 1)
				{
					m_MyWriter01->SetInput(m_RayBoxFilter->GetOutput());
					m_MyWriter01->SetFileName("Q01Input.mhd");
					m_MyWriter01->Update();
					m_MyWriter02->SetInput(m_QuantRBIFilter->GetOutput());
					m_MyWriter02->SetFileName("Q01_Output.mhd");
					m_MyWriter02->Update();
				}
				// When we reach the number of projections per subset:
				// - plug the output of the pipeline back into the Forward projection filter
				// - set the input of the Back projection filter to zero
				// - reset the projectionsProcessedInSubset to zero
				if (projectionsProcessedInSubset == m_NumberOfProjectionsPerSubset)
				{
					if (m_EnforcePositivity)
						pimg = m_ThresholdFilter->GetOutput();
					else
						pimg = m_AddFilter->GetOutput();

					pimg->DisconnectPipeline();
					m_QuantFPInputFilter->SetInput(pimg);// it only has one input
					// m_ForwardProjectionFilter->SetInput(1, pimg );
					m_AddFilter->SetInput2(pimg);
					m_BackProjectionFilter->SetInput(0, m_ConstantVolumeSource->GetOutput());

					projectionsProcessedInSubset = 0;
				}

				// Otherwise, just plug the output of the back projection filter
				// back as its input
				else
				{
					if (i)
					{
						// use the quantized result. Prevent the use of full precision data
						pimg = m_QuantAdderInput->GetOutput();
						pimg->DisconnectPipeline();
						m_BackProjectionFilter->SetInput(0, pimg);
						
						
						//pimg = m_BackProjectionFilter->GetOutput();
						//pimg->DisconnectPipeline();
						//m_BackProjectionFilter->SetInput(0, pimg);
					}
					else
					{
						m_BackProjectionFilter->SetInput(0, m_ConstantVolumeSource->GetOutput());
					}
				}

				// Change projection subset
				subsetRegion.SetIndex(Dimension - 1, projOrder[i]);
				m_ExtractFilter->SetExtractionRegion(subsetRegion);
				m_ExtractFilterRayBox->SetExtractionRegion(subsetRegion);

				// Set gating weight for the current projection
				if (m_IsGated)
				{
					m_GatingWeightsFilter->SetConstant2(m_GatingWeights[i]);
				}

				// This is required to reset the full pipeline
				if (m_MyVerbose)
				{
					std::cout << "Before m_BackProjectionFilter->UpdateOutputInformation(); in GenerateOutput" << std::endl; std::cin.get();
				}
				
				//m_BackProjectionFilter->GetOutput()->UpdateOutputInformation();
				//m_BackProjectionFilter->GetOutput()->PropagateRequestedRegion();
				
				m_QuantAdderInput->GetOutput()->UpdateOutputInformation();
				m_QuantAdderInput->GetOutput()->PropagateRequestedRegion();

				m_ExtractProbe.Start();
				if (m_MyVerbose)
				{
					std::cout << "Before RBI, Type Enter key" << std::endl; std::cin.get();
				}


				m_ExtractFilter->Update();
				m_ExtractFilterRayBox->Update();
				m_ExtractProbe.Stop();

				m_ZeroMultiplyProbe.Start();
				m_ZeroMultiplyFilter->Update();
				m_ZeroMultiplyProbe.Stop();
				if (m_MyVerbose)
				{
					std::cout << "Before FP, Type Enter key" << std::endl; std::cin.get();
				}


				m_ForwardProjectionProbe.Start();
				m_ForwardProjectionFilter->Update();
				m_ForwardProjectionProbe.Stop();
				if (m_MyVerbose)
				{
					std::cout << "After FP, Type Enter key" << std::endl; std::cin.get();
				}

				m_SubtractProbe.Start();
				m_SubtractFilter->Update();
				m_SubtractProbe.Stop();

				m_MultiplyProbe.Start();
				m_MultiplyFilter->Update();
				m_MultiplyProbe.Stop();

				m_RayBoxProbe.Start();
				m_RayBoxFilter->Update();
				m_RayBoxProbe.Stop();


				//m_RayBoxFilter
				m_Subtractor->SetInput1(m_RayBoxFilter->GetOutput());
				m_Subtractor->SetInput2(m_QuantRBIFilter->GetOutput());


				m_QuantRBIFilter->Update();


				m_MinMaxCalculator->SetImage(m_Subtractor->GetOutput());

				m_MinMaxCalculator->Compute();


				double max_pos_Diff = std::pow(2.0, -1.0*m_QuantRBIFilter->GetNumberOfBits());
				assert(m_MinMaxCalculator->GetMaximum() < max_pos_Diff);
				assert(std::abs(m_MinMaxCalculator->GetMinimum()) < max_pos_Diff);

				std::cout << "The maximum difference is " << m_MinMaxCalculator->GetMaximum() << ". 1/2^n = " << max_pos_Diff << std::endl;

				m_Subtractor = SubtractImageFilterType::New();
				m_MinMaxCalculator = MinMaxCalcType::New();

				m_DivideProbe.Start();
				m_DivideFilter->Update();
				m_DivideProbe.Stop();


			//	m_MyWriter->SetInput(m_QuantRBIFilter->GetOutput());
			//	m_MyWriter->Update();

				if (m_IsGated)
				{
					m_GatingProbe.Start();
					m_GatingWeightsFilter->Update();
					m_GatingProbe.Stop();
				}

				m_DisplacedDetectorProbe.Start();
				m_DisplacedDetectorFilter->Update();
				m_DisplacedDetectorProbe.Stop();
				if (m_MyVerbose)
				{
					std::cout << "Before BP, Type Enter key" << std::endl; std::cin.get();
				}

				m_BackProjectionProbe.Start();
				m_BackProjectionFilter->Update();
				m_BackProjectionProbe.Stop();
				
				m_QuantAdderInput->Update();
				if (m_MyVerbose)
				{
					std::cout << "After BP, Type Enter key" << std::endl; std::cin.get();
				}

				projectionsProcessedInSubset++;
				if ((projectionsProcessedInSubset == m_NumberOfProjectionsPerSubset) || (i == nProj - 1))
				{
					//m_AddFilter->SetInput1(m_BackProjectionFilter->GetOutput());
					m_AddFilter->SetInput1(m_QuantAdderInput->GetOutput());
					if (m_MyVerbose)
					{
						std::cout << "Before add, Type Enter key" << std::endl; std::cin.get();
					}

					m_AddProbe.Start();
					m_AddFilter->Update();
					m_AddProbe.Start();
					if (m_MyVerbose)
					{
						std::cout << "after add, Type Enter key" << std::endl; std::cin.get();
					}


					if (m_EnforcePositivity)
					{
						m_ThresholdProbe.Start();
						m_ThresholdFilter->Update();
						m_ThresholdProbe.Stop();
					}
				}

			}
		}
		if (m_EnforcePositivity)
		{
			this->GraftOutput(m_ThresholdFilter->GetOutput());
		}
		else
		{
			this->GraftOutput(m_AddFilter->GetOutput());
		}
	}

	template<class TInputImage, class TOutputImage>
	void
		SartWithZipmlFilter<TInputImage, TOutputImage>
		::PrintTiming(std::ostream & os) const
	{
		os << "SartWithZipmlFilter timing:" << std::endl;
		os << "  Extraction of projection sub-stacks: " << m_ExtractProbe.GetTotal()
			<< ' ' << m_ExtractProbe.GetUnit() << std::endl;
		os << "  Multiplication by zero: " << m_ZeroMultiplyProbe.GetTotal()
			<< ' ' << m_ZeroMultiplyProbe.GetUnit() << std::endl;
		os << "  Forward projection: " << m_ForwardProjectionProbe.GetTotal()
			<< ' ' << m_ForwardProjectionProbe.GetUnit() << std::endl;
		os << "  Subtraction: " << m_SubtractProbe.GetTotal()
			<< ' ' << m_SubtractProbe.GetUnit() << std::endl;
		os << "  Multiplication by lambda: " << m_MultiplyProbe.GetTotal()
			<< ' ' << m_MultiplyProbe.GetUnit() << std::endl;
		os << "  Ray box intersection: " << m_RayBoxProbe.GetTotal()
			<< ' ' << m_RayBoxProbe.GetUnit() << std::endl;
		os << "  Division: " << m_DivideProbe.GetTotal()
			<< ' ' << m_DivideProbe.GetUnit() << std::endl;
		os << "  Multiplication by the gating weights: " << m_GatingProbe.GetTotal()
			<< ' ' << m_GatingProbe.GetUnit() << std::endl;
		os << "  Displaced detector: " << m_DisplacedDetectorProbe.GetTotal()
			<< ' ' << m_DisplacedDetectorProbe.GetUnit() << std::endl;
		os << "  Back projection: " << m_BackProjectionProbe.GetTotal()
			<< ' ' << m_BackProjectionProbe.GetUnit() << std::endl;
		os << "  Volume update: " << m_AddProbe.GetTotal()
			<< ' ' << m_AddProbe.GetUnit() << std::endl;
		if (m_EnforcePositivity)
		{
			os << "  Positivity enforcement: " << m_ThresholdProbe.GetTotal()
				<< ' ' << m_ThresholdProbe.GetUnit() << std::endl;
		}
	}

} // end namespace rtk

#endif // rtkSartWithZipmlFilter_hxx
