#ifndef __itkQuantizationFilter_h
#define __itkQuantizationFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
template< class TImage>
class QuantizationFilter : public ImageToImageFilter< TImage, TImage >
{
public:
  /** Standard class typedefs. */
  typedef QuantizationFilter             Self;
  typedef ImageToImageFilter< TImage, TImage > Superclass;
  typedef SmartPointer< Self >        Pointer;

  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageFilter, ImageToImageFilter);

  //itkSetMacro(Pitch, typename TImage::ValueType);
  itkSetMacro(NumberOfBits, int);
	itkGetMacro(NumberOfBits, int);
  

protected:
	QuantizationFilter();// do not use the default constructor  QuantizationFilter(){}  -  QuantizationFilter();
  ~QuantizationFilter(){}

  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData(const OutputImageRegionType &,
                                    ThreadIdType);

private:
  QuantizationFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  typename TImage::ValueType m_Maximum, m_Minimum;
  double m_Pitch; // higher precision
  int m_NumberOfBits;



};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "QuantizationFilter.hxx"
#endif


#endif // __itkQuantizationFilter_h
