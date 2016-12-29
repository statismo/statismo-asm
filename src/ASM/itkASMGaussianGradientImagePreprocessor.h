/*
 * This file is part of the statismo library.
 *
 *
 * Copyright (c) 2011-2015 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H
#define STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H

#include "ASMImagePreprocessor.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "HDF5Utils.h"

namespace itk {

    template <typename TPointSet, typename TImage>
    class ASMGaussianGradientPreprocessedImage: public statismo::ASMPreprocessedImage<TPointSet> {
    private:
        typedef DiscreteGaussianImageFilter<TImage, TImage> GaussianFilterType;
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;
        //typedef typename InterpolatedImageType::CovariantVectorType CovariantVectorType;
        typedef CovariantVector<float, TImage::ImageDimension> CovariantVectorType;
        typedef vnl_vector<statismo::ScalarType> VectorType;


        const TImage* m_inputImage;
        const InterpolatedImagePointerType m_interpolatedImage;

        ASMGaussianGradientPreprocessedImage(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage): m_inputImage(inputImage), m_interpolatedImage(interpolatedImage) {}

        static statismo::VectorType fromVnlVector(const VectorType& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

        }
    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;



        virtual bool IsDefined(const PointType& point) const {
            typename TImage::IndexType index;
            // we don't care about the actual index, just whether it's present
            return m_inputImage->TransformPhysicalPointToIndex(point, index);
        }

        virtual statismo::VectorType Evaluate(const PointType& point) const {
            CovariantVectorType cv = m_interpolatedImage->EvaluateDerivative(point);
            return fromVnlVector(cv.GetVnlVector());
        }

        static ASMGaussianGradientPreprocessedImage* Create(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage) {
            return new ASMGaussianGradientPreprocessedImage(inputImage, interpolatedImage);
        }

        virtual typename statismo::ASMPreprocessedImage<TPointSet>::ImageType GetType() const {
                return statismo::ASMPreprocessedImage<TPointSet>::Gradient;
        };

    };

    template <typename TPointSet, typename TImage>
    class ASMGaussianGradientImagePreprocessor: public statismo::ASMImagePreprocessor<TPointSet, TImage> {
    private:
        typedef DiscreteGaussianImageFilter<TImage, TImage> GaussianFilterType;
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;

        const float m_sigma;

        InterpolatedImagePointerType Interpolate(const TImage* const image) const {
            InterpolatedImagePointerType inter = InterpolatedImageType::New();
            inter->SetSplineOrder(1);
            if (m_sigma != 0) {
                typename GaussianFilterType::Pointer smooth = GaussianFilterType::New();
                smooth->SetVariance(m_sigma * m_sigma);
                // FIXME: smooth->SetMaximumKernelWidth ???
                smooth->SetInput(image);
                smooth->Update();
                inter->SetInputImage(smooth->GetOutput());
            } else {
                inter->SetInputImage(image);
            }

            return inter;
        }

    public:
        typedef ASMGaussianGradientPreprocessedImage<TPointSet, TImage> PreprocessedImplType;

        ASMGaussianGradientImagePreprocessor(float sigma): m_sigma(sigma) {}

        virtual ASMGaussianGradientImagePreprocessor<TPointSet, TImage>* Clone() const {
            return new ASMGaussianGradientImagePreprocessor(m_sigma);
        };

        virtual PreprocessedImplType* Preprocess(const TImage* image) const {
            return PreprocessedImplType::Create(image, Interpolate(image));
        };

        virtual void Save(const H5::Group &h5Group) const {
            statismo::HDF5Utils::writeStringAttribute(h5Group, "identifier", "builtin::GaussianGradient");
            statismo::HDF5Utils::writeIntAttribute(h5Group, "majorVersion", 1);
            statismo::HDF5Utils::writeIntAttribute(h5Group, "minorVersion", 0);

            statismo::HDF5Utils::writeFloat(h5Group, "stddev", m_sigma);
        }
    };

    template<typename TPointSet, typename TImage>
    class ASMGaussianGradientImagePreprocessorFactory : public statismo::ASMImagePreprocessorFactory<TPointSet, TImage> {

        typedef ASMGaussianGradientImagePreprocessor<TPointSet, TImage> InstanceType;

    public:

        static const ASMGaussianGradientImagePreprocessorFactory *GetInstance() {
            static ASMGaussianGradientImagePreprocessorFactory *instance = new ASMGaussianGradientImagePreprocessorFactory();
            return instance;
        }

        virtual std::string GetDescriptor() const {
            return "builtin::GaussianGradient";
        }

        virtual const statismo::ASMImagePreprocessor<TPointSet, TImage> *Instantiate(
                const H5::Group &h5Group) const {

            float sigma = statismo::HDF5Utils::readFloat(h5Group, "stddev");
            InstanceType* instance = new InstanceType(sigma);
            return instance;
        }
    };


}
#endif //STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H
