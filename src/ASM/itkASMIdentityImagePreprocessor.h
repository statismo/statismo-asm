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

#ifndef STATISMO_ITKASMIDENTITYIMAGEPREPROCESSOR_H
#define STATISMO_ITKASMIDENTITYIMAGEPREPROCESSOR_H

#include "ASMImagePreprocessor.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "HDF5Utils.h"

namespace itk {

    template <typename TPointSet, typename TImage>
    class ASMIdentityPreprocessedImage: public statismo::ASMPreprocessedImage<TPointSet> {
    private:
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;
        //typedef typename InterpolatedImageType::CovariantVectorType CovariantVectorType;
        typedef CovariantVector<float, TImage::ImageDimension> CovariantVectorType;
        typedef vnl_vector<statismo::ScalarType> VectorType;


        const TImage* m_inputImage;
        const InterpolatedImagePointerType m_interpolatedImage;

        ASMIdentityPreprocessedImage(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage): m_inputImage(inputImage), m_interpolatedImage(interpolatedImage) {}

        static statismo::VectorType fromVnlVector(const VectorType& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

        }
    public:

        virtual typename statismo::ASMPreprocessedImage<TPointSet>::ImageType GetType() const  {
            return statismo::ASMPreprocessedImage<TPointSet>::ImageType::Intensity;
        };


        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        virtual bool IsDefined(const PointType& point) const {
            typename TImage::IndexType index;
            // we don't care about the actual index, just whether it's present
            return m_inputImage->TransformPhysicalPointToIndex(point, index);
        }

        virtual statismo::VectorType Evaluate(const PointType& point) const {
            CovariantVectorType cv = m_interpolatedImage->Evaluate(point);
            return fromVnlVector(cv.GetVnlVector());
        }

        static ASMIdentityPreprocessedImage* Create(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage) {
            return new ASMIdentityPreprocessedImage(inputImage, interpolatedImage);
        }


    };

    template <typename TPointSet, typename TImage>
    class ASMIdentityImagePreprocessor: public statismo::ASMImagePreprocessor<TPointSet, TImage> {
    private:
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;


        InterpolatedImagePointerType Interpolate(const TImage* const image) const {
            InterpolatedImagePointerType inter = InterpolatedImageType::New();
            inter->SetSplineOrder(1);

            inter->SetInputImage(image);
            return inter;
        }

    public:
        typedef ASMIdentityPreprocessedImage<TPointSet, TImage> PreprocessedImplType;

        ASMIdentityImagePreprocessor() {}


        virtual ASMIdentityImagePreprocessor<TPointSet, TImage>* Clone() const {
            return new ASMIdentityImagePreprocessor();
        };

        virtual PreprocessedImplType* Preprocess(const TImage* image) const {
            return PreprocessedImplType::Create(image, Interpolate(image));
        };

        void Save(const H5::Group& h5Group) const {}

    };

    template<typename TPointSet, typename TImage>
    class ASMIdentityImagePreprocessorFactory : public statismo::ASMImagePreprocessorFactory<TPointSet, TImage> {

        typedef ASMIdentityImagePreprocessor<TPointSet, TImage> InstanceType;

    public:

        static const ASMIdentityImagePreprocessorFactory *GetInstance() {
            static ASMIdentityImagePreprocessorFactory *instance = new ASMIdentityImagePreprocessorFactory();
            return instance;
        }

        virtual std::string GetDescriptor() const {
            return "builtin::Identity";
        }

        virtual const statismo::ASMImagePreprocessor<TPointSet, TImage> *Instantiate(
                const H5::Group &h5Group) const {

            InstanceType* instance = new InstanceType();
            return instance;
        }
    };


}
#endif //STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H
