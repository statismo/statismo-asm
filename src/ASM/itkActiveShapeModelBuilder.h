/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Langguth (christoph.langguth@unibas.ch)
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

#ifndef ITKACTIVESHAPEMODELBUILDER_H_
#define ITKACTIVESHAPEMODELBUILDER_H_

#include <itkObject.h>
#include <itkObjectFactory.h>

#include "itkActiveShapeModel.h"
#include "ActiveShapeModelBuilder.h"
#include "itkASMOps.h"

namespace itk {

/**
 * \brief ITK Wrapper for the statismo::PCAModelBuilder class.
 * \see statismo::PCAModelBuilder for detailed documentation.
 */
    template <typename MeshType, typename ImageType>
    class ActiveShapeModelBuilder : public Object {
    public:

        typedef ActiveShapeModelBuilder            Self;
        typedef Object	Superclass;
        typedef SmartPointer<Self>                Pointer;
        typedef SmartPointer<const Self>          ConstPointer;

        typedef ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;

        itkNewMacro( Self );
        itkTypeMacro( ActiveShapeModelBuilder, Object );


        typedef statismo::ActiveShapeModelBuilder<MeshType, ImageType> ImplType;
        typedef typename ImplType::FeatureExtractorPointerType FeatureExtractorPointerType;
        typedef typename ImplType::ImagePreprocessorPointerType ImagePreprocessorPointerType;
        typedef typename ImplType::DataPairType DataPairType;
        typedef typename ImplType::DataPairsType DataPairsType;

        ActiveShapeModelBuilder() : m_impl(ImplType::Create()) {}

        virtual ~ActiveShapeModelBuilder() {
            if (m_impl) {
                delete m_impl;
                m_impl = 0;
            }
        }

        typename ActiveShapeModelType::Pointer BuildNewModel(const typename ActiveShapeModelType::ImplType::RepresenterType* representer, const DataPairsType& datapairs, const FeatureExtractorPointerType featureExtractor, const ImagePreprocessorPointerType imagePreprocessor, const statismo::ASMPointIdSampler* profileSampler) {
            statismo::ActiveShapeModel<MeshType, ImageType>* inner = m_impl->BuildNewModel(representer, datapairs, featureExtractor, imagePreprocessor, profileSampler);
            typename ActiveShapeModelType::Pointer outer = ActiveShapeModelType::New();
            outer->SetstatismoImplObj(inner);
            typename ActiveShapeModelType::StatisticalModelPointerType ssm = ActiveShapeModelType::StatisticalModelType::New();
            statismo::StatisticalModel<MeshType>* issm = const_cast<statismo::StatisticalModel<MeshType>*>(inner->GetStatisticalModel());
            ssm->SetstatismoImplObj(issm);
            outer->SetStatisticalModel(ssm);
            return outer;
        }




    private:
        ActiveShapeModelBuilder(const ActiveShapeModelBuilder& orig);
        ActiveShapeModelBuilder& operator=(const ActiveShapeModelBuilder& rhs);

        ImplType* m_impl;
    };


}

#endif //ITKACTIVESHAPEMODELBUILDER_H_

