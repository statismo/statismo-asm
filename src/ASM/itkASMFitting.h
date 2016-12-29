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

#ifndef STATISMO_ITKASMFITTING_H
#define STATISMO_ITKASMFITTING_H

#include <itkObject.h>
#include <itkMacro.h>
#include "ASMFitting.h"
#include "itkActiveShapeModel.h"
#include "itkASMPointSampler.h"
#include "itkASMOps.h"

namespace itk {

    template<typename TPointSet, typename TImage>
    class ASMFittingResult : public Object {
    public:
        typedef ASMFittingResult Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelType;
        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelPointerType;
        typedef typename ASMMeshOperations::RigidTransformPointerType RigidTransformPointerType;

        itkNewMacro(Self);
        itkTypeMacro(Self, Object);

        ~ASMFittingResult() {}

        typedef statismo::ASMFittingResult<RigidTransformPointerType> ImplType;
        typedef vnl_vector<statismo::ScalarType> VectorType;

        void SetInternalData(std::shared_ptr<ASMMeshOperations> asmMeshOperations, ImplType statismoResult, ModelPointerType model) {
            m_meshOperations = asmMeshOperations;
            m_model = model;
            m_statismoResult = statismoResult;
        }

        bool IsValid() {
            return m_statismoResult.IsValid();
        }

        VectorType GetCoefficients() {
            return toVnlVector(m_statismoResult.GetCoefficients());
        }

        RigidTransformPointerType GetRigidTransformation() {
            return m_statismoResult.GetRigidTransform();
        }

        typename TPointSet::Pointer GetMesh() {
            typename TPointSet::Pointer instance = m_model->GetStatisticalModel()->DrawSample(GetCoefficients());
            return m_meshOperations->TransformMesh(instance, GetRigidTransformation());
        }


    private:
        std::shared_ptr<ASMMeshOperations> m_meshOperations;
        ImplType m_statismoResult;
        ModelPointerType m_model;

        VectorType toVnlVector(const statismo::VectorType& v) {
            return VectorType(v.data(), v.rows());

        }
    };

    template<typename TPointSet, typename TImage>
    class ASMFittingStep : public Object {
    public:
        typedef ASMFittingStep Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(Self, Object);

        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelPointerType;
        typedef typename ASMMeshOperations::RigidTransformPointerType RigidTransformPointerType;
        typedef typename TPointSet::Pointer PointSetPointerType;
        typedef typename statismo::ASMPreprocessedImage<TPointSet> *ImagePointerType;
        typedef statismo::ASMFittingConfiguration ConfigurationType;
        typedef typename ASMPointSampler<TPointSet, TImage>::Pointer SamplerPointerType;
        typedef statismo::ASMFittingStep<TPointSet, TImage> ImplType;
        typedef ASMFittingResult<TPointSet, TImage> ResultType;


        ASMFittingStep() : m_asmOperations(new ASMMeshOperations()), m_model(0), m_target(0), m_configuration(0, 0, 0), m_transform(0) { }

        void SetModel(ModelPointerType model) {
            m_model = model;
        }

        void SetCoefficients(statismo::VectorType coeffs) {
            m_coeffs = coeffs;
        }

        void SetRigidTransformation(RigidTransformPointerType transform) {
            m_transform = transform;
        }

        void SetTarget(ImagePointerType target) {
            m_target = target;
        }

        void SetSampler(SamplerPointerType sampler) {
            m_sampler = sampler;
        }

        void SetConfiguration(const ConfigurationType &configuration) {
            m_configuration = configuration;
        }

        void Update() {
            ImplType *impl = ImplType::Create(m_asmOperations.get(), m_configuration, m_model->GetstatismoImplObj(), m_coeffs, m_transform,
                                              m_target, m_sampler);

            statismo::ASMFittingResult<RigidTransformPointerType> result = impl->Perform();
            m_result = ResultType::New();
            m_result->SetInternalData(m_asmOperations, result, m_model);

            delete impl;
        }

        typename ResultType::Pointer GetOutput() {
            return m_result;
        }

    private:
        std::shared_ptr<ASMMeshOperations> m_asmOperations;
        ModelPointerType m_model;
        statismo::VectorType m_coeffs;
        RigidTransformPointerType m_transform;
        ImagePointerType m_target;
        SamplerPointerType m_sampler;
        ConfigurationType m_configuration;
        typename ResultType::Pointer m_result;
    };

    template<typename TPointSet, typename TImage>
    class ASMFitting : public Object {
    public:
        typedef ASMFitting Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(Self, Object);
    };
}
#endif //STATISMO_ITKASMFITTING_H
