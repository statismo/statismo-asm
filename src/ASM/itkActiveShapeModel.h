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

#ifndef STATISMO_ITKACTIVESHAPEMODEL_H
#define STATISMO_ITKACTIVESHAPEMODEL_H

#include <itkObject.h>
#include <itkMacro.h>
#include "itkStatisticalModel.h"
#include "ActiveShapeModel.h"


namespace itk {

    template <typename TPointSet, typename TImage>
    class ActiveShapeModel: public Object {
    public:
        /* standard ITK typedefs and macros */
        typedef ActiveShapeModel Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        itkNewMacro( Self );
        itkTypeMacro( Self, Object);

        typedef statismo::ActiveShapeModel<TPointSet, TImage> ImplType;
        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef typename StatisticalModelType::PointType PointType;
        typedef typename StatisticalModelType::Pointer StatisticalModelPointerType;
        typedef typename StatisticalModelType::ImplType StatisticalModelImplType;

        template <class F>
        typename boost::result_of<F()>::type callstatismoImpl(F f) const {
            try {
                return f();
            } catch (statismo::StatisticalModelException& s) {
                itkExceptionMacro(<< s.what());
            }
        }

        virtual void SetstatismoImplObj(ImplType* impl) {
            if (m_impl && m_impl != impl) {
                delete m_impl;
            }
            m_impl = impl;
        }

        virtual ImplType* GetstatismoImplObj() const {
            return m_impl;
        }

        ActiveShapeModel() : m_impl(0), m_ssm(0) {}

        virtual ~ActiveShapeModel() {
            if (m_impl) {
                delete m_impl;
            }
        }

        const StatisticalModelPointerType GetStatisticalModel() {
            return m_ssm;
        }

        void SetStatisticalModel(StatisticalModelPointerType ssm) {
            m_ssm = ssm;
            // update implementation to contain the new (statismo) SSM.
            // We must use a copy because otherwise we'd have two objects sharing a pointer,
            // both wanting to delete the pointer in their own destructor.
            StatisticalModelImplType* m = ssm->GetstatismoImplObj();
            StatisticalModelImplType* impl_copy = StatisticalModelImplType::Create(m->GetRepresenter(), m->GetMeanVector(), m->GetOrthonormalPCABasisMatrix(), m->GetPCAVarianceVector(), m->GetNoiseVariance());
            ImplType* newImpl = GetstatismoImplObj()->CloneWithStatisticalModel(impl_copy);
            SetstatismoImplObj(newImpl);
        }

        void Load(typename ImplType::RepresenterType* representer, const char* filename) {
            try {
                SetstatismoImplObj(ImplType::Load(representer, filename));
                m_ssm = StatisticalModelType::New();
                // same logic as above, we need a copy, not a shared pointer.
                StatisticalModelImplType* m = const_cast<StatisticalModelImplType*>(m_impl->GetStatisticalModel());
                StatisticalModelImplType* impl_copy = StatisticalModelImplType::Create(m->GetRepresenter(), m->GetMeanVector(), m->GetOrthonormalPCABasisMatrix(), m->GetPCAVarianceVector(), m->GetNoiseVariance());
                m_ssm->SetstatismoImplObj(impl_copy);
            } catch (statismo::StatisticalModelException& s) {
                itkExceptionMacro(<< s.what());
            }
        }

        void Save(const std::string& filename) {
            GetstatismoImplObj()->Save(filename);
        }


    private:
        ImplType* m_impl;
        StatisticalModelPointerType m_ssm;
    };

}

#endif //STATISMO_ITKACTIVESHAPEMODEL_H
