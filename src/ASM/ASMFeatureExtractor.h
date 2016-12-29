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

#ifndef STATISMO_ASMFEATUREEXTRACTOR_H
#define STATISMO_ASMFEATUREEXTRACTOR_H

#include "CommonTypes.h"
#include "Representer.h"
#include "ASMImagePreprocessor.h"
#include "ASMOps.h"
//#include "ActiveShapeModel.h"

namespace H5 {
    class Group;
}

namespace statismo {
    //forward declaration
    template <typename TPointSet, typename TImage> class ActiveShapeModel;

    template <typename TPointSet, typename TImage>
    class ASMFeatureExtractor {
            typedef typename Representer<TPointSet>::DatasetType MeshType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
    typedef typename Representer<TPointSet>::PointType PointType;
    typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef typename statismo::ASMOperations<TPointSet>::RigidTransformPointerType RigidTransformPointerType;
    typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;


    public:
        virtual ~ASMFeatureExtractor(){}
        virtual void Delete() = 0;
        virtual FeatureExtractorType* Clone() const = 0;
        virtual FeatureExtractorType* CloneForTarget(const ActiveShapeModelType* const model, const VectorType& coefficients, const MeshType* mesh) const = 0;
        virtual FeatureExtractorType* CloneForMesh(TPointSet* mesh) const = 0;
        virtual bool ExtractFeatures(statismo::VectorType& output, const PreprocessedImageType* const image, const PointType& point) const = 0;
        virtual void Save(const H5::Group& h5Group) const = 0;
    };

    template<typename TPointSet, typename TImage>
    class ASMFeatureExtractorFactory {
    typedef ASMFeatureExtractorFactory<TPointSet, TImage> ASMFeatureExtractorFactoryType;
    private:
        static std::vector<const ASMFeatureExtractorFactoryType*> *implementations() {
            static std::vector<const ASMFeatureExtractorFactoryType*> impls;
            return &impls;
        }
        ASMFeatureExtractorFactory(const ASMFeatureExtractorFactoryType& o) { }
        ASMFeatureExtractorFactory& operator=(const ASMFeatureExtractorFactoryType& o) {}

    protected:
        ASMFeatureExtractorFactory() {}

    public:
        virtual std::string GetDescriptor() const = 0;
        virtual const ASMFeatureExtractor<TPointSet, TImage>* Instantiate(const H5::Group& h5Group) const = 0;

        static std::vector<const ASMFeatureExtractorFactoryType*> GetImplementations() {
            return *implementations();
        }

        static void RegisterImplementation(const ASMFeatureExtractorFactoryType* impl) {
            if (GetImplementation(impl->GetDescriptor())) {
                //ignoring already-registered implementation
                return;
            }
            implementations()->push_back(impl);
        }

        static const ASMFeatureExtractorFactoryType* GetImplementation(std::string descriptor) {
            std::vector<const ASMFeatureExtractorFactoryType* > impls = GetImplementations();
            for (typename std::vector<const ASMFeatureExtractorFactoryType* >::iterator impl = impls.begin(); impl != impls.end(); ++impl) {
                if ((*impl)->GetDescriptor() == descriptor) {
                    return *impl;
                }
            }
            return 0;
        }
    };
}

#endif //STATISMO_ASMFEATUREEXTRACTOR_H
