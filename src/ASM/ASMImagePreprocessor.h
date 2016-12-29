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

#ifndef STATISMO_ASMIMAGEPREPROCESSOR_H
#define STATISMO_ASMIMAGEPREPROCESSOR_H

#include "CommonTypes.h"
#include "Representer.h"

namespace statismo {

    template <typename TPointSet>
    class ASMPreprocessedImage {


    public:

        enum ImageType {
            Intensity = 0,
            Gradient = 1
        };

        typedef typename Representer<TPointSet>::PointType PointType;
        virtual ~ASMPreprocessedImage() {}
        virtual bool IsDefined(const PointType& point) const = 0;
        virtual VectorType Evaluate(const PointType& point) const = 0;
        virtual ImageType GetType() const = 0;
    };

    template <typename TPointSet, typename TImage>
    class ASMImagePreprocessor {
    public:
        virtual ~ASMImagePreprocessor() {}
        virtual ASMImagePreprocessor<TPointSet, TImage>* Clone() const = 0;
        virtual ASMPreprocessedImage<TPointSet>* Preprocess(const TImage* image) const = 0;
        virtual void Save(const H5::Group& h5Group) const = 0;
    };


    template<typename TPointSet, typename TImage>
    class ASMImagePreprocessorFactory {
        typedef ASMImagePreprocessorFactory<TPointSet, TImage> ASMImagePreprocessorFactoryType;
    private:
        static std::vector<const ASMImagePreprocessorFactoryType*> *implementations() {
            static std::vector<const ASMImagePreprocessorFactoryType*> impls;
            return &impls;
        }
        ASMImagePreprocessorFactory(const ASMImagePreprocessorFactoryType& o) { }
        ASMImagePreprocessorFactory& operator=(const ASMImagePreprocessorFactoryType& o) {}

    protected:
        ASMImagePreprocessorFactory() {}

    public:
        virtual std::string GetDescriptor() const = 0;
        virtual const ASMImagePreprocessor<TPointSet, TImage>* Instantiate(const H5::Group& h5Group) const = 0;

        static std::vector<const ASMImagePreprocessorFactoryType*> GetImplementations() {
            return *implementations();
        }

        static void RegisterImplementation(const ASMImagePreprocessorFactoryType* impl) {
            if (GetImplementation(impl->GetDescriptor())) {
                //ignoring already-registered implementation
                return;
            }
            implementations()->push_back(impl);
        }

        static const ASMImagePreprocessorFactoryType* GetImplementation(std::string descriptor) {
            std::vector<const ASMImagePreprocessorFactoryType* > impls = GetImplementations();
            for (typename std::vector<const ASMImagePreprocessorFactoryType* >::iterator impl = impls.begin(); impl != impls.end(); ++impl) {
                if ((*impl)->GetDescriptor() == descriptor) {
                    return *impl;
                }
            }
            return 0;
        }
    };

}
#endif //STATISMO_ASMIMAGEPREPROCESSOR_H
