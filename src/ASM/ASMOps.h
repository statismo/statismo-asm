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

#ifndef __ASM_OPERATIONS
#define __ASM_OPERATIONS

#include "Representer.h"

namespace statismo {


    template<class T>
    struct ASMOperationsTraits {
    };

    template <class T>
    class ASMOperations {
    public:

        typedef typename Representer<T>::DatasetType  MeshType;
        typedef typename Representer<T>::DatasetPointerType  MeshPointerType;
        typedef typename Representer<T>::PointType PointType;
        typedef typename ASMOperationsTraits<T>::RigidTransformType RigidTransformType;
        typedef typename ASMOperationsTraits<T>::RigidTransformPointerType RigidTransformPointerType;


        virtual PointType
        TransformPoint(PointType &point, const RigidTransformPointerType transform, bool inverse) const = 0;

        virtual MeshPointerType
        TransformMesh(MeshType* mesh, const RigidTransformPointerType transform) const = 0;

        virtual RigidTransformPointerType
        ComputeRigidTransformFromLandmarks(const std::vector <PointType> &fixedLandmarks,
                                           const std::vector <PointType> &movingLandmarks) const = 0;

    };

}


#endif