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


#ifndef __ITK_ASM_MESH_OPERATIONS
#define __ITK_ASM_MESH_OPERATIONS

#include "Representer.h"
#include "ASMOps.h"

#include <itkIdentityTransform.h>
#include <itkPoint.h>
#include <itkTransformMeshFilter.h>
#include <itkVector.h>
#include <itkLandmarkBasedTransformInitializer.h>
#include <itkVersorRigid3DTransform.h>
#include <itkImage.h>


namespace statismo {


    template<>
    struct ASMOperationsTraits<itk::Mesh <float, 3u> > {
            typedef itk::VersorRigid3DTransform<float> RigidTransformType;
            typedef RigidTransformType::Pointer RigidTransformPointerType;
    };

}

namespace itk {


    class ASMMeshOperations : public statismo::ASMOperations<itk::Mesh<float, 3u> > {
    public:

        static const unsigned MeshDimension = 3u;
        typedef itk::Mesh<float, MeshDimension> MeshType;
        typedef statismo::ASMOperationsTraits<MeshType>::RigidTransformType RigidTransformType;
        typedef statismo::ASMOperationsTraits<MeshType>::RigidTransformPointerType RigidTransformPointerType;


        virtual PointType
        TransformPoint(PointType &point, const RigidTransformPointerType transform, bool inverse) const {
            if (transform) {
                if (inverse) {
                    return transform->GetInverseTransform()->TransformPoint(point);
                } else {
                    return transform->TransformPoint(point);
                }
            } else {
                return point;
            }
        }

        virtual MeshPointerType
        TransformMesh(MeshType* mesh, const RigidTransformPointerType transform) const {
            if (transform) {
                typedef itk::TransformMeshFilter<MeshType, MeshType, RigidTransformType> TransformMeshFilterType;

                typename TransformMeshFilterType::Pointer tf = TransformMeshFilterType::New();
                tf->SetInput(mesh);
                tf->SetTransform(transform);
                tf->Update();

                typename MeshType::Pointer output = tf->GetOutput();
                output->DisconnectPipeline();
                return output;
            } else {
                return mesh;
            }
        }

        virtual RigidTransformPointerType
        ComputeRigidTransformFromLandmarks(const std::vector <PointType> &fixedLandmarks,
                                           const std::vector <PointType> &movingLandmarks) const {

            // initialize the rigid transform
            typedef itk::Image<float, MeshDimension> DistanceImageType;
            typedef itk::VersorRigid3DTransform<float> RigidTransformType;
            typedef itk::LandmarkBasedTransformInitializer<RigidTransformType, DistanceImageType, DistanceImageType> LandmarkTransformInitializerType;
            typedef itk::Point<double, MeshDimension> LandmarkPointType;

            //FIXME: YUCK. There must be a better way around this.
            std::vector<LandmarkPointType> fixedDoubleLandmarks;
            std::vector<LandmarkPointType> movingDoubleLandmarks;

            for (typename std::vector<PointType>::const_iterator lm = fixedLandmarks.begin(); lm != fixedLandmarks.end(); ++lm) {
                LandmarkPointType lmd;
                for (int i=0; i < (*lm).Dimension; ++i) { lmd[i] = (*lm)[i];}
                fixedDoubleLandmarks.push_back(lmd);
            }

            for (typename std::vector<PointType>::const_iterator lm = movingLandmarks.begin(); lm != movingLandmarks.end(); ++lm) {
                LandmarkPointType lmd;
                for (int i=0; i < (*lm).Dimension; ++i) { lmd[i] = (*lm)[i];}
                movingDoubleLandmarks.push_back(lmd);
            }


            RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
            typename LandmarkTransformInitializerType::Pointer initializer = LandmarkTransformInitializerType::New();
            initializer->SetFixedLandmarks(fixedDoubleLandmarks);
            initializer->SetMovingLandmarks(movingDoubleLandmarks);
            initializer->SetTransform(rigidTransform);
            initializer->InitializeTransform();

            return RigidTransformPointerType(rigidTransform.GetPointer());
        }


    };

}


#endif