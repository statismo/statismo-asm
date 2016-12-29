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

#ifndef STATISMO_ACTIVESHAPEMODELBUILDER_H
#define STATISMO_ACTIVESHAPEMODELBUILDER_H

#include "DeferredItem.h"
#include "ASMFeatureExtractor.h"
#include "ASMPointSampler.h"
#include "ASMImagePreprocessor.h"
#include "PCAModelBuilder.h"
#include "ASMPointIdSampler.h"

namespace statismo {

    template<typename MeshType, typename ImageType>
    class ActiveShapeModelBuilder {

    public:

        typedef std::pair<std::string,std::pair<DeferredItem<MeshType*>*,DeferredItem<ImageType*>*> > DataPairType;
        typedef std::vector<DataPairType> DataPairsType;

        typedef ASMFeatureExtractor<MeshType, ImageType>* FeatureExtractorPointerType;
        typedef ASMImagePreprocessor<MeshType, ImageType>* ImagePreprocessorPointerType;

        typedef ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;

        /**
         * Factory method to create a new ActiveShapeModelBuilder
         */
        static ActiveShapeModelBuilder *Create() {
            return new ActiveShapeModelBuilder();
        }

        /**
         * Destroy the object.
         * The same effect can be achieved by deleting the object in the usual
         * way using the c++ delete keyword.
         */
        void Delete() {
            delete this;
        }


        /**
         * The destructor
         */
        virtual ~ActiveShapeModelBuilder() { }

        ActiveShapeModelType* BuildNewModel(const typename ActiveShapeModelType::RepresenterType* representer, const DataPairsType& datapairs, const FeatureExtractorPointerType featureExtractor, const ImagePreprocessorPointerType imagePreprocessor, const ASMPointIdSampler* profileSampler) {
            std::vector<MeshType*> meshes;

            typedef std::vector<VectorType> ProfilesDataType;
            typedef std::pair<unsigned int, ProfilesDataType> ProfilesIdAndDataType;
            std::vector<ProfilesIdAndDataType> profileData;

            DataManager<MeshType>* dataManager = DataManager<MeshType>::Create(representer);

            // FIXME!!! This is ITK-specific code, and must be replaced by a real (uniform) sampler anyway
            std::vector<unsigned int> profilePoints = profileSampler->SamplePoints();

            for (std::vector<unsigned int>::const_iterator it = profilePoints.begin(); it != profilePoints.end(); it++) {
                ProfilesDataType data;
                profileData.push_back(std::make_pair((*it), data));
            }

            std::cout << "Building Active Shape Model from " << datapairs.size() << " datasets, using " << profilePoints.size() << " profile points. This can take a long time." << std::endl;

            for (typename DataPairsType::const_iterator pair = datapairs.begin(); pair != datapairs.end(); pair++) {
                std::cout << "Processing: " << (*pair).first << " ..." << std::endl;

                MeshType* mesh = (*pair).second.first->Get();
                ImageType* image = (*pair).second.second->Get();

                ASMPreprocessedImage<MeshType>* preprocessed = imagePreprocessor->Preprocess(image);

                FeatureExtractorPointerType fe = featureExtractor->CloneForMesh(mesh);

                dataManager->AddDataset(mesh, (*pair).first);
                meshes.push_back(mesh);

                for (std::vector<ProfilesIdAndDataType>::iterator pd = profileData.begin(); pd != profileData.end(); pd++) {
                    VectorType features;
                    typename ActiveShapeModelType::PointType point = mesh->GetPoint((*pd).first);
                    bool ok = fe->ExtractFeatures(features, preprocessed, point);
                    if (ok) {
                        (*pd).second.push_back(features);
                    }
                }

                image->Delete();
                delete preprocessed;
                fe->Delete();
            }

            unsigned int dataCount = datapairs.size();

            std::vector<ASMProfile> profiles;
            for (std::vector<ProfilesIdAndDataType>::const_iterator pd = profileData.begin(); pd != profileData.end(); pd++) {
                if ((*pd).second.size() != dataCount) {
                    std::cout << "Not creating profile for point id " << (*pd).first << ": not all input data produced usable features" <<std::endl;
                } else {
                    MultiVariateNormalDistribution distribution = MultiVariateNormalDistribution::EstimateFromData((*pd).second);
                    ASMProfile profile((*pd).first, distribution);
                    profiles.push_back(profile);
                }
            }

            PCAModelBuilder<MeshType>* ssmBuilder = PCAModelBuilder<MeshType>::Create();
            StatisticalModel<MeshType>* ssm = ssmBuilder->BuildNewModel(dataManager->GetData(), 0);

            dataManager->Delete();
            ssmBuilder->Delete();

            for(typename std::vector<MeshType*>::const_iterator it = meshes.begin(); it != meshes.end(); it++) {
                (*it)->Delete();
            }

            return ActiveShapeModelType::Create(ssm, imagePreprocessor, featureExtractor, profiles);
        }

    private:
        // to prevent (outside) use
        ActiveShapeModelBuilder() {
        }


        ActiveShapeModelBuilder(const ActiveShapeModelBuilder &orig);

        ActiveShapeModelBuilder &operator=(const ActiveShapeModelBuilder &rhs);

    };

};

#endif //STATISMO_ACTIVESHAPEMODELBUILDER_H

