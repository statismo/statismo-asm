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

#ifndef STATISMO_ACTIVESHAPEMODEL_H
#define STATISMO_ACTIVESHAPEMODEL_H

#include "MultiVariateNormalDistribution.h"
#include "ASMOps.h"
#include "StatisticalModel.h"
#include "StatismoIO.h"
#include "ASMProfile.h"
#include "ASMFeatureExtractor.h"
#include "ASMImagePreprocessor.h"
#include "ASMOps.h"


namespace statismo {

    template<typename TPointSet, typename TImage>
    class ActiveShapeModel {

    public:
        typedef Representer<TPointSet> RepresenterType;
        typedef typename ASMOperations<TPointSet>::RigidTransformType RigidTransformationType;
        typedef typename ASMOperations<TPointSet>::RigidTransformPointerType RigidTransformationPointerType;

        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef typename RepresenterType::PointType PointType;
        typedef ASMFeatureExtractorFactory<TPointSet, TImage> FeatureExtractorFactoryType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
        typedef ASMImagePreprocessor<TPointSet, TImage> ImagePreprocessorType;
        typedef ASMImagePreprocessorFactory<TPointSet, TImage> ImagePreprocessorFactoryType;


        static ActiveShapeModel<TPointSet, TImage> *Create(const StatisticalModelType *statisticalModel, const ImagePreprocessorType* preprocessor, const FeatureExtractorType *fe,
                                                           std::vector<ASMProfile> &profiles) {
            return new ActiveShapeModel<TPointSet, TImage>(statisticalModel, preprocessor, fe, profiles);
        };

        virtual ~ActiveShapeModel() {
            if (m_statisticalModel) {
                delete m_statisticalModel;
                m_statisticalModel = 0;
            }
            if (m_featureExtractor) {
                delete m_featureExtractor;
                m_featureExtractor = 0;
            }
            if (m_preprocessor) {
                delete m_preprocessor;
                m_preprocessor = 0;
            }
        }


        std::vector<ASMProfile> GetProfiles() const {
            return m_profiles;
        }

        const StatisticalModelType *GetStatisticalModel() const {
            return m_statisticalModel;
        }

        const RepresenterType* GetRepresenter() const {
            return (const RepresenterType*) GetStatisticalModel()->GetRepresenter();
        }


        const FeatureExtractorType *GetFeatureExtractor() const {
            return m_featureExtractor;
        }

        const ImagePreprocessorType *GetImagePreprocessor() const {
            return m_preprocessor;
        }

        statismo::MultiVariateNormalDistribution GetMarginalAtPointId(unsigned pointId) const {
            PointType imean = m_statisticalModel->DrawMeanAtPoint(pointId);
            MatrixType icovariances = m_statisticalModel->GetCovarianceAtPoint(pointId, pointId);

            unsigned int dimensions = m_statisticalModel->GetRepresenter()->GetDimensions();

            statismo::VectorType mean(dimensions);
            statismo::MatrixType covariances(dimensions, dimensions);
            for (int i = 0; i < dimensions; ++i) {
                mean[i] = imean[i];
                for (int j = 0; j < dimensions; ++j) {
                    covariances(i, j) = icovariances(i, j);
                }
            }

            return statismo::MultiVariateNormalDistribution(mean, covariances);
        }

        ActiveShapeModel<TPointSet, TImage>* CloneWithStatisticalModel(const StatisticalModelType* statisticalModel) const {
            const ImagePreprocessorType* preprocessor = m_preprocessor->Clone();
            const FeatureExtractorType* featureExtractor = m_featureExtractor->Clone();
            std::vector<ASMProfile> profiles = m_profiles;
            ActiveShapeModel<TPointSet, TImage>* r =  new ActiveShapeModel(statisticalModel, preprocessor, featureExtractor, profiles);
            return r;
        }


        // TODO refactor in separate IO method
        void Save(const std::string& filename) {
            // Save the contained SSM
            statismo::IO<TPointSet>::SaveStatisticalModel(GetStatisticalModel(), filename);

            H5::H5File file;

            try {
                file = H5::H5File(filename, H5F_ACC_RDWR);
            } catch (H5::Exception &e) {
                std::string msg(std::string("could not open HDF5 file \n") + e.getCDetailMsg());
                throw StatisticalModelException(msg.c_str());
            }

            H5::Group rootGroup = file.openGroup("/");

            H5::Group asmGroup = rootGroup.createGroup("activeShapeModel");

            HDF5Utils::writeIntAttribute(asmGroup, "majorVersion", 1);
            HDF5Utils::writeIntAttribute(asmGroup, "minorVersion", 0);

            // write feature extractor and image preprocessor metadata
            H5::Group feGroup = asmGroup.createGroup("featureExtractor");
            m_featureExtractor->Save(feGroup);
            feGroup.close();

            H5::Group ppGroup = asmGroup.createGroup("imagePreprocessor");
            m_preprocessor->Save(ppGroup);
            ppGroup.close();

            // write profiles
            H5::Group profilesGroup = asmGroup.createGroup("profiles");

            std::vector<int> pointIds;
            unsigned int profileLength = 0;

            for (std::vector<ASMProfile>::const_iterator it = m_profiles.begin(); it != m_profiles.end(); it++) {
                pointIds.push_back((*it).GetPointId());
                if (profileLength == 0) {
                    profileLength = (*it).GetDistribution().GetCovariance().rows();
                }
            }

            MatrixType means(m_profiles.size(), profileLength);
            MatrixType covariances(m_profiles.size() * profileLength, profileLength);

            unsigned int index = 0;
            for (std::vector<ASMProfile>::const_iterator it = m_profiles.begin(); it != m_profiles.end(); it++) {
                VectorType mean = (*it).GetDistribution().GetMean();
                means.block(index, 0, 1, profileLength) = mean.transpose().block(0, 0, 1, profileLength);

                MatrixType cov = (*it).GetDistribution().GetCovariance();
                covariances.block(index * profileLength, 0, profileLength, profileLength) = cov.block(0, 0, profileLength, profileLength);
                ++index;
            }

            HDF5Utils::writeIntAttribute(profilesGroup, "numberOfPoints", m_profiles.size());
            HDF5Utils::writeIntAttribute(profilesGroup, "profileLength", profileLength);

            HDF5Utils::writeArray(profilesGroup, "pointIds", pointIds);
            HDF5Utils::writeMatrix(profilesGroup, "means", means);
            HDF5Utils::writeMatrix(profilesGroup, "covariances", covariances);

            profilesGroup.close();

            // clean up
            asmGroup.close();
            rootGroup.close();
            file.close();


        }


        static ActiveShapeModel<TPointSet, TImage> *Load(RepresenterType *representer,
                                                               const std::string &filename) {
            H5::H5File file;

            try {
                file = H5::H5File(filename, H5F_ACC_RDONLY);
            } catch (H5::Exception &e) {
                std::string msg(std::string("could not open HDF5 file \n") + e.getCDetailMsg());
                throw StatisticalModelException(msg.c_str());
            }

            H5::Group rootGroup = file.openGroup("/");

            StatisticalModelType *statisticalModel = statismo::IO<TPointSet>::LoadStatisticalModel(representer, rootGroup);

            H5::Group asmGroup = rootGroup.openGroup("activeShapeModel");
            H5::Group feGroup = asmGroup.openGroup("featureExtractor");
            H5::Group ppGroup = asmGroup.openGroup("imagePreprocessor");
            H5::Group profilesGroup = asmGroup.openGroup("profiles");

            std::vector<int> pointIds;
            statismo::MatrixType means;
            statismo::MatrixType covariances;

            unsigned int numPoints = (unsigned int) statismo::HDF5Utils::readIntAttribute(profilesGroup,
                                                                                          "numberOfPoints");
            unsigned int profileLength = (unsigned int) statismo::HDF5Utils::readIntAttribute(profilesGroup,
                                                                                              "profileLength");

            statismo::HDF5Utils::readMatrix(profilesGroup, "covariances", covariances);
            statismo::HDF5Utils::readArray(profilesGroup, "pointIds", pointIds);
            statismo::HDF5Utils::readMatrix(profilesGroup, "means", means);

            std::vector<ASMProfile> profiles;
            profiles.reserve(numPoints);

            for (unsigned int i = 0; i < numPoints; ++i) {
                unsigned int covOffset = i * profileLength;

                statismo::VectorType mean = means.row(i);

                statismo::MatrixType cov(profileLength, profileLength);
                cov.block(0, 0, profileLength, profileLength) = covariances.block(covOffset, 0, profileLength,
                                                                                  profileLength);


                statismo::MultiVariateNormalDistribution mvd(mean, cov);

                profiles.push_back(statismo::ASMProfile(pointIds[i], mvd));
            }

            std::string ppId = HDF5Utils::readStringAttribute(ppGroup, "identifier");
            const ImagePreprocessorFactoryType *ppFactory = ImagePreprocessorFactoryType::GetImplementation(ppId);
            if (!ppFactory) {
                std::string msg(std::string("No image preprocessor implementation found for identifier: ") + ppId);
                throw StatisticalModelException(msg.c_str());
            }
            const ImagePreprocessorType *preprocessor = ppFactory->Instantiate(ppGroup);

            std::string feId = HDF5Utils::readStringAttribute(feGroup, "identifier");
            const FeatureExtractorFactoryType *feFactory = FeatureExtractorFactoryType::GetImplementation(feId);
            if (!feFactory) {
                std::string msg(std::string("No feature extractor implementation found for identifier: ") + feId);
                throw StatisticalModelException(msg.c_str());
            }
            const FeatureExtractorType *featureExtractor = feFactory->Instantiate(feGroup);

            ActiveShapeModel *am = new ActiveShapeModel(statisticalModel, preprocessor, featureExtractor, profiles);

            feGroup.close();
            asmGroup.close();
            rootGroup.close();
            file.close();

            return am;
        }

    protected:
        ActiveShapeModel(const StatisticalModelType *statisticalModel, const ImagePreprocessorType* preprocessor, const FeatureExtractorType *fe,
                         std::vector<ASMProfile> &profiles)
                : m_statisticalModel(statisticalModel),
                  m_preprocessor(preprocessor),
                  m_featureExtractor(fe),
                  m_profiles(profiles) { }

    private:
        std::vector<ASMProfile> m_profiles;
        const StatisticalModelType *m_statisticalModel;
        const ImagePreprocessorType *m_preprocessor;
        const FeatureExtractorType *m_featureExtractor;
    };
};

#endif //STATISMO_ACTIVESHAPEMODEL_H

