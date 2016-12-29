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

#ifndef STATISMO_ASMFITTING_H
#define STATISMO_ASMFITTING_H

#include "ActiveShapeModel.h"
#include "ASMPointSampler.h"
#include "ASMOps.h"
#include <boost/thread.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread/future.hpp>

namespace statismo {

    class ASMFittingConfiguration {
    private:
        float m_featureDistanceThreshold;
        float m_pointDistanceThreshold;
        float m_modelCoefficientBounds;

    public:
        ASMFittingConfiguration(float featureDistanceThreshold, float pointDistanceThreshold,
                                float modelCoefficientBounds) :
                m_featureDistanceThreshold(featureDistanceThreshold),
                m_pointDistanceThreshold(pointDistanceThreshold),
                m_modelCoefficientBounds(modelCoefficientBounds) { }


        float GetFeatureDistanceThreshold() const {
            return m_featureDistanceThreshold;
        }

        float GetPointDistanceThreshold() const {
            return m_pointDistanceThreshold;
        }

        float GetModelCoefficientBounds() const {
            return m_modelCoefficientBounds;
        }
    };

    template<typename RigidTransformPointerType>
    class ASMFittingResult {
    public:
        ASMFittingResult() : m_isValid(false) {
        }

        ASMFittingResult(bool isValid, VectorType coefficients, RigidTransformPointerType rigidTransform) :
                m_isValid(isValid),
                m_coefficients(coefficients),
                m_rigidTransform(rigidTransform) { }

        bool IsValid() {
            return m_isValid;
        }

        VectorType GetCoefficients() {
            return m_coefficients;
        }

        RigidTransformPointerType GetRigidTransform() {
            return m_rigidTransform;
        }

    private:
        bool m_isValid;
        VectorType m_coefficients;
        RigidTransformPointerType m_rigidTransform;
    };

    template<typename TPointSet, typename TImage>
    class ASMFittingStep {
        typedef ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef ASMPointSampler<TPointSet, TImage> PointSamplerType;
        typedef typename Representer<TPointSet>::PointType PointType;
        typedef StatisticalModel<TPointSet> StatisticalModelType;
        typedef ASMFeatureExtractor<TPointSet, TImage> FeatureExtractorType;
        typedef ASMPreprocessedImage<TPointSet> PreprocessedImageType;
        typedef typename ASMOperations<TPointSet>::RigidTransformPointerType RigidTransformPointerType;


        ASMFittingStep(const ASMOperations<TPointSet>* asmOperations, const ASMFittingConfiguration &configuration, const ActiveShapeModelType *activeShapeModel,
                       const VectorType sourceCoefficients, const RigidTransformPointerType sourceTransform,
                       const PreprocessedImageType *targetImage, const PointSamplerType *sampler)
                :
                m_asmOperations(asmOperations),
                m_configuration(configuration),
                m_model(activeShapeModel),
                m_sourceCoefficients(sourceCoefficients),
                m_sourceTransform(sourceTransform),
                m_target(targetImage),
                m_sampler(sampler) { }

    private:




        class ProfileResult {
        public:
            unsigned int pointId;
            PointType candidatePoint;
            PointType transformedCandidatePoint;
        };

        class ChunkResult {
        public:
            ChunkResult() { }

            std::vector<ProfileResult> results;

            // emulate move semantics, as boost::async seems to depend on it.
            ChunkResult &operator=(BOOST_COPY_ASSIGN_REF(ChunkResult)rhs) { // Copy assignment
                if (&rhs != this) {
                    copyMembers(rhs);
                }
                return *this;
            }

            ChunkResult(BOOST_RV_REF(ChunkResult)that) { //Move constructor
                copyMembers(that);
            }

            ChunkResult &operator=(BOOST_RV_REF(ChunkResult)rhs) { //Move assignment
                if (&rhs != this) {
                    copyMembers(rhs);
                }
                return *this;
            }

        private:
        BOOST_COPYABLE_AND_MOVABLE(ChunkResult)

            void copyMembers(const ChunkResult &that) {
                results = that.results;
            }
        };

        class ChunkRequest {
        public:
            ChunkRequest(std::vector<ASMProfile> *_profiles, unsigned int _dimensions) : profiles(_profiles), dimensions(_dimensions) { }

            unsigned int startInclusive;
            unsigned int endExclusive;
            std::vector<ASMProfile> *profiles;
            unsigned int dimensions;
        };


        ChunkResult PerformChunk(ChunkRequest in) const {
            ChunkResult out;
            unsigned int index = in.startInclusive;


            typename statismo::ASMOperations<TPointSet>::MeshPointerType transformedMesh = m_asmOperations->TransformMesh(m_model->GetStatisticalModel()->DrawSample(m_sourceCoefficients), m_sourceTransform);

            PointSamplerType *sampler = m_sampler->CloneForTarget(m_model, m_sourceCoefficients, transformedMesh);
            FeatureExtractorType *fe = m_model->GetFeatureExtractor()->CloneForTarget(m_model, m_sourceCoefficients, transformedMesh);

            for (std::vector<ASMProfile>::const_iterator profile = (in.profiles->begin() + index);
                 index < in.endExclusive; ++index, ++profile) {
                ProfileResult result;
                unsigned pointId = (*profile).GetPointId();
                PointType profilePoint = m_model->GetStatisticalModel()->DrawSampleAtPoint(m_sourceCoefficients, pointId, false);
                PointType transformedProfilePoint = m_asmOperations->TransformPoint(profilePoint, m_sourceTransform, false);
                PointType transformedCandidatePoint;
                float featureDistance = FindBestMatchingPointForProfile(transformedCandidatePoint, sampler, fe,
                                                                        transformedProfilePoint,
                                                                        (*profile).GetDistribution());
                //std::cout << profilePoint << " -> " << transformedProfilePoint << " " << featureDistance << std::endl;
                if (featureDistance >= 0 && featureDistance <= m_configuration.GetFeatureDistanceThreshold()) {
                    PointType candidatePoint =m_asmOperations->TransformPoint(transformedCandidatePoint, m_sourceTransform, true);
                    statismo::VectorType point(in.dimensions);
                    for (int i = 0; i < in.dimensions; ++i) {
                        point[i] = candidatePoint[i];
                    }
                    statismo::MultiVariateNormalDistribution marginal = m_model->GetMarginalAtPointId(pointId);
                    float pointDistance = marginal.MahalanobisDistance(point);

                    if (pointDistance <= m_configuration.GetPointDistanceThreshold()) {
                        result.pointId = pointId;
                        result.candidatePoint = candidatePoint;
                        result.transformedCandidatePoint = transformedCandidatePoint;
                        out.results.push_back(result);
                    } else {
                        // shape distance is larger than threshold
                    }
                } else {
                    // feature distance is larger than threshold
                }
            }
            sampler->Delete();
            fe->Delete();
            return out;
        }

    public:


        static ASMFittingStep *Create(const ASMOperations<TPointSet>* asmOperations,
                                        const ASMFittingConfiguration &configuration,
                                      const ActiveShapeModelType *activeShapeModel,
                                      const VectorType sourceCoefficients,
                                      const RigidTransformPointerType sourceTransform,
                                      const PreprocessedImageType *targetImage, const PointSamplerType *sampler) {
            return new ASMFittingStep(asmOperations, configuration, activeShapeModel, sourceCoefficients, sourceTransform, targetImage, sampler);
        }

        ~ASMFittingStep() {
        }

        ASMFittingResult<RigidTransformPointerType> Perform() const {
            typename StatisticalModelType::PointValueListType deformations;
            std::vector<PointType> referencePoints;
            std::vector<PointType> targetPoints;

            std::vector<ASMProfile> profiles = m_model->GetProfiles();

            unsigned int dimensions = m_model->GetStatisticalModel()->GetRepresenter()->GetDimensions();

            std::vector<PointType> domainPoints = m_model->GetStatisticalModel()->GetRepresenter()->GetDomain().GetDomainPoints();

            unsigned int profilesCount = profiles.size();
            unsigned int chunksCount = boost::thread::hardware_concurrency();
            if (chunksCount > profilesCount) chunksCount = profilesCount; // anybody have a 2048-core machine? ;-)
            if (chunksCount == 0) chunksCount = 1; // just in case
            unsigned int chunkSize = profilesCount / chunksCount;

            std::vector<ChunkRequest> requests;
            std::vector<boost::future<ChunkResult> *> results;

            for (unsigned int i = 0; i < chunksCount; ++i) {
                ChunkRequest request(&profiles, dimensions);
                request.startInclusive = i * chunkSize;
                if (i == chunksCount - 1) {
                    request.endExclusive = profilesCount;
                } else {
                    request.endExclusive = (i + 1) * chunkSize;
                }
                requests.push_back(request);
            }

            for (typename std::vector<ChunkRequest>::const_iterator r = requests.begin(); r != requests.end(); ++r) {
                ChunkRequest req = (*r);
                boost::future<ChunkResult> *fut = new boost::future<ChunkResult>(boost::async(boost::launch::async,boost::bind(&ASMFittingStep<TPointSet, TImage>::PerformChunk,this, req)));
                results.push_back(fut);
            }

            for (typename std::vector<boost::future<ChunkResult> *>::const_iterator future = results.begin();
                 future != results.end(); ++future) {
                ChunkResult chunk = (*future)->get();
                for (typename std::vector<ProfileResult>::const_iterator p = chunk.results.begin();
                     p != chunk.results.end(); ++p) {
                    PointType refPoint = domainPoints[(*p).pointId];
                    deformations.push_back(std::make_pair(refPoint, (*p).candidatePoint));
                    referencePoints.push_back(refPoint);
                    targetPoints.push_back((*p).transformedCandidatePoint);
                }
                delete (*future);
            }

            if (deformations.size() > 0) {
                //std::cout << "constraints: " << deformations.size() << std::endl;
                VectorType coeffs = m_model->GetStatisticalModel()->ComputeCoefficientsForPointValues(deformations, 1e-6);

                float maxCoeff = m_configuration.GetModelCoefficientBounds();
                for (size_t i = 0; i < coeffs.size(); ++i) {
                    coeffs(i) = std::min(maxCoeff, std::max(-maxCoeff, coeffs(i)));
                }
                RigidTransformPointerType adjust = m_asmOperations->ComputeRigidTransformFromLandmarks( referencePoints, targetPoints);
                return ASMFittingResult<RigidTransformPointerType>(true, coeffs, adjust);

            } else {
                // Invalid result: coefficients with no content
                //std::cout << "returning invalid result" << std::endl;
                VectorType coeffs;
                return ASMFittingResult<RigidTransformPointerType>(false, coeffs, 0);
            }
        }

    private:
        const ASMOperations<TPointSet>* m_asmOperations;
        const ActiveShapeModelType *m_model;
        const VectorType m_sourceCoefficients;
        const RigidTransformPointerType m_sourceTransform;
        const PreprocessedImageType *m_target;
        const PointSamplerType *m_sampler;
        ASMFittingConfiguration m_configuration;

        float FindBestMatchingPointForProfile(PointType &result, const PointSamplerType *sampler,
                                              const FeatureExtractorType *const fe, const PointType &profilePoint,
                                              const statismo::MultiVariateNormalDistribution &profile) const {
            float bestFeatureDistance = -1;

            std::vector<PointType> samples = sampler->Sample(profilePoint);

            for (typename std::vector<PointType>::const_iterator sample = samples.begin();
                 sample != samples.end(); ++sample) {
                statismo::VectorType features;
                bool ok = fe->ExtractFeatures(features, m_target, *sample);
                if (ok) {
                    float featureDistance = profile.MahalanobisDistance(features);
                    if (bestFeatureDistance == -1 || bestFeatureDistance > featureDistance) {
                        result = *sample;
                        bestFeatureDistance = featureDistance;
                    }
                }
            }
            return bestFeatureDistance;
        }
    };
}

#endif //STATISMO_ASMFITTING_H

