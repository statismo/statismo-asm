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
#ifndef STATISMO_MULTIVARIATENORMALDISTRIBUTION_H
#define STATISMO_MULTIVARIATENORMALDISTRIBUTION_H

#include "CommonTypes.h"
#include <cmath>
#include "StatismoUtils.h"

namespace statismo {
    class MultiVariateNormalDistribution {

    private:
        MatrixType covInv;

    public:

        typedef std::vector<VectorType> DataType;

        VectorType mean;
        MatrixType covariance;


        MultiVariateNormalDistribution() {}

        MultiVariateNormalDistribution(VectorType _mean, MatrixType _covariance): mean(_mean), covariance(_covariance) {
            if (!pseudoInverse(covariance, covInv)) {
                //FIXME: throw some exception
            }
        }

        float logpdf(VectorType data) const {
            double mhDist = MahalanobisDistance(data);
            double detCov = 1;
            for (unsigned i = 0; i < covariance.rows(); ++i) {
                detCov *= covariance(i,i);
            }
            double normalizer = -0.5*mean.size() * std::log(2*M_PI) -0.5 * detCov;
            return -0.5 * mhDist + normalizer;
        }

        float MahalanobisDistance(VectorType data) const {
            VectorType x0 = data - mean;

            float d = sqrt(x0.dot((covInv * x0)));
            return d;
        }


        const VectorType& GetMean() const {
            return mean;
        }

        const MatrixType& GetCovariance() const {
            return covariance;
        }

        /**
            * @brief Calculate the pseudo-inverse of a matrix.
            * @param a the matrix to calculate the pseudo-inverse of.
            * @param result OUT parameter where the result is stored.
            * @return a boolean indicating whether the operation was successful.
        */
        static bool pseudoInverse(const MatrixType &a, MatrixType &result, double epsilon = std::numeric_limits<typename MatrixType::Scalar>::epsilon())
        {
            if(a.rows()<a.cols()) {
                return false;
            }

            Eigen::JacobiSVD< MatrixType > svd = a.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

            MatrixType::Scalar tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs().maxCoeff();

            result = svd.matrixV() * MatrixType( (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().
                    array().inverse(), 0) ).asDiagonal() * svd.matrixU().adjoint();
            return true;
        }


        static MultiVariateNormalDistribution EstimateFromData(const DataType& data) {
            unsigned int numSamples = data.size();

            VectorType mean = data[0];
            for (unsigned int i= 1; i < numSamples; ++i) {
                mean += data[i];
            }
            mean /= numSamples;

            MatrixType covariance(mean.rows(), mean.rows());
            for (int r = covariance.rows() - 1; r >= 0; --r) {
                for (int c =covariance.cols() - 1; c >= 0; --c ) {
                    covariance(r,c) = 0;
                }
            }

            for (DataType::const_iterator item = data.begin(); item != data.end(); item++) {
                VectorType diff = (*item) - mean;
                MatrixType outer = diff * diff.transpose();
                covariance = covariance + outer;
            }

            if (numSamples > 1) {
                covariance /= (numSamples - 1);
            }
            return MultiVariateNormalDistribution(mean, covariance);
        }

    };
}
#endif //STATISMO_MULTIVARIATENORMALDISTRIBUTION_H
