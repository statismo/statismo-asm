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

#ifndef STATISMO_ITKASMRANDOMPOINTIDSAMPLER_H
#define STATISMO_ITKASMRANDOMPOINTIDSAMPLER_H

#include <itkObject.h>
#include <itkObjectFactory.h>

#include "ASMPointIdSampler.h"

namespace itk {
    template <typename MeshType>
    class ASMRandomPointIdSampler : public statismo::ASMPointIdSampler, public Object {
    public:
        /* standard ITK typedefs and macros */
        typedef ASMRandomPointIdSampler Self;
        typedef Object Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(Self, Object);

        ASMRandomPointIdSampler(): m_mesh(0), m_numberOfPoints(0) {}

        void SetNumberOfPoints(unsigned int numberOfPoints) {
            m_numberOfPoints = numberOfPoints;
        }

        void SetMesh(typename MeshType::Pointer mesh) {
            m_mesh = mesh;
        }

        virtual std::vector<unsigned int> SamplePoints() const {
            assert(m_numberOfPoints > 0 && m_pointIds.size() == m_numberOfPoints);
            return m_pointIds;
        }

        void Update() {
            assert (m_numberOfPoints > 0 && m_mesh);

            if (m_pointIds.size() == 0) {
                unsigned int maximum = m_mesh->GetNumberOfPoints();
                if (maximum < m_numberOfPoints) {
                    std::cout << "Warning: requested number of points (" << m_numberOfPoints <<
                              ") is larger than number of points on mesh (" << maximum << "), capping." << std::endl;
                    m_numberOfPoints = maximum;
                }

                std::vector<bool> selected(maximum);

                for (unsigned int i = 0; i < maximum; ++i) {
                    selected[i] = false;
                }

                unsigned int targetCount = m_numberOfPoints;
                bool keep_filter = true;

                // if we're supposed to sample more than half of the points, we invert the logic,
                // i.e., we're going to randomly draw points which will *not* be used.
                if (targetCount > maximum / 2) {
                    targetCount = maximum - targetCount;
                    keep_filter = false;
                }

                boost::random::mt19937 prng;
                prng.seed(std::time(0));
                boost::random::uniform_int_distribution<> draw(0, maximum - 1);

                for (unsigned int count = 0; count < targetCount;) {
                    unsigned int index = draw(prng);
                    if (!selected[index]) {
                        selected[index] = true;
                        ++count;
                    }
                };

                for (unsigned int i = 0; i < maximum; ++i) {
                    if (selected[i] == keep_filter) {
                        m_pointIds.push_back(i);
                    }
                }
            }
        }

    private:

        std::vector<unsigned int> m_pointIds;
        typename MeshType::Pointer m_mesh;
        unsigned int m_numberOfPoints;
    };

}

#endif //STATISMO_ITKASMRANDOMPOINTIDSAMPLER_H

