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


#ifndef STATISMO_ITKDEFERREDITEMS_H
#define STATISMO_ITKDEFERREDITEMS_H

#include "DeferredItem.h"
#include <itkImageFileReader.h>
#include <itkMeshFileReader.h>



namespace itk {
    template<typename ImageType>
    class DeferredImage: public statismo::DeferredItem<ImageType*> {
    public:
        typedef itk::ImageFileReader<ImageType> ImageReaderType;

        DeferredImage(std::string filename) : m_filename(filename) {}

        virtual ImageType* Get() {
            typename ImageReaderType::Pointer reader = ImageReaderType::New();
            reader->SetFileName(m_filename.c_str());
            reader->Update();
            typename ImageType::Pointer image = reader->GetOutput();
            // this is needed because we're returning a raw pointer
            image->Register();
            return image;
        }

    private:
        std::string m_filename;
    };


    template<typename MeshType>
    class DeferredMesh: public statismo::DeferredItem<MeshType*> {
    public:
        typedef itk::MeshFileReader<MeshType> MeshReaderType;

        DeferredMesh(std::string filename) : m_filename(filename) {}

        virtual MeshType* Get() {
            typename MeshReaderType::Pointer reader = MeshReaderType::New();
            reader->SetFileName(m_filename.c_str());
            reader->Update();
            typename MeshType::Pointer mesh = reader->GetOutput();
            // this is needed because we're returning a raw pointer
            mesh->Register();
            return mesh;
        }

    private:
        std::string m_filename;
    };

//  This was an attempt to unify things, but it gives crazy warnings when compiling with GCC.
//  "warning: ‘itk::DeferredMesh<itk::Mesh<float, 3u> >’ declared with greater visibility than the type of its field ‘itk::DeferredMesh<itk::Mesh<float, 3u> >::<anonymous>’ [-Wattributes]"

//    template<typename DataType, typename ReaderType>
//    class DeferredLoad: public statismo::DeferredItem<DataType*> {
//    public:
//        DeferredLoad(std::string filename) : m_filename(filename) {}
//
//        virtual DataType* Get() {
//            typename ReaderType::Pointer reader = ReaderType::New();
//            reader->SetFileName(m_filename.c_str());
//            reader->Update();
//            typename DataType::Pointer data = reader->GetOutput();
//            // this is needed because we're returning a raw pointer
//            data->Register();
//            return data;
//        }
//
//    private:
//        std::string m_filename;
//    };
//
//
//    template<typename MeshType>
//    class DeferredMesh: public itk::DeferredLoad<MeshType, itk::MeshFileReader<MeshType> > {
//    public:
//        typedef itk::MeshFileReader<MeshType> MeshReaderType;
//        DeferredMesh(std::string filename) : DeferredLoad<MeshType, MeshReaderType>(filename) {}
//    };
}

#endif //STATISMO_ITKDEFERREDITEMS_H

