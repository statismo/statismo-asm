/*
 * This file is part of the statismo library.
 *
 * Author: Anita Lerch
 *
 * Copyright (c) 2008-2015 University of Basel
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

#ifndef __ITKTRIANGLEMESH_H_
#define __ITKTRIANGLEMESH_H_

#include <itkObject.h>
#include <itkTriangleCell.h>
#include <itkCovariantVector.h>
#include <itkMesh.h>
#include <vector>
#include <set>

namespace itk {

template <
typename TPixelType,
unsigned int VDimension = 3,
typename TMeshTraits = DefaultStaticMeshTraits< TPixelType , VDimension, VDimension >
>
/*! The TriangleMeshAdapter class. */
/*!
 * The TriangleMeshAdapter extends the mesh class, with specific function for the triangle mesh.
 */
class ITK_EXPORT TriangleMeshAdapter : public Object 
{
public:

  /* Standard typedefs. */
  typedef TriangleMeshAdapter  Self;
  typedef itk::Mesh<TPixelType, VDimension, TMeshTraits>  Mesh;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TriangleMeshAdapter, Object );

  /* Hold on to the type information specified by the template parameters. */
  typedef typename Mesh::PixelType     PixelType;
  typedef typename Mesh::CellPixelType CellPixelType;

  /* Convenient typedefs obtained from TMeshTraits template parameter. */
  typedef typename Mesh::Pointer         MeshPointer;
  typedef typename Mesh::PointIdentifier         PointIdentifier;
  typedef typename Mesh::CellIdentifier          CellIdentifier;
  typedef typename Mesh::PointType               PointType;
  typedef typename Mesh::PointsContainer         PointsContainer;
  typedef typename Mesh::CellsContainer          CellsContainer;
  typedef typename Mesh::CellLinksContainer      CellLinksContainer;
  typedef typename Mesh::PointDataContainer      PointDataContainer;
  typedef typename Mesh::CellDataContainer       CellDataContainer;
  typedef typename Mesh::PointCellLinksContainer PointCellLinksContainer;

  /* Create types that are pointers to each of the container types. */
  typedef typename PointsContainer::Pointer       PointsContainerPointer;
  typedef typename CellsContainer::Pointer        CellsContainerPointer;
  typedef typename CellLinksContainer::Pointer    CellLinksContainerPointer;
  typedef typename PointDataContainer::Pointer    PointDataContainerPointer;
  typedef typename CellDataContainer::Pointer     CellDataContainerPointer;

  /* Create types for the TriangleMeshAdapter. */
  typedef typename itk::CovariantVector    <PixelType,VDimension>  PointNormalType;
  typedef typename itk::VectorContainer< PointIdentifier,PointNormalType >     PointNormalsContainer;
  typedef typename PointNormalsContainer::Pointer     PointNormalsContainerPointer;
  typedef typename itk::CovariantVector<CellPixelType, VDimension>  CellNormalType;
  typedef typename itk::VectorContainer< CellIdentifier,CellNormalType >     CellNormalsContainer;
  typedef typename CellNormalsContainer::Pointer     CellNormalsContainerPointer;

  /* Create types that are iterators for each of the container types. */
  typedef typename
  PointsContainer::ConstIterator        PointsContainerConstIterator;
  typedef typename
  PointsContainer::Iterator             PointsContainerIterator;
  typedef typename
  CellsContainer::ConstIterator         CellsContainerConstIterator;
  typedef typename
  CellsContainer::Iterator              CellsContainerIterator;
  typedef typename
  CellLinksContainer::ConstIterator     CellLinksContainerIterator;
  typedef typename
  PointDataContainer::ConstIterator     PointDataContainerIterator;
  typedef typename
  CellDataContainer::ConstIterator      CellDataContainerIterator;
  typedef typename 
  PointCellLinksContainer::iterator    PointCellLinksContainerIterator;

  /* The base cell type for cells in this mesh. */
  typedef typename Mesh::CellType             MeshCellType;
  typedef typename itk::TriangleCell< MeshCellType >    TriangleCellType;


  /* Triangle Mesh Adapter Method. */
  void ComputeCellAreas();
  CellDataContainerPointer GetCellAreas();
  void ComputePointAreas();
  PointDataContainerPointer GetPointAreas();
  void ComputeCellNormals();
  CellNormalsContainerPointer GetCellNormals();
  void ComputePointNormals();
  PointNormalsContainerPointer GetPointNormals();
  PointType GetCenterOfGravity();
  bool SetMesh(Mesh *mesh);
  void PrintPointData();


protected:
  TriangleMeshAdapter();
  ~TriangleMeshAdapter();


private:

  CellDataContainerPointer m_cell_areas;
    PointDataContainerPointer m_point_areas;
  CellNormalsContainerPointer m_cell_normals;
  PointNormalsContainerPointer m_point_normals;
  MeshPointer m_mesh;

  TriangleMeshAdapter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTriangleMeshAdapter.hxx"
#endif
#endif /*TRIANGLEMESH_H_*/
