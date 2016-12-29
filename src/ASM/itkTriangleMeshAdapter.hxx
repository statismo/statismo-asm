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

#ifndef __itkTriangleMeshAdapter_hxx
#define __itkTriangleMeshAdapter_hxx

#include "itkTriangleMeshAdapter.h"


namespace itk { 
/**
* A protected default constructor allows the New() routine to create an
* instance of Mesh.  All the containers are initialized to empty
* containers.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::TriangleMeshAdapter():
  m_cell_areas(CellDataContainer::New()),
  m_point_areas(PointDataContainer::New()),
  m_cell_normals(CellNormalsContainer::New()),
  m_point_normals(PointNormalsContainer::New()),
  m_mesh(0)
  {
  }

/**
* TriangleMeshAdapter Destructor takes care of releasing the memory of Cells
* and CellBoundaries objects for which normal pointers are stored.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::~TriangleMeshAdapter()
{
}

/*!
  Set the mesh to the TriangleMeshAdapter.
  \param mesh The triangle mesh, on which the triangle mesh adapter operated on.
  \return True if the triangle mesh was set.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
bool 
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::SetMesh(Mesh *mesh)
{
  CellsContainerIterator cellIter = mesh->GetCells()->Begin();  
  CellsContainerIterator end = mesh->GetCells()->End();
  while( cellIter != end )
  {
    MeshCellType * cellptr = cellIter.Value();
    if(dynamic_cast<TriangleCellType *>( cellptr ) == 0)        
    {      
      itkDebugMacro(<< "itk::TriangleMeshAdapter::SetMesh() the input mesh" 
          << typeid(mesh).name()
          << " is not a triangle Mesh.");
      itkExceptionMacro(<< "itk::TriangleMeshAdapter::SetMesh() the input mesh" 
          << typeid(mesh).name()
          << " is not a triangle Mesh.");      
      return false;
    }
    ++cellIter;
  }

  this->m_mesh = mesh;
  return true; // return true only if all mesh cells are triangle cells
}

/*!
  Compute the area of the triangel cells.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
void 
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::ComputeCellAreas()
{
  itkDebugMacro(<< "itk::TriangleMeshAdapter::ComputeCellAreas()");

  this->ComputeCellNormals();
  CellsContainerIterator  cellIterator = this->m_mesh->GetCells()->Begin(); 

  CellIdentifier i=0;
  while( cellIterator != this->m_mesh->GetCells()->End() ) 
  {
    const CellIdentifier *pointIds = cellIterator.Value()->GetPointIds();
    PointType p1, p2, p3;
    this->m_mesh->GetPoint (pointIds[0], &p1);
    this->m_mesh->GetPoint (pointIds[1], &p2);
    this->m_mesh->GetPoint (pointIds[2], &p3);

    CellNormalType va, vb;
    va = p2-p1;
    vb = p3-p1;
    CellPixelType vp1 = va[1]*vb[2]-va[2]*vb[1];
    CellPixelType vp2 = va[2]*vb[0]-va[0]*vb[2];
    CellPixelType vp3 = va[0]*vb[1]-va[1]*vb[0];

    CellNormalType vp;
    vp[0] = vp1;
    vp[1] = vp2;
    vp[2] = vp3;
    this->m_cell_normals->InsertElement (i, vp);

    CellPixelType a = sqrt ((vp1*vp1)+(vp2*vp2)+(vp3*vp3))/2;  

    this->m_cell_areas->InsertElement (i, a);
    ++cellIterator;  
    ++i;
  }
}

/*! 
  Get the area of the triangel cells.
  \return A CellDataContainerPointer, which contain the areas of the cells.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
typename TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>::CellDataContainerPointer
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::GetCellAreas()
{
  itkDebugMacro(<< "itk::TriangleMeshAdapter::GetCellAreas()");
  this->ComputeCellAreas();
  return this->m_cell_areas;
}

/*!
  Compute the corresponding area of the mesh points.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
void 
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::ComputePointAreas()
{
  this->m_mesh->BuildCellLinks();
  this->ComputeCellAreas();

  PointsContainerIterator pointIterator = this->m_mesh->GetPoints()->Begin(); 
  CellLinksContainerPointer cellLinksContainer = this->m_mesh->GetCellLinks(); 

  PointIdentifier j=0;
  while( pointIterator != this->m_mesh->GetPoints()->End() ) 
  {
    PointType point = pointIterator.Value();
    PointCellLinksContainer ptCellLink = cellLinksContainer->GetElement(j);
    //PointCellLinksContainer ptCellLink = (*cellLinksContainer)[j];
    PixelType sum = 0;
    PointCellLinksContainerIterator it;
    for (it = ptCellLink.begin(); it != ptCellLink.end(); it++) { 
      sum += this->m_cell_areas->GetElement(*it)/3;
    }  
    this->m_point_areas->InsertElement (j, sum);

    ++j;
    ++pointIterator;  
  }
}

/*! 
  Get the corresponding area of the mesh points.
  \return A PointDataContainerPointer, which contain the corresponding area of the mesh points.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
typename TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>::PointDataContainerPointer
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::GetPointAreas()
{
  this->ComputePointAreas();  
  return this->m_point_areas;
}

/*!
  Compute the normal of the triangl cells.
  Note that the results are NOT normalized to unit vectors.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
void 
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::ComputeCellNormals()
{
  itkDebugMacro(<< "itk::TriangleMeshAdapter::ComputeCellNormals()");

  CellsContainerIterator  cellIterator = this->m_mesh->GetCells()->Begin(); 

  CellIdentifier i=0;
  while( cellIterator != this->m_mesh->GetCells()->End() ) 
  {
    const CellIdentifier *pointIds = cellIterator.Value()->GetPointIds();
    PointType p1, p2, p3;
    this->m_mesh->GetPoint (pointIds[0], &p1);
    this->m_mesh->GetPoint (pointIds[1], &p2);
    this->m_mesh->GetPoint (pointIds[2], &p3);

    CellNormalType va, vb, vp;
    va = p2-p1;
    vb = p3-p1;
    vp[0] = va[1]*vb[2]-va[2]*vb[1];
    vp[1] = va[2]*vb[0]-va[0]*vb[2];
    vp[2] = va[0]*vb[1]-va[1]*vb[0];

    this->m_cell_normals->InsertElement (i, vp);
    ++cellIterator;  
    ++i;
  }
}

/*! 
  Get the normal of the triangle cells.
  \return A CellNormalsContainerPointer, which contain the normals of the triangle cells.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
typename TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>::CellNormalsContainerPointer
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::GetCellNormals()
{
  this->ComputeCellNormals();
  return this->m_cell_normals;
}

/*!
  Compute the normal of the mesh point.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
void 
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::ComputePointNormals()
{
  this->m_mesh->BuildCellLinks();
  this->ComputeCellNormals();

  PointsContainerIterator pointIterator = this->m_mesh->GetPoints()->Begin(); 
  PointsContainerIterator pointIteratorEnd = this->m_mesh->GetPoints()->End();
  CellLinksContainerPointer cellLinksContainer = this->m_mesh->GetCellLinks(); 

  PointIdentifier j=0;
  while( pointIterator != pointIteratorEnd )
  {
    PointType point = pointIterator.Value();
    PointCellLinksContainer ptCellLink = cellLinksContainer->GetElement(j);
    PointNormalType sum;
    sum.Fill(0.0);
    PointCellLinksContainerIterator it;
    for (it = ptCellLink.begin(); it != ptCellLink.end(); it++) {
      sum += this->m_cell_normals->GetElement(*it);
    }  
    this->m_point_normals->InsertElement (j, sum);

    ++j;
    ++pointIterator;  
  }

}

/*! 
  Get the normal of the mesh points.
  \return A PointNormalsContainerPointer, which contain the normals the mesh points.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
typename TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>::PointNormalsContainerPointer
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::GetPointNormals()
{
  this->ComputePointNormals();
  return this->m_point_normals;
}

/*!
  Compute the center of gravity.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
typename TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>::PointType
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::GetCenterOfGravity()
{
  this->ComputePointAreas();

  PointsContainerIterator pointIterator = this->m_mesh->GetPoints()->Begin();
  PointsContainerIterator pointIteratorEnd = this->m_mesh->GetPoints()->End(); 
  PointDataContainerIterator pointDataIterator = this->m_point_areas->Begin(); 

  PointType sum;
  sum.Fill(0);
  PixelType area = 0;
  while( pointIterator != pointIteratorEnd ) 
  {
    PointType p = pointIterator.Value();
    PixelType w = pointDataIterator.Value();
    sum[0] += p[0]*w;
    sum[1] += p[1]*w;
    sum[2] += p[2]*w;
    area += w;
    ++pointDataIterator;
    ++pointIterator;
  }

  sum[0] = sum[0]/area;
  sum[1] = sum[1]/area;
  sum[2] = sum[2]/area;
  return sum;
}

/*!
  Print information of the triangle mesh adapter.
*/
template <typename TPixelType, unsigned int VDimension, typename TMeshTraits>
void 
TriangleMeshAdapter<TPixelType, VDimension, TMeshTraits>
::PrintPointData()
{
  std::cout << "Adapter Mesh Container size: " << m_point_areas->Size() << std::endl;
  std::cout << "Adapter Mesh Container address: " << m_point_areas.GetPointer() << std::endl;  
  std::cout << "Adapter Mesh Point Data - Adress: " << this->m_mesh->GetPointData() 
  << "\nAdapter Mesh Point Data - Size: " << this->m_mesh->GetPointData()->Size() << std::endl;
  std::cout << "Container size: " << m_point_areas->Size() << std::endl;
}


} // namespace itk


#endif //#ifndef __itkTriangleMeshAdapter_hxx
