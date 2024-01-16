#pragma once

#ifndef KERNEL_H
#define KERNEL_H

#include<vector>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/io/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::PolyMesh_ArrayKernelT<> T_PolygonMesh;  //quad mesh
typedef OpenMesh::TriMesh_ArrayKernelT<> PolygonMesh;
typedef OpenMesh::Vec3d POINT3D;



class NOR_POINT
{
	POINT3D pos;
	POINT3D nor;
};

class SURFACE
{
public:
	PolygonMesh mesh;
	PolygonMesh nor;
};



#endif