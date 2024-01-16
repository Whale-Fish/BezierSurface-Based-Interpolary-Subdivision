#pragma once
#include<vector>
#include "Kernel.h"
using namespace std;

class subdiv
{
public:
	subdiv();
	~subdiv();

	vector<PolygonMesh::VertexHandle> vhandle;
	vector<T_PolygonMesh::VertexHandle> Tvhandle;

	int getindex(PolygonMesh::VertexHandle vh);
	void create_new_point_and_nor(POINT3D P1,POINT3D P2,POINT3D N1,POINT3D N2,POINT3D &newP,POINT3D &newN);

	void subdivsion(vector<POINT3D> points, vector<POINT3D> nors, vector<POINT3D>&resp, vector<POINT3D>&resn);
	//void sub_for_mesh(PolygonMesh* mesh, PolygonMesh &subvMesh);
	void sub_for_mesh_new(PolygonMesh* mesh, PolygonMesh &subvMesh);
	void sub_for_quad_mesh_new(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh);

	void sub_for_mesh_by_sqrt3(PolygonMesh* mesh, PolygonMesh &subvMesh);
	void sub_for_mesh_by_sqrt3_PN_method(PolygonMesh* mesh, PolygonMesh &subvMesh);

	void sub_for_mesh_by_sqrt2(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh);

	void sqrt3(PolygonMesh* mesh, PolygonMesh &subvMesh);
	void sqrt2(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh);
	void Isqrt3(PolygonMesh* mesh, PolygonMesh &subvMesh);
	void Isqrt2(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh);
	//void sub_for_mesh_by_sqrt3_test(PolygonMesh* mesh, PolygonMesh &subvMesh);
};

