#include "subdiv.h"
#include<math.h>
#include<cstdlib>
#include<set>
#include<map>
#include<string.h>
#include<cmath>
#include "eigenSolution.h"

#define PI 3.1415926525

double dot(POINT3D a, POINT3D b);
POINT3D cross(POINT3D a, POINT3D b);
PolygonMesh tmpMesh;

double two_bezier_triangle_cos = 1.000;
int error_normal_cnt = 0;

bool testFlg = false;

int cnttt = 0;

typedef struct from_to_id_pair
{
	int id1;
	int id2;
	bool operator==(const from_to_id_pair a) const
	{
		return ((this->id1 == a.id1) && (this->id2 == a.id2));
	}
	bool operator<(const from_to_id_pair a) const
	{
		if (this->id1 < a.id1)
		{
			return true;
		}
		else if (this->id1 == a.id1)
		{
			return this->id2 < a.id2;
		}
		else
		{
			return false;
		}	
	}

}id_pair;

double dot(POINT3D a, POINT3D b)
{
	return sqrt(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

POINT3D cross(POINT3D a, POINT3D b)
{
	POINT3D res;
	res[0] = a[1]*b[2] - a[2]*b[1];
	res[1] = -a[0]*b[2] + a[2]*b[0];
	res[2] = a[0]*b[1] - a[1]*b[0];
	return res;
}


void get_bezier_ctrlpoint(POINT3D P1, POINT3D P2, POINT3D N1, POINT3D N2, POINT3D& resp1, POINT3D& resp2)
{
	POINT3D T0, T1, V;
	V = P2 - P1;
	T0 = V - (N1 | V)*N1;
	T1 = V - (N2 | V)*N2;
	double cos0, cos1, k0, k1, sin0, sin1;
	cos0 = (T0 | V) / (T0.norm()*V.norm());
	cos1 = (T1 | V) / (T1.norm()*V.norm());
	sin0 = sqrt(1 - cos0 * cos0);
	sin1 = sqrt(1 - cos1 * cos1);
	k0 = 2.0 / (3 * cos0*(cos0 + 1));
	k1 = 2.0 / (3 * cos1*(cos1 + 1));
	if (abs(abs(cos1) - abs(cos0)) > 0.95)
	{
		k0 = 2 * sin1 * k0 / (sin0 + sin1);
		k1 = 2 * sin0 * k1 / (sin0 + sin1);
	}

	resp1 = P1 + k0 * T0;
	resp2 = P2 - k1 * T1;
}

void get_bezier_ctrlpoint_PN(POINT3D P1, POINT3D P2, POINT3D N1, POINT3D N2, POINT3D& resp1, POINT3D& resp2)
{
	double w12, w21;
	w12 = (P2 - P1) | N1;
	w21 = (P1 - P2) | N2;
	resp1 = (2 * P1 + P2 - w12 * N1) / 3;
	resp2 = (2 * P2 + P1 - w21 * N2) / 3;
}

void get_nor_PN(POINT3D P1, POINT3D P2, POINT3D N1, POINT3D N2, POINT3D& resn) {
	double cof = 2 * ((P2 - P1) | (N1 + N2)) / ((P2 - P1) | (P2 - P1));
	resn = N1 + N2 - cof * (P2 - P1);
	resn = resn.normalize();
}

void get_bezier_point_by_u(POINT3D P1, POINT3D P2, POINT3D P3, POINT3D P4, float u, POINT3D& resp)
{
	resp = P1 * (1 - u)*(1 - u)*(1 - u) + 3 * P2*(1 - u)*(1 - u)*u + 3 * P3*(1 - u)*u*u + P4 * u*u*u;
}

//输入三点及其法向，以及uvw参数，得到输出uvw对应的Bezier面片上的点和法向
void cac_point_and_nor_on_tri_bezier_surface(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D nor1, POINT3D nor2, POINT3D nor3, float u, float v, float w, POINT3D& resp, POINT3D& resn)
{

	POINT3D b003, b012, b021, b030, b120, b210, b300, b201, b102, b111;
	b300 = p1;
	b030 = p2;
	b003 = p3;
	get_bezier_ctrlpoint(b300, b030, nor1, nor2, b210, b120);
	get_bezier_ctrlpoint(b030, b003, nor2, nor3, b021, b012);
	get_bezier_ctrlpoint(b003, b300, nor3, nor1, b102, b201);
	//get_bezier_ctrlpoint_PN(b300, b030, nor1, nor2, b210, b120);
	//get_bezier_ctrlpoint_PN(b030, b003, nor2, nor3, b021, b012);
	//get_bezier_ctrlpoint_PN(b003, b300, nor3, nor1, b102, b201);
	b111 = (b012 + b021 + b120 + b210 + b201 + b102) / 4.0 - (b300 + b030 + b003) / 6.0;

	resp = (b300*u*u*u + b030*v*v*v + b003*w*w*w) + 3 * (b012*v*w*w + b021*v*v*w + b120*u*v*v + b210*u*u*v + b201*u*u*w + b102*u*w*w) + 6 * b111*u*v*w;


	POINT3D derivative1, derivative2;
	// use two test directional derivative to calc nor of resp; 
	float lambda1, mu1, gamma1;
	float lambda2, mu2, gamma2;

	//test data
	lambda1 = -0.5; mu1 = -0.5; gamma1 = 1;
	lambda2 = 0; mu2 = -1; gamma2 = 1;


	//derivative1 = 3 * (b300*u*u*lambda1 + b030 * v*v*mu1 + b003 * w*w*gamma1 + 2 * b111*u*v*gamma1 + 2 * b111*u*mu1*w + 2 * b111*lambda1*v*w + 2 * b210*u*v*lambda1 + b210 * u*u*mu1 + 2 * b120*u*v*mu1 + b120 * v*v*lambda1 + 2 * b012*v*w*gamma1 + b012 * w*w*mu1 + 2 * b021*v*w*mu1 + b021 * v*v*gamma1 + 2 * b201*u*w*lambda1 + b201 * u*u*gamma1 + 2 * b102*u*w*gamma1 + b102 * w*w*lambda1);

	//derivative2 = 3 * (b300*u*u*lambda2 + b030 * v*v*mu2 + b003 * w*w*gamma2 + 2 * b111*u*v*gamma2 + 2 * b111*u*mu2*w + 2 * b111*lambda2*v*w + 2 * b210*u*v*lambda2 + b210 * u*u*mu2 + 2 * b120*u*v*mu2 + b120 * v*v*lambda2 + 2 * b012*v*w*gamma2 + b012 * w*w*mu2 + 2 * b021*v*w*mu2 + b021 * v*v*gamma2 + 2 * b201*u*w*lambda2 + b201 * u*u*gamma2 + 2 * b102*u*w*gamma2 + b102 * w*w*lambda2);


	derivative1 = 1 / 3.0*(b300 + 2 * b210 + b120 + b201 - b003 - b102 - 2 * b012 - b021);
	derivative2 = 1 / 3.0*(b030 + 2 * b120 + b021 + b210 - b003 - b012 - 2 * b102 - b201);

	resn = derivative1.cross(derivative2);
	resn /= resn.norm();


	//derivative1 = 3 * (b300*u*u - b003 * w*w - 2 * b012*v*w - b021 * v*v + b120 * v*v + 2 * b210*u*v + b201 * (2u * w - u * u)
	//	+ b102 * (w*w - 2u * w) + 2 * b111*(v*w - u * v));

	//derivative2 = 3 * (b030*v*v - b003 * w*w + b012*(w*w-2*v*w) + b021*(2*v*w-v*v) + 2*b120 * u*v + b210*u*u - b201 *u*u - 2*b102 * u*w + 2 * b111*(u*w - u * v));


	//resn = derivative1.cross(derivative2);
	//resn /= resn.norm();
	//cout << resn << " "<<endl;
}

void cac_point_and_nor_on_tri_bezier_surface_PN(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D nor1, POINT3D nor2, POINT3D nor3, float u, float v, float w, POINT3D& resp, POINT3D& resn)
{

	POINT3D b003, b012, b021, b030, b120, b210, b300, b201, b102, b111;
	b300 = p1;
	b030 = p2;
	b003 = p3;
	//get_bezier_ctrlpoint(b300, b030, nor1, nor2, b210, b120);
	//get_bezier_ctrlpoint(b030, b003, nor2, nor3, b021, b012);
	//get_bezier_ctrlpoint(b003, b300, nor3, nor1, b102, b201);
	get_bezier_ctrlpoint_PN(b300, b030, nor1, nor2, b210, b120);
	get_bezier_ctrlpoint_PN(b030, b003, nor2, nor3, b021, b012);
	get_bezier_ctrlpoint_PN(b003, b300, nor3, nor1, b102, b201);
	b111 = (b012 + b021 + b120 + b210 + b201 + b102) / 4.0 - (b300 + b030 + b003) / 6.0;

	resp = (b300*u*u*u + b030 * v*v*v + b003 * w*w*w) + 3 * (b012*v*w*w + b021 * v*v*w + b120 * u*v*v + b210 * u*u*v + b201 * u*u*w + b102 * u*w*w) + 6 * b111*u*v*w;


	POINT3D derivative1, derivative2;
	// use two test directional derivative to calc nor of resp; 
	float lambda1, mu1, gamma1;
	float lambda2, mu2, gamma2;

	//test data
	lambda1 = -0.5; mu1 = -0.5; gamma1 = 1;
	lambda2 = 0; mu2 = -1; gamma2 = 1;


	derivative1 = 3 * (b300*u*u*lambda1 + b030 * v*v*mu1 + b003 * w*w*gamma1 + 2 * b111*u*v*gamma1 + 2 * b111*u*mu1*w + 2 * b111*lambda1*v*w + 2 * b210*u*v*lambda1 + b210 * u*u*mu1 + 2 * b120*u*v*mu1 + b120 * v*v*lambda1 + 2 * b012*v*w*gamma1 + b012 * w*w*mu1 + 2 * b021*v*w*mu1 + b021 * v*v*gamma1 + 2 * b201*u*w*lambda1 + b201 * u*u*gamma1 + 2 * b102*u*w*gamma1 + b102 * w*w*lambda1);

	derivative2 = 3 * (b300*u*u*lambda2 + b030 * v*v*mu2 + b003 * w*w*gamma2 + 2 * b111*u*v*gamma2 + 2 * b111*u*mu2*w + 2 * b111*lambda2*v*w + 2 * b210*u*v*lambda2 + b210 * u*u*mu2 + 2 * b120*u*v*mu2 + b120 * v*v*lambda2 + 2 * b012*v*w*gamma2 + b012 * w*w*mu2 + 2 * b021*v*w*mu2 + b021 * v*v*gamma2 + 2 * b201*u*w*lambda2 + b201 * u*u*gamma2 + 2 * b102*u*w*gamma2 + b102 * w*w*lambda2);

	resn = derivative1.cross(derivative2);
	resn /= resn.norm();


	//derivative1 = 3 * (b300*u*u - b003 * w*w - 2 * b012*v*w - b021 * v*v + b120 * v*v + 2 * b210*u*v + b201 * (2u * w - u * u)
	//	+ b102 * (w*w - 2u * w) + 2 * b111*(v*w - u * v));

	//derivative2 = 3 * (b030*v*v - b003 * w*w + b012*(w*w-2*v*w) + b021*(2*v*w-v*v) + 2*b120 * u*v + b210*u*u - b201 *u*u - 2*b102 * u*w + 2 * b111*(u*w - u * v));


	//resn = derivative1.cross(derivative2);
	//resn /= resn.norm();
	//cout << resn << " "<<endl;
}


void cac_two_tri_bezier_surfaces_normals_continuity(PolygonMesh::EdgeIter e_it, PolygonMesh &subvMesh)
{
	PolygonMesh::HalfedgeHandle halfe_it1 = subvMesh.halfedge_handle(e_it, 0);
	PolygonMesh::HalfedgeHandle halfe_it2 = subvMesh.halfedge_handle(e_it, 1);

	POINT3D p1, p2, p3, nor1, nor2, nor3, newp1, newnor1, newp2, newnor2;
	p1 = subvMesh.point(subvMesh.from_vertex_handle(halfe_it1));
	p2 = subvMesh.point(subvMesh.to_vertex_handle(halfe_it1));
	p3 = subvMesh.point(subvMesh.opposite_vh(halfe_it1));
	nor1 = subvMesh.normal(subvMesh.from_vertex_handle(halfe_it1));
	nor2 = subvMesh.normal(subvMesh.to_vertex_handle(halfe_it1));
	nor3 = subvMesh.normal(subvMesh.opposite_vh(halfe_it1));
	cac_point_and_nor_on_tri_bezier_surface(p1, p2, p3, nor1, nor2, nor3, 0.5, 0.5, 0, newp1, newnor1);

	p1 = subvMesh.point(subvMesh.from_vertex_handle(halfe_it2));
	p2 = subvMesh.point(subvMesh.to_vertex_handle(halfe_it2));
	p3 = subvMesh.point(subvMesh.opposite_vh(halfe_it2));
	nor1 = subvMesh.normal(subvMesh.from_vertex_handle(halfe_it2));
	nor2 = subvMesh.normal(subvMesh.to_vertex_handle(halfe_it2));
	nor3 = subvMesh.normal(subvMesh.opposite_vh(halfe_it2));
	cac_point_and_nor_on_tri_bezier_surface(p1, p2, p3, nor1, nor2, nor3, 0.5, 0.5, 0, newp2, newnor2);

	//cout << newp1<<"  " << newnor1<< endl;
	//cout << newp2<<"  " << newnor2<< endl;
	if (abs((newnor1 | newnor2)) < 0.5) {
		cout << "-----法向夹角：" << abs((newnor1 | newnor2)) << endl;
		cout << "point:"<<p1<<" " << p2 << endl;
		cout << "nor:" << nor1 << " " << nor2 << endl;
		cout << "边上法向夹角：" << abs(nor1|nor2)<<endl;
	}
	//POINT3D tmp = (newnor1 | newnor2);

	//cout << "id:"<<(*e_it).idx()<<(newnor1 | newnor2)<< endl;
	two_bezier_triangle_cos = (two_bezier_triangle_cos > abs((newnor1 | newnor2))) ? abs((newnor1 | newnor2)) : two_bezier_triangle_cos;
	if (abs((newnor1 | newnor2)) < 0.9)
	{
		error_normal_cnt++;
	}
	
}

bool find_extraordinary_triangle(PolygonMesh::EdgeIter e_it, PolygonMesh &subvMesh)
{
	PolygonMesh::HalfedgeHandle halfe_it1 = subvMesh.halfedge_handle(e_it, 0);
	PolygonMesh::HalfedgeHandle halfe_it2 = subvMesh.halfedge_handle(e_it, 1);

	POINT3D p1, p2, p3, nor1, nor2, nor3, newp1, newnor1, newp2, newnor2;
	POINT3D N1, N2;
	p1 = subvMesh.point(subvMesh.from_vertex_handle(halfe_it1));
	p2 = subvMesh.point(subvMesh.to_vertex_handle(halfe_it1));
	p3 = subvMesh.point(subvMesh.opposite_vh(halfe_it1));
	double cos1 = (p3 - p1) | (p3 - p2) / ((p3 - p2).norm()*(p3 - p1).norm());
	N1 = (p1 - p2) % (p1 - p3);
	//nor1 = subvMesh.normal(subvMesh.from_vertex_handle(halfe_it1));
	//nor2 = subvMesh.normal(subvMesh.to_vertex_handle(halfe_it1));
	//nor3 = subvMesh.normal(subvMesh.opposite_vh(halfe_it1));
	//N1 = ((nor1 + nor2 + nor3) / 3);

	p1 = subvMesh.point(subvMesh.from_vertex_handle(halfe_it2));
	p2 = subvMesh.point(subvMesh.to_vertex_handle(halfe_it2));
	p3 = subvMesh.point(subvMesh.opposite_vh(halfe_it2));
	double cos2 = (p3 - p1) | (p3 - p2) / ((p3 - p2).norm()*(p3 - p1).norm());
	N2 = (p1 - p2) % (p1 - p3);
	//nor1 = subvMesh.normal(subvMesh.from_vertex_handle(halfe_it2));
	//nor2 = subvMesh.normal(subvMesh.to_vertex_handle(halfe_it2));
	//nor3 = subvMesh.normal(subvMesh.opposite_vh(halfe_it2));
	//N2 = ((nor1 + nor2 + nor3) / 3);

	//cout << "id:"<<(*e_it).idx()<<(newnor1 | newnor2)<< endl;
	double cosN1N2 = abs(N1 | N2) / N1.norm() / N2.norm();
	
	if (cosN1N2<0.1&&(abs(cos1)>0.95||abs(cos2)>0.95))
	{
		if (subvMesh.is_flip_ok(e_it))
		{
			return true;	
		}
		
	}

	return false;	


}

int is_unregular_face(POINT3D P1, POINT3D P2, POINT3D P3)
{
	double cos1,cos2,cos3;
	cos1 = (P1 - P2) | (P1 - P3) / (P1 - P2).norm() / (P1 - P3).norm();
	cos2 = (P2 - P1) | (P2 - P3) / (P2 - P1).norm() / (P2 - P3).norm();
	cos3 = (P3 - P1) | (P3 - P2) / (P3 - P1).norm() / (P3 - P2).norm();

	if (cos1 < -0.3)
	{
		return 1;
	}
	else if (cos2 < -0.3)
	{
		return 2;
	}
	else if (cos3 < -0.3)
	{
		return 3;
	}
	else
		return 0;

}

//use (1/36)*[1,4,1;4,16,4;1,4,1]
void get_quad_bezier_suface_ctrlpoint_and_cac(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn)
{
	POINT3D b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
	b00 = p1;
	b03 = p2;
	b33 = p3;
	b30 = p4;
	// get boundary
	get_bezier_ctrlpoint(b00, b03, nor1, nor2, b01, b02);
	get_bezier_ctrlpoint(b03, b33, nor2, nor3, b13, b23);
	get_bezier_ctrlpoint(b30, b33, nor4, nor3, b31, b32);
	get_bezier_ctrlpoint(b00, b30, nor1, nor4, b10, b20);
	//get_bezier_ctrlpoint_PN(b00, b03, nor1, nor2, b01, b02);
	//get_bezier_ctrlpoint_PN(b03, b33, nor2, nor3, b13, b23);
	//get_bezier_ctrlpoint_PN(b30, b33, nor4, nor3, b31, b32);
	//get_bezier_ctrlpoint_PN(b00, b30, nor1, nor4, b10, b20);

	POINT3D const_b1, const_b2, const_b3, const_b4;
	const_b1 = b00 + 4 * b01 + b02 + 4 * b10 + b20;
	const_b2 = b01 + 4 * b02 + b03 + 4 * b13 + b23;
	const_b3 = b10 + 4 * b20 + b30 + 4 * b31 + b32;
	const_b4 = b13 + 4 * b23 + b31 + 4 * b32 + b33;

	b11 = 116.0 / 2079 * const_b1 + 4.0 / 297 * const_b2 + 4.0 / 297 * const_b3 + 17.0 / 2079 * const_b4;
	b12 = 4.0 / 297 * const_b1 + 116.0 / 2079 * const_b2 + 17.0 / 2079 * const_b3 + 4.0 / 297 * const_b4;
	b21 = 4.0 / 297 * const_b1 + 17.0 / 2079 * const_b2 + 116.0 / 2079 * const_b3 + 4.0 / 297 * const_b4;
	b22 = 17.0 / 2079 * const_b1 + 4.0 / 297 * const_b2 + 4.0 / 297 * const_b3 + 116.0 / 2079 * const_b4;

	resp = 1.0 / 64 * (b00 + 3 * b01 + 3 * b02 + b03 + 3 * (b10 + 3 * b11 + 3 * b12 + b13) + 3 * (b20 + 3 * b21 + 3 * b22 + b23) + b30 + 3 * b31 + 3 * b32 + b33);
	
	POINT3D derivative_u = 3.0 / 32 * (-(b00 + 3 * b01 + 3 * b02 + b03) - 3 * (b10 + 3 * b11 + 3 * b12 + b13) + 3 * (b20 + 3 * b21 + 3 * b22 + b23) + b30 + 3 * b31 + 3 * b32 + b33);

	POINT3D derivative_v = 3.0 / 32 * ((-b00 - 3 * b01 + 3 * b02 + b03) + 3 * (-b10 - 3 * b11 + 3 * b12 + b13) + 3 * (-b20 - 3 * b21 + 3 * b22 + b23) + (-b30 - 3 * b31 + 3 * b32 + b33));

	resn = -(derivative_u.cross(derivative_v));
	resn /= resn.norm();
}

void get_quad_bezier_suface_ctrlpoint_and_cac_PN_quads(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn)
{
	POINT3D b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
	b00 = p1;
	b03 = p2;
	b33 = p3;
	b30 = p4;
	// get boundary
	//get_bezier_ctrlpoint(b00, b03, nor1, nor2, b01, b02);
	//get_bezier_ctrlpoint(b03, b33, nor2, nor3, b13, b23);
	//get_bezier_ctrlpoint(b30, b33, nor4, nor3, b31, b32);
	//get_bezier_ctrlpoint(b00, b30, nor1, nor4, b10, b20);
	get_bezier_ctrlpoint_PN(b00, b03, nor1, nor2, b01, b02);
	get_bezier_ctrlpoint_PN(b03, b33, nor2, nor3, b13, b23);
	get_bezier_ctrlpoint_PN(b30, b33, nor4, nor3, b31, b32);
	get_bezier_ctrlpoint_PN(b00, b30, nor1, nor4, b10, b20);

	POINT3D q, E0, E1, E2, E3, V0, V1, V2, V3;
	q = b01 + b02 + b31 + b32 + b10 + b20 + b13 + b23;
	E0 = 1.0 / 18 * (2 * b10 + 2 * b01 + 2 * q - b32 - b23);
	E1 = 1.0 / 18 * (2 * b20 + 2 * b31 + 2 * q - b02 - b13);
	E2 = 1.0 / 18 * (2 * b32 + 2 * b23 + 2 * q - b10 - b01);
	E3 = 1.0 / 18 * (2 * b02 + 2 * b13 + 2 * q - b31 - b20);
	V0 = 1.0 / 9 * (4 * b00 + 2 * (b03 + b30) + b33);
	V1 = 1.0 / 9 * (4 * b30 + 2 * (b00 + b33) + b03);
	V2 = 1.0 / 9 * (4 * b33 + 2 * (b30 + b03) + b00);
	V3 = 1.0 / 9 * (4 * b03 + 2 * (b00 + b33) + b30);

	//b11 = 3.0 / 2 * E0 - 1.0 / 2 * V0;
	//b21 = 3.0 / 2 * E1 - 1.0 / 2 * V1;
	//b22 = 3.0 / 2 * E2 - 1.0 / 2 * V2;
	//b12 = 3.0 / 2 * E3 - 1.0 / 2 * V3;

	float k = 1;
	b11 = E0 + k * (E0 - V0) / 2.0;
	b21 = E1 + k * (E1 - V1) / 2.0;
	b22 = E2 + k * (E2 - V2) / 2.0;
	b12 = E3 + k * (E3 - V3) / 2.0;

	resp = 1.0 / 64 * (b00 + 3 * b01 + 3 * b02 + b03 + 3 * (b10 + 3 * b11 + 3 * b12 + b13) + 3 * (b20 + 3 * b21 + 3 * b22 + b23) + b30 + 3 * b31 + 3 * b32 + b33);

	POINT3D derivative_u = 3.0 / 32 * (-(b00 + 3 * b01 + 3 * b02 + b03) - 3 * (b10 + 3 * b11 + 3 * b12 + b13) + 3 * (b20 + 3 * b21 + 3 * b22 + b23) + b30 + 3 * b31 + 3 * b32 + b33);

	POINT3D derivative_v = 3.0 / 32 * ((-b00 - 3 * b01 + 3 * b02 + b03) + 3 * (-b10 - 3 * b11 + 3 * b12 + b13) + 3 * (-b20 - 3 * b21 + 3 * b22 + b23) + (-b30 - 3 * b31 + 3 * b32 + b33));

	resn = -(derivative_u.cross(derivative_v));
	resn /= resn.norm();
}

void get_quad_bezier_suface_ctrlpoint_and_cac_degree2(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn)
{
	POINT3D b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
	b00 = p1;
	b03 = p2;
	b33 = p3;
	b30 = p4;
	// get boundary
	get_bezier_ctrlpoint(b00, b03, nor1, nor2, b01, b02);
	get_bezier_ctrlpoint(b03, b33, nor2, nor3, b13, b23);
	get_bezier_ctrlpoint(b30, b33, nor4, nor3, b31, b32);
	get_bezier_ctrlpoint(b00, b30, nor1, nor4, b10, b20);

	POINT3D B1, B2, B3, B4, B5;
	get_bezier_point_by_u(b00, b01, b02, b03, 0.5, B1);
	get_bezier_point_by_u(b03, b13, b23, b33, 0.5, B2);
	get_bezier_point_by_u(b30, b31, b32, b33, 0.5, B3);
	get_bezier_point_by_u(b00, b10, b20, b30, 0.5, B4);
	B5 = (B1 + B2 + B3 + B4) / 4;

	resp = 1.0 / 16 * (p1 + p2 + p3 + p4 + 2 * B1 + 2 * B2 + 2 * B3 + 2 * B4 + 4 * B5);
	POINT3D d1, d2;
	d1 = 1.0 / 4 * (-p1 - p2 + p3 + p4);
	d2 = 1.0 / 4 * (-p1 + p2 + p3 - p4);
	resn = d1.cross(d2);
	resn = -resn / resn.norm();

}

void minimize_diagonal_energy(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn) {
	POINT3D b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
	b00 = p1;
	b03 = p2;
	b33 = p3;
	b30 = p4;
	// get boundary
	get_bezier_ctrlpoint(b00, b03, nor1, nor2, b01, b02);
	get_bezier_ctrlpoint(b03, b33, nor2, nor3, b13, b23);
	get_bezier_ctrlpoint(b30, b33, nor4, nor3, b31, b32);
	get_bezier_ctrlpoint(b00, b30, nor1, nor4, b10, b20);
	//get_bezier_ctrlpoint_PN(b00, b03, nor1, nor2, b01, b02);
	//get_bezier_ctrlpoint_PN(b03, b33, nor2, nor3, b13, b23);
	//get_bezier_ctrlpoint_PN(b30, b33, nor4, nor3, b31, b32);
	//get_bezier_ctrlpoint_PN(b00, b30, nor1, nor4, b10, b20);

	b11 = 1.0 / 27 * (-14 * b00 + 21 * b01 - 3 * b02 - 4 * b03 + 21 * b10 + 6 * b13 - 3 * b20 + 3 * b23 - 4 * b30 + 6 * b31 + 3 * b32 - 5 * b33);
	b12 = 1.0 / 27 * (-4 * b00 - 3 * b01 + 21 * b02 - 14 * b03 + 6 * b10 + 21 * b13 + 3 * b20 - 3 * b23 - 5 * b30 + 3 * b31 + 6 * b32 - 4 * b33);
	b21 = 1.0 / 27 * (-4 * b00 + 6 * b01 + 3 * b02 - 5 * b03 - 3 * b10 + 3 * b13 + 21 * b20 + 6 * b23 - 14 * b30 + 21 * b31 - 3 * b32 - 4 * b33);
	b22 = 1.0 / 27 * (-5 * b00 + 3 * b01 + 6 * b02 - 4 * b03 + 3 * b10 - 3 * b13 + 6 * b20 + 21 * b23 - 4 * b30 - 3 * b31 + 21 * b32 - 14 * b33);

	resp = 1.0 / 64 * (b00 + 3 * b01 + 3 * b02 + b03 + 3 * (b10 + 3 * b11 + 3 * b12 + b13) + 3 * (b20 + 3 * b21 + 3 * b22 + b23) + b30 + 3 * b31 + 3 * b32 + b33);

	POINT3D derivative_u = 3.0 / 32 * (-(b00 + 3 * b01 + 3 * b02 + b03) - 3 * (b10 + 3 * b11 + 3 * b12 + b13) + 3 * (b20 + 3 * b21 + 3 * b22 + b23) + b30 + 3 * b31 + 3 * b32 + b33);

	POINT3D derivative_v = 3.0 / 32 * ((-b00 - 3 * b01 + 3 * b02 + b03) + 3 * (-b10 - 3 * b11 + 3 * b12 + b13) + 3 * (-b20 - 3 * b21 + 3 * b22 + b23) + (-b30 - 3 * b31 + 3 * b32 + b33));

	resn = -(derivative_u.cross(derivative_v));
	resn /= resn.norm();
}




/* Local Method */

vector<vector<POINT3D>> get_classic_PN_quads_ctrlP(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4) {
	POINT3D b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
	b00 = p1;
	b03 = p2;
	b33 = p3;
	b30 = p4;
	// get boundary
	get_bezier_ctrlpoint(b00, b03, nor1, nor2, b01, b02);
	get_bezier_ctrlpoint(b03, b33, nor2, nor3, b13, b23);
	get_bezier_ctrlpoint(b30, b33, nor4, nor3, b31, b32);
	get_bezier_ctrlpoint(b00, b30, nor1, nor4, b10, b20);
	//get_bezier_ctrlpoint_PN(b00, b03, nor1, nor2, b01, b02);
	//get_bezier_ctrlpoint_PN(b03, b33, nor2, nor3, b13, b23);
	//get_bezier_ctrlpoint_PN(b30, b33, nor4, nor3, b31, b32);
	//get_bezier_ctrlpoint_PN(b00, b30, nor1, nor4, b10, b20);

	POINT3D q, E0, E1, E2, E3, V0, V1, V2, V3;
	q = b01 + b02 + b31 + b32 + b10 + b20 + b13 + b23;
	E0 = 1.0 / 18 * (2 * b10 + 2 * b01 + 2 * q - b32 - b23);
	E1 = 1.0 / 18 * (2 * b20 + 2 * b31 + 2 * q - b02 - b13);
	E2 = 1.0 / 18 * (2 * b32 + 2 * b23 + 2 * q - b10 - b01);
	E3 = 1.0 / 18 * (2 * b02 + 2 * b13 + 2 * q - b31 - b20);
	V0 = 1.0 / 9 * (4 * b00 + 2 * (b03 + b30) + b33);
	V1 = 1.0 / 9 * (4 * b30 + 2 * (b00 + b33) + b03);
	V2 = 1.0 / 9 * (4 * b33 + 2 * (b30 + b03) + b00);
	V3 = 1.0 / 9 * (4 * b03 + 2 * (b00 + b33) + b30);

	b11 = 3.0 / 2 * E0 - 1.0 / 2 * V0;
	b21 = 3.0 / 2 * E1 - 1.0 / 2 * V1;
	b22 = 3.0 / 2 * E2 - 1.0 / 2 * V2;
	b12 = 3.0 / 2 * E3 - 1.0 / 2 * V3;

	vector<vector<POINT3D>>ret(4);
	ret[0].push_back(b00); ret[0].push_back(b01); ret[0].push_back(b02); ret[0].push_back(b03);
	ret[1].push_back(b10); ret[1].push_back(b11); ret[1].push_back(b12); ret[1].push_back(b13);
	ret[2].push_back(b20); ret[2].push_back(b21); ret[2].push_back(b22); ret[2].push_back(b23);
	ret[3].push_back(b30); ret[3].push_back(b31); ret[3].push_back(b32); ret[3].push_back(b33);

	return ret;
}

vector<POINT3D> get_degree_elevation_ctrlP_curve(vector<POINT3D>& p, int deg) {
	vector<POINT3D>ret;
	ret.push_back(p[0]);
	for (int i = 1; i <= deg; ++i) {
		ret.push_back(i*1.0 / (deg + 1)*p[i - 1] + (1 - i*1.0 / (deg + 1))*p[i]);
	}
	ret.push_back(p[deg]);

	return ret;
}

vector<vector<POINT3D>> get_degree_elevation_ctrlP_surface(vector<vector<POINT3D>>& p, int deg) {
	vector<vector<POINT3D>>tmpv;
	
	// 横向进行升阶
	for (auto &v : p) {
		tmpv.push_back(get_degree_elevation_ctrlP_curve(v, deg));
	}

	vector<vector<POINT3D>>ret(deg + 2, vector<POINT3D>(deg + 2));

	// 纵向进行升阶
	for (int j = 0; j < deg + 2; ++j) {
		vector<POINT3D>tmpColVec,tmpEleVec;
		for (int i = 0; i < deg + 1; ++i) {
			tmpColVec.push_back(tmpv[i][j]);
		}
		tmpEleVec = get_degree_elevation_ctrlP_curve(tmpColVec, deg);

		for (int k = 0; k < deg + 2; ++k) {
			ret[k][j] = tmpEleVec[k];
		}
	}

	return ret;
}

void cacCurve(POINT3D p1, POINT3D p2, POINT3D n1, POINT3D n2, vector<POINT3D>&p) {
	Eigen::MatrixXd H(12, 12);
	Eigen::Matrix4d block(4, 4);
	Eigen::Matrix4d zero = Eigen::Matrix4d::Zero(4, 4);
	block << 2, -1, 0, 0,
		-1, 2, -1, 0,
		0, -1, 2, -1,
		0, 0, -1, 2;
	H << block, zero, zero,
		zero, block, zero,
		zero, zero, block;

	H = H * 2;

	double p1x = p1[0], p2x = p2[0], n1x = n1[0], n2x = n2[0];
	double p1y = p1[1], p2y = p2[1], n1y = n1[1], n2y = n2[1];
	double p1z = p1[2], p2z = p2[2], n1z = n1[2], n2z = n2[2];

	Eigen::MatrixXd A(4, 12);
	A << n1x, 0, 0, 0, n1y, 0, 0, 0, n1z, 0, 0, 0,
		0, 0, 0, n2x, 0, 0, 0, n2y, 0, 0, 0, n2z,
		-n1x, n1x, 0, 0, -n1y, n1y, 0, 0, -n1z, n1z, 0, 0,
		0, 0, -n2x, n2x, 0, 0, -n2y, n2y, 0, 0, -n2z, n2z;


	Eigen::MatrixXd c(12, 1);
	Eigen::MatrixXd b(4, 1);

	c << -2 * p1x, 0, 0, -2 * p2x, -2 * p1y, 0, 0, -2 * p2y, -2 * p1z, 0, 0, -2 * p2z;
	b << p1x * n1x + p1y * n1y + p1z * n1z, p2x*n2x + p2y * n2y + p2z * n2z, 0, 0;

	solution s;
	//vector<Eigen::MatrixXd>ret;
	//ret = s.eigen_elimination(A, b);
	//Eigen::MatrixXd A1, b1;
	//A1 = ret[0]; b1 = ret[1];

	Eigen::MatrixXd x, lambda; //x1,x2,x3,x4,...z3,z4;
	s.qp_lagrange(H, c, A, b, x, lambda, 12, A.rows());


	p.push_back(POINT3D(p1x, p1y, p1z));
	p.push_back(POINT3D(x(0, 0), x(4, 0), x(8, 0)));
	p.push_back(POINT3D(x(1, 0), x(5, 0), x(9, 0)));
	p.push_back(POINT3D(x(2, 0), x(6, 0), x(10, 0)));
	p.push_back(POINT3D(x(3, 0), x(7, 0), x(11, 0)));
	p.push_back(POINT3D(p2x, p2y, p2z));
}

void cacNorField(POINT3D p1, POINT3D p2, POINT3D n1, POINT3D n2, vector<POINT3D>&p, vector<POINT3D>&n) {
	Eigen::MatrixXd H = Eigen::MatrixXd::Identity(7, 7);
	H = H * 2;

	double p1x = p1[0], p2x = p2[0], n1x = n1[0], n2x = n2[0];
	double p1y = p1[1], p2y = p2[1], n1y = n1[1], n2y = n2[1];
	double p1z = p1[2], p2z = p2[2], n1z = n1[2], n2z = n2[2];

	double tx = 0.5*(n1x + n2x);
	double ty = 0.5*(n1y + n2y);
	double tz = 0.5*(n1z + n2z);

	//POINT3D n1(n1x, n1y, n1z);
	//POINT3D n2(n2x, n2y, n2z);



	Eigen::MatrixXd A(5, 7);
	A << (p[1][0] - p[0][0]), (p[1][1] - p[0][1]), (p[1][2] - p[0][2]), (p[3] - p[2]) | n1, 0, 0, 0,
		6 * (p[2][0] - p[1][0]), 6 * (p[2][1] - p[1][1]), 6 * (p[2][2] - p[1][2]), (p[4] - p[3]) | n1, (p[3] - p[2]) * 6 | n1, (p[1] - p[0]) | n2, 0,
		36 * (p[3][0] - p[2][0]), 36 * (p[3][1] - p[2][1]), 36 * (p[3][2] - p[2][2]), (p[5] - p[4]) | n1, (p[4] - p[3]) * 16 | n1, (p[2] - p[1]) * 16 | n2, (p[1] - p[0]) | n2,
		6 * (p[4][0] - p[3][0]), 6 * (p[4][1] - p[3][1]), 6 * (p[4][2] - p[3][2]), 0, (p[5] - p[4]) | n1, (p[3] - p[2]) * 6 | n2, (p[2] - p[1]) | n2,
		(p[5][0] - p[4][0]), (p[5][1] - p[4][1]), (p[5][2] - p[4][2]), 0, 0, 0, (p[3] - p[2]) | n2;

	//if ((n1 | n2) > 0.95) {
	//	cout << A;
	//	

	//}


	Eigen::MatrixXd c(7, 1);
	Eigen::MatrixXd b(5, 1);

	c << -2 * tx, -2 * ty, -2 * tz, -2, -2, -2, -2;
	b << 0, 0, 0, 0, 0;

	solution s;
	//vector<Eigen::MatrixXd>ret;
	//ret = s.eigen_elimination(A, b);
	//Eigen::MatrixXd A1, b1;
	//A1 = ret[0]; b1 = ret[1];


	Eigen::MatrixXd x, lambda; //x,y,z,s0,s1,s3,s4;
	s.qp_lagrange(H, c, A, b, x, lambda, 7, A.rows()); //5

	n.push_back(n1*x(3, 0));
	n.push_back(n1*x(4, 0));
	n.push_back(POINT3D(x(0, 0), x(1, 0), x(2, 0)));
	n.push_back(n2*x(5, 0));
	n.push_back(n2*x(6, 0));
}

void slightDisturb(POINT3D &n) {
	double theta = std::acos(n[2]);
	double pha = n[0]!=0?std::atan(n[1] / n[0]):-1;

	double offset = 0.0005;
	int num = std::rand() % 10;
	int flg = num % 2 ? 1 : -1;

	if (num == 0)num = 1;
	double scale = num / 10.0;
	theta += flg * scale*offset;

	cout << n << endl;
	if (pha == -1) {
		n[0] = 0;
		n[1] = (n[1] > 0) ? std::sin(theta) : -std::sin(theta);
		n[2] = std::cos(theta);
	}
	else {
		num = std::rand() % 10;
		flg = num % 2 ? 1 : -1;
		scale = num / 10.0;

		pha += flg * scale*offset;

		n[0] = std::sin(theta)*std::cos(pha);
		n[1] = std::sin(theta)*std::sin(pha);
		n[2] = std::cos(theta);
	}

	cout << n << endl;
	cout << "===========" << endl;
}

void preNormalDisjustion(T_PolygonMesh *subvMesh) {
	//for (auto it = subvMesh->edges_begin(); it != subvMesh->edges_end(); ++it) {
	//	auto tmp_heh1 = subvMesh->halfedge_handle(it, 0);
	//	auto tmp_vh1 = subvMesh->from_vertex_handle(tmp_heh1);

	//	auto nor1 = subvMesh->normal(tmp_vh1);

	//	POINT3D p1, p2, n1, n2;
	//	n1[0] = nor1[0]; 
	//	n1[1] = nor1[1]; 
	//	n1[2] = nor1[2];

	//	auto tmp_heh2 = subvMesh->halfedge_handle(it, 1);
	//	auto tmp_vh2 = subvMesh->from_vertex_handle(tmp_heh2);

	//	auto nor2 = subvMesh->normal(tmp_vh2);
	//	n2[0] = nor2[0];
	//	n2[1] = nor2[1];
	//	n2[2] = nor2[2];

	//	if (n1 == n2) {
	//		slightDisturb(n1);
	//	}

	//	T_PolygonMesh::Normal  n;
	//	n[0] = n1[0]; n[1] = n1[1]; n[2] = n1[2];
	//	subvMesh->set_normal(tmp_vh1, n);
	//}

	int cnt = 0;
	for (auto v_it = subvMesh->vertices_begin(); v_it != subvMesh->vertices_end(); v_it++) {
		POINT3D tmpp;
		T_PolygonMesh::Normal tt = subvMesh->normal(v_it);
		if (cnt % 2) {
			tmpp[0] = 0.5*tt[0] + 0.866025*tt[1];
			tmpp[1] = -0.866025*tt[0] + 0.5*tt[1];
			tmpp[2] = tt[2];
		}
		else {
			tmpp[0] = 0.5*tt[0] - 0.866025*tt[1];
			tmpp[1] = 0.866025*tt[0] + 0.5*tt[1];
			tmpp[2] = tt[2];
		}
		cnt++;
		T_PolygonMesh::Normal n;
		n[0] = tmpp[0];
		n[1] = tmpp[1];
		n[2] = tmpp[2];
		subvMesh->set_normal(v_it, n);
	}

	cout << "++++++++++++++" << endl;


}


void cacBoundaryInfo(POINT3D p1, POINT3D p2, POINT3D n1, POINT3D n2, vector<POINT3D>&p, vector<POINT3D>&n) {
	if ((n1 | n2)>0.995) {
		testFlg = true;
		POINT3D resp1,resp2;
		get_bezier_ctrlpoint_PN(p1, p2, n1, n2, resp1, resp2);
		std::vector<POINT3D>bezierVec;
		bezierVec.push_back(p1);
		bezierVec.push_back(resp1);
		bezierVec.push_back(resp2);
		bezierVec.push_back(p2);
		bezierVec = get_degree_elevation_ctrlP_curve(bezierVec, 3);
		bezierVec = get_degree_elevation_ctrlP_curve(bezierVec, 4);

		p = bezierVec;
		
		POINT3D nor;
		get_nor_PN(p1, p2, n1, n2, nor);
		std::vector<POINT3D>bezierNorVec;
		bezierNorVec.push_back(n1);
		bezierNorVec.push_back(nor);
		bezierNorVec.push_back(n2);

		bezierNorVec = get_degree_elevation_ctrlP_curve(bezierNorVec, 2);
		bezierNorVec = get_degree_elevation_ctrlP_curve(bezierNorVec, 3);

		n = bezierNorVec;
	}
	else {
		cacCurve(p1, p2, n1, n2, p);
		cacNorField(p1, p2, n1, n2, p, n);
	}
}

// p10 为边界左上一层的point,p15则是右上一层的边界point
vector<Eigen::MatrixXd> initConstraintsEigen(vector<POINT3D>&bdPoints, vector<POINT3D>&bdNors, POINT3D p10, POINT3D p15) {
	POINT3D p0, p1, p2, p3, p4, p5;
	POINT3D n0, n1, n2, n3, n4;

	p0 = bdPoints[0]; p1 = bdPoints[1]; p2 = bdPoints[2]; p3 = bdPoints[3]; p4 = bdPoints[4]; p5 = bdPoints[5];
	n0 = bdNors[0]; n1 = bdNors[1]; n2 = bdNors[2]; n3 = bdNors[3]; n4 = bdNors[4];

	Eigen::MatrixXd block(8, 12);
	Eigen::MatrixXd b(8, 1);

	block << n0[0], n0[1], n0[2], 0, 0, 0, 0, 0, 0, 0, 0, 0,
		10 * n1[0], 10 * n1[1], 10 * n1[2], 5 * n0[0], 5 * n0[1], 5 * n0[2], 0, 0, 0, 0, 0, 0,
		15 * n2[0], 15 * n2[1], 15 * n2[2], 20 * n1[0], 20 * n1[1], 20 * n1[2], 5 * n0[0], 5 * n0[1], 5 * n0[2], 0, 0, 0,
		20 * n3[0], 20 * n3[1], 20 * n3[2], 60 * n2[0], 60 * n2[1], 60 * n2[2], 40 * n1[0], 40 * n1[1], 40 * n1[2], 5 * n0[0], 5 * n0[1], 5 * n0[2],
		5 * n4[0], 5 * n4[1], 5 * n4[2], 40 * n3[0], 40 * n3[1], 40 * n3[2], 60 * n2[0], 60 * n2[1], 60 * n2[2], 20 * n1[0], 20 * n1[1], 20 * n1[2],
		0, 0, 0, 5 * n4[0], 5 * n4[1], 5 * n4[2], 20 * n3[0], 20 * n3[1], 20 * n3[2], 15 * n2[0], 15 * n2[1], 15 * n2[2],
		0, 0, 0, 0, 0, 0, 5 * n4[0], 5 * n4[1], 5 * n4[2], 10 * n3[0], 10 * n3[1], 10 * n3[2],
		0, 0, 0, 0, 0, 0, 0, 0, 0, n4[0], n4[1], n4[2];

	b << (p1 | n0), 5 * (p2 | n0) + 10 * (p1 | n1) - 3 * ((p10 - p0) | n2), 5 * (p3 | n0) + 20 * (p2 | n1) + 15 * (p1 | n2) - 2 * ((p10 - p0) | n3),
		5 * (p4 | n0) + 40 * (p3 | n1) + 60 * (p2 | n2) + 20 * (p1 | n3) - ((p10 - p0) | n4), -((p15 - p5) | n0) + 20 * (p4 | n1) + 60 * (p3 | n2) + 40 * (p2 | n3) + 5 * (p1 | n4), -2 * ((p15 - p5) | n1) + 15 * (p4 | n2) + 20 * (p3 | n3) + 5 * (p2 | n4), -3 * ((p15 - p5) | n2) + 10 * (p4 | n3) + 5 * (p3 | n4), (p4 | n4);

	vector<Eigen::MatrixXd>ret;
	ret.push_back(block);
	ret.push_back(b);

	return ret;
}

void get_quad_bezier_ctrlpoint_local_method(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn) {

	POINT3D b50, b51, b52, b53, b54, b55;
	POINT3D b40, b41, b42, b43, b44, b45;
	POINT3D b30, b31, b32, b33, b34, b35;
	POINT3D b20, b21, b22, b23, b24, b25;
	POINT3D b10, b11, b12, b13, b14, b15;
	POINT3D b00, b01, b02, b03, b04, b05;

	POINT3D n40, n41, n42, n43, n44;
	POINT3D n30,                n34;
	POINT3D n20,				n24;
	POINT3D n10,				n14;
	POINT3D n00, n01, n02, n03, n04;

	b00 = p1; b05 = p2; b55 = p3; b50 = p4;
	n00 = nor1; n04 = nor2; n44 = nor3; n40 = nor4;

	Eigen::MatrixXd block1(8, 12), block2(8, 12), block3(8, 12), block4(8, 12);
	Eigen::MatrixXd b1(8, 1), b2(8, 1), b3(8, 1), b4(8, 1);
	vector<Eigen::MatrixXd>tmpvec;

	//caculate boundary curve and normal fields
	vector<POINT3D>bdPsB, bdNorsB;
	vector<POINT3D>bdPsR, bdNorsR;
	vector<POINT3D>bdPsU, bdNorsU;
	vector<POINT3D>bdPsL, bdNorsL;

	//if (abs(p1[0] - 0.046591) < 0.0001&&abs(p2[0] - 0.023227) < 0.0001&&abs(p3[0] - 0.022918) < 0.0001) {
	//	cout << endl;
	//	testFlg = true;
	//}

	//bottom
	cacBoundaryInfo(p1, p2, nor1, nor2, bdPsB, bdNorsB);
	b01 = bdPsB[1]; b02 = bdPsB[2]; b03 = bdPsB[3]; b04 = bdPsB[4];
	n00 = bdNorsB[0]; n01 = bdNorsB[1]; n02 = bdNorsB[2]; n03 = bdNorsB[3];

	//right
	cacBoundaryInfo(p2, p3, nor2, nor3, bdPsR, bdNorsR);
	b15 = bdPsR[1]; b25 = bdPsR[2]; b35 = bdPsR[3]; b45 = bdPsR[4];
	n04 = bdNorsR[0]; n14 = bdNorsR[1]; n24 = bdNorsR[2]; n34 = bdNorsR[3];

	//up
	cacBoundaryInfo(p3, p4, nor3, nor4, bdPsU, bdNorsU);
	b54 = bdPsU[1]; b53 = bdPsU[2]; b52 = bdPsU[3]; b51 = bdPsU[4];
	n44 = bdNorsU[0]; n43 = bdNorsU[1]; n42 = bdNorsU[2]; n41 = bdNorsU[3];

	//left
	cacBoundaryInfo(p4, p1, nor4, nor1, bdPsL, bdNorsL);
	b40 = bdPsL[1]; b30 = bdPsL[2]; b20 = bdPsL[3]; b10 = bdPsL[4];
	n40 = bdNorsL[0]; n30 = bdNorsL[1]; n20 = bdNorsL[2]; n10 = bdNorsL[3];


	tmpvec = initConstraintsEigen(bdPsB, bdNorsB, b10, b15);
	block1 = tmpvec[0]; b1 = tmpvec[1];
	tmpvec = initConstraintsEigen(bdPsR, bdNorsR, b04, b54);
	block2 = tmpvec[0]; b2 = tmpvec[1];
	tmpvec = initConstraintsEigen(bdPsU, bdNorsU, b45, b40);
	block3 = tmpvec[0]; b3 = tmpvec[1];
	tmpvec = initConstraintsEigen(bdPsL, bdNorsL, b51, b01);
	block4 = tmpvec[0]; b4 = tmpvec[1];

	//cac classical PN quads
	vector<vector<POINT3D>>t;

	t = get_classic_PN_quads_ctrlP(p1, p2, p3, p4, nor1, nor2, nor3, nor4);
	t = get_degree_elevation_ctrlP_surface(t, 3);
	t = get_degree_elevation_ctrlP_surface(t, 4);


	// init A H b c
	Eigen::MatrixXd H = Eigen::MatrixXd::Identity(36, 36); // b11x,b11y,b11z...b14x,b14y,b14z,b24x...b44z,b43x,b43y,...b41z,b31x,..b21x,b21y,b21z;

	H = H * 2;

	Eigen::MatrixXd A(32, 36);


	Eigen::MatrixXd zero24 = Eigen::MatrixXd::Zero(8, 24);
	Eigen::MatrixXd zero9 = Eigen::MatrixXd::Zero(8, 9);
	Eigen::MatrixXd zero15 = Eigen::MatrixXd::Zero(8, 15);
	Eigen::MatrixXd zero18 = Eigen::MatrixXd::Zero(8, 18);
	Eigen::MatrixXd zero6 = Eigen::MatrixXd::Zero(8, 6);

	Eigen::MatrixXd block21 = block2.block<8, 3>(0, 0);
	Eigen::MatrixXd block22 = block2.block<8, 3>(0, 3);
	Eigen::MatrixXd block23 = block2.block<8, 3>(0, 6);
	Eigen::MatrixXd block24 = block2.block<8, 3>(0, 9);

	Eigen::MatrixXd block31 = block3.block<8, 3>(0, 0);
	Eigen::MatrixXd block32 = block3.block<8, 3>(0, 3);
	Eigen::MatrixXd block33 = block3.block<8, 3>(0, 6);
	Eigen::MatrixXd block34 = block3.block<8, 3>(0, 9);

	Eigen::MatrixXd block41 = block4.block<8, 3>(0, 0);
	Eigen::MatrixXd block42 = block4.block<8, 3>(0, 3);
	Eigen::MatrixXd block43 = block4.block<8, 3>(0, 6);
	Eigen::MatrixXd block44 = block4.block<8, 3>(0, 9);

	A << block1, zero24,
		zero9, block2, zero15,
		zero18, block3, zero6,
		block44, zero24, block41, block42, block43;

	//A << n00[0], n00[1], n00[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	10 * n01[0], 10 * n01[1], 10 * n01[2], 5 * n00[0], 5 * n00[1], 5 * n00[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	15 * n02[0], 15 * n02[1], 15 * n02[2], 20 * n01[0], 20 * n01[1], 20 * n01[2], 5 * n00[0], 5 * n00[1], 5 * n00[2], 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	20 * n03[0], 20 * n03[1], 20 * n03[2], 60 * n02[0], 60 * n02[1], 60 * n02[2], 40 * n01[0], 40 * n01[1], 40 * n01[2], 5 * n00[0], 5 * n00[1], 5 * n00[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	5 * n04[0], 5 * n04[1], 5 * n04[2], 40 * n03[0], 40 * n03[1], 40 * n03[2], 60 * n02[0], 60 * n02[1], 60 * n02[2], 20 * n01[0], 20 * n01[1], 20 * n01[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 5 * n04[0], 5 * n04[1], 5 * n04[2], 20 * n03[0], 20 * n03[1], 20 * n03[2], 15 * n02[0], 15 * n02[1], 15 * n02[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 5 * n04[0], 5 * n04[1], 5 * n04[2], 10 * n03[0], 10 * n03[1], 10 * n03[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, n04[0], n04[1], n04[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, n04[0], n04[1], n04[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 10 * n14[0], 10 * n14[1], 10 * n14[2], 5 * n04[0], 5 * n04[1], 5 * n04[2], 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 15 * n24[0], 15 * n24[1], 15 * n24[2], 20 * n14[0], 20 * n14[1], 20 * n14[2], 5 * n04[0], 5 * n04[1], 5 * n04[2],
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 20 * n34[0], 20 * n34[1], 20 * n34[2], 60 * n24[0], 60 * n24[1], 60 * n24[2], 40 * n14[0], 40 * n14[1], 40 * n14[2],
	//	5 * n04[0], 5 * n04[1], 5 * n04[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 5 * n44[0], 5 * n44[1], 5 * n44[2], 40 * n34[0], 40 * n34[1], 40 * n34[2], 60 * n24[0], 60 * n24[1], 60 * n24[2],
	//	20 * n14[0], 20 * n14[1], 20 * n14[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5 * n44[0], 5 * n44[1], 5 * n44[2], 20 * n34[0], 20 * n34[1], 20 * n34[2],
	//	15 * n24[0], 15 * n24[1], 15 * n24[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5 * n44[0], 5 * n44[1], 5 * n44[2],
	//	10 * n34[0], 10 * n34[1], 10 * n34[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	n44[0], n44[1], n44[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	n44[0], n44[1], n44[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	10 * n43[0], 10 * n43[1], 10 * n43[2], 5 * n44[0], 5 * n44[1], 5 * n44[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	15 * n42[0], 15 * n42[1], 15 * n42[2], 20 * n43[0], 20 * n43[1], 20 * n43[2], 5 * n44[0], 5 * n44[1], 5 * n44[2], 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	20 * n41[0], 20 * n41[1], 20 * n41[2], 60 * n42[0], 60 * n42[1], 60 * n42[2], 40 * n43[0], 40 * n43[1], 40 * n43[2], 5 * n44[0], 5 * n44[1], 5 * n44[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	5 * n40[0], 5 * n40[1], 5 * n40[2], 40 * n41[0], 40 * n41[1], 40 * n41[2], 60 * n42[0], 60 * n42[1], 60 * n42[2], 20 * n43[0], 20 * n43[1], 20 * n43[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 5 * n40[0], 5 * n40[1], 5 * n40[2], 20 * n41[0], 20 * n41[1], 20 * n41[2], 15 * n42[0], 15 * n42[1], 15 * n42[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 5 * n40[0], 5 * n40[1], 5 * n40[2], 10 * n41[0], 10 * n41[1], 10 * n41[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, n40[0], n40[1], n40[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, n40[0], n40[1], n40[2], 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 10 * n30[0], 10 * n30[1], 10 * n30[2], 5 * n40[0], 5 * n40[1], 5 * n40[2], 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 15 * n20[0], 15 * n20[1], 15 * n20[2], 20 * n30[0], 20 * n30[1], 20 * n30[2], 5 * n40[0], 5 * n40[1], 5 * n40[2],
	//	5 * n40[0], 5 * n40[1], 5 * n40[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 20 * n10[0], 20 * n10[1], 20 * n10[2], 60 * n20[0], 60 * n20[1], 60 * n20[2], 40 * n30[0], 40 * n30[1], 40 * n30[2],
	//	20 * n30[0], 20 * n30[1], 20 * n30[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 5 * n00[0], 5 * n00[1], 5 * n00[2], 40 * n10[0], 40 * n10[1], 40 * n10[2], 60 * n20[0], 60 * n20[1], 60 * n20[2],
	//	15 * n20[0], 15 * n20[1], 15 * n20[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5 * n00[0], 5 * n00[1], 5 * n00[2], 20 * n10[0], 20 * n10[1], 20 * n10[2],
	//	10 * n10[0], 10 * n10[1], 10 * n10[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5 * n00[0], 5 * n00[1], 5 * n00[2],
	//	n00[0], n00[1], n00[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;


	Eigen::MatrixXd c(36, 1);
	Eigen::MatrixXd b(32, 1);

	c << -2 * t[1][1][0], -2 * t[1][1][1], -2 * t[1][1][2], -2 * t[1][2][0], -2 * t[1][2][1], -2 * t[1][2][2],
		-2 * t[1][3][0], -2 * t[1][3][1], -2 * t[1][3][2], -2 * t[1][4][0], -2 * t[1][4][1], -2 * t[1][4][2],
		-2 * t[2][4][0], -2 * t[2][4][1], -2 * t[2][4][2], -2 * t[3][4][0], -2 * t[3][4][1], -2 * t[3][4][2],
		-2 * t[4][4][0], -2 * t[4][4][1], -2 * t[4][4][2], -2 * t[4][3][0], -2 * t[4][3][1], -2 * t[4][3][2],
		-2 * t[4][2][0], -2 * t[4][2][1], -2 * t[4][2][2], -2 * t[4][1][0], -2 * t[4][1][1], -2 * t[4][1][2],
		-2 * t[3][1][0], -2 * t[3][1][1], -2 * t[3][1][2], -2 * t[2][1][0], -2 * t[2][1][1], -2 * t[2][1][2];


	b << b1, b2, b3, b4;

	if (abs(p1[0] - 0.046591) < 0.0001&&abs(p2[0] - 0.023227) < 0.0001&&abs(p3[0] - 0.022918) < 0.0001) {
		cout << A << endl;
		cout << b << endl;
		cout << c << endl;
	}


	// 满秩处理
	Eigen::MatrixXd Ablock1 = A.block<7, 36>(0, 0);
	Eigen::MatrixXd Ablock2 = A.block<7, 36>(8, 0);
	Eigen::MatrixXd Ablock3 = A.block<7, 36>(16, 0);
	Eigen::MatrixXd Ablock4 = A.block<7, 36>(24, 0);

	Eigen::MatrixXd bblock1 = b.block<7, 1>(0, 0);
	Eigen::MatrixXd bblock2 = b.block<7, 1>(8, 0);
	Eigen::MatrixXd bblock3 = b.block<7, 1>(16, 0);
	Eigen::MatrixXd bblock4 = b.block<7, 1>(24, 0);

	Eigen::MatrixXd Alast(28, 36), blast(28, 1);
	Alast << Ablock1, Ablock2, Ablock3, Ablock4;
	blast << bblock1, bblock2, bblock3, bblock4;

	// cac point & normal
	solution s;

	vector<Eigen::MatrixXd>rett;
	rett = s.eigen_elimination(Alast, blast);

	Eigen::MatrixXd Alast1, blast1;
	Alast1 = rett[0];
	blast1 = rett[1];

	Alast1 = Alast;
	blast1 = blast;

	Eigen::MatrixXd x, lambda;
	s.qp_lagrange(H, c, Alast1, blast1, x, lambda, 36, 28);

	b11[0] = x(0, 0); b11[1] = x(1, 0); b11[2] = x(2, 0);
	b12[0] = x(3, 0); b12[1] = x(4, 0); b12[2] = x(5, 0);
	b13[0] = x(6, 0); b13[1] = x(7, 0); b13[2] = x(8, 0);
	b14[0] = x(9, 0); b14[1] = x(10, 0); b14[2] = x(11, 0);
	b24[0] = x(12, 0); b24[1] = x(13, 0); b24[2] = x(14, 0);
	b34[0] = x(15, 0); b34[1] = x(16, 0); b34[2] = x(17, 0);
	b44[0] = x(18, 0); b44[1] = x(19, 0); b44[2] = x(20, 0);
	b43[0] = x(21, 0); b43[1] = x(22, 0); b43[2] = x(23, 0);
	b42[0] = x(24, 0); b42[1] = x(25, 0); b42[2] = x(26, 0);
	b41[0] = x(27, 0); b41[1] = x(28, 0); b41[2] = x(29, 0);
	b31[0] = x(30, 0); b31[1] = x(31, 0); b31[2] = x(32, 0);
	b21[0] = x(33, 0); b21[1] = x(34, 0); b21[2] = x(35, 0);

	if (testFlg) {
		b11 = t[1][1]; b12 = t[1][2]; b13 = t[1][3]; b14 = t[1][4];
		b24 = t[2][4]; b34 = t[3][4]; b44 = t[4][4]; b43 = t[4][3];
		b42 = t[4][2]; b41 = t[4][1]; b31 = t[3][1]; b21 = t[2][1];
		testFlg = false;
	}


	b22 = t[2][2]; b23 = t[2][3]; b32 = t[3][2]; b33 = t[3][3];

	//cac newp & newNor
	double f = 1.0 / 1024;
	resp = f * (b00 + 5 * b01 + 10 * b02 + 10 * b03 + 5 * b04 + b05 + 5 * (b10 + 5 * b11 + 10 * b12 + 10 * b13 + 5 * b14 + b15) + 10 * (b20 + 5 * b21 + 10 * b22 + 10 * b23 + 5 * b24 + b25) + 10 * (b30 + 5 * b31 + 10 * b32 + 10 * b33 + 5 * b34 + b35) + 5 * (b40 + 5 * b41 + 10 * b42 + 10 * b43 + 5 * b44 + b45) + (b50 + 5 * b51 + 10 * b52 + 10 * b53 + 5 * b54 + b55));


	POINT3D du = f * (-10 * (b00 + 5 * b01 + 10 * b02 + 10 * b03 + 5 * b04 + b05) - 30 * (b10 + 5 * b11 + 10 * b12 + 10 * b13 + 5 * b14 + b15) - 20 * (b20 + 5 * b21 + 10 * b22 + 10 * b23 + 5 * b24 + b25) + 20 * (b30 + 5 * b31 + 10 * b32 + 10 * b33 + 5 * b34 + b35) + 30 * (b40 + 5 * b41 + 10 * b42 + 10 * b43 + 5 * b44 + b45) + 10 * (b50 + 5 * b51 + 10 * b52 + 10 * b53 + 5 * b54 + b55));

	POINT3D dv = f * (-10 * b00 - 30 * b01 - 20 * b02 + 20 * b03 + 30 * b04 + 10 * b05 + 5 * (-10 * b10 - 30 * b11 - 20 * b12 + 20 * b13 + 30 * b14 + 10 * b15) + 10 * (-10 * b20 - 30 * b21 - 20 * b22 + 20 * b23 + 30 * b24 + 10 * b25) + 10 * (-10 * b30 - 30 * b31 - 20 * b32 + 20 * b33 + 30 * b34 + 10 * b35) + 5 * (-10 * b40 - 30 * b41 - 20 * b42 + 20 * b43 + 30 * b44 + 10 * b45) + (-10 * b50 - 30 * b51 - 20 * b52 + 20 * b53 + 30 * b54 + 10 * b55));

	resn = du % dv; 
	resn = -resn.normalize();

	if (((resp-p1)| (resp - p1))>1) {
		cout << "------------" << endl;
		cout << "p1:" << p1 << endl;
		cout << "p2:" << p2 << endl;
		cout << "p3:" << p3 << endl;
		cout << "p4:" << p4 << endl;
		cout << "nor1:" << nor1 << endl;
		cout << "nor2:" << nor2 << endl;
		cout << "nor3:" << nor3 << endl;
		cout << "nor4:" << nor4 << endl;
		cout << "resp:" << resp << endl;
		cout << "resn:" << resn << endl;
	}
}

void create_new_T(POINT3D P1, POINT3D P2, POINT3D N1, POINT3D N2, POINT3D &newT)
{
	POINT3D T0, T1, V, B1, B2, F0, F1, F;
	V = P2 - P1;
	T0 = V - (N1 | V)*N1;
	T1 = V - (N2 | V)*N2;
	double cos0, cos1, k0, k1, sin0, sin1;
	cos0 = (T0 | V) / (T0.norm()*V.norm());
	cos1 = (T1 | V) / (T1.norm()*V.norm());
	sin0 = sqrt(1 - cos0 * cos0);
	sin1 = sqrt(1 - cos1 * cos1);
	k0 = 2.0 / (3 * cos0*(cos0 + 1));
	k1 = 2.0 / (3 * cos1*(cos1 + 1));
	//if (abs(abs(cos1) - abs(cos0)) > 0.95)
	//{
	//	k0 = 2 * sin1 * k0 / (sin0 + sin1);
	//	k1 = 2 * sin0 * k1 / (sin0 + sin1);
	//}

	B1 = P1 + k0 * T0;
	B2 = P2 - k1 * T1;

	newT = 3.0 / 4 * (B2 + P2 - P1 - B1);
}

void mao_method(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn) {

	if (1) {
		POINT3D p12, p23, p34, p41;
		POINT3D nor12, nor23, nor34, nor41;
		subdiv sub;
		sub.create_new_point_and_nor(p1, p2, nor1, nor2, p12, nor12);
		sub.create_new_point_and_nor(p2, p3, nor2, nor3, p23, nor23);
		sub.create_new_point_and_nor(p3, p4, nor3, nor4, p34, nor34);
		sub.create_new_point_and_nor(p4, p1, nor4, nor1, p41, nor41);

		POINT3D M0, M1, Mn0, Mn1;
		sub.create_new_point_and_nor(p12, p34, nor12, nor34, M0, Mn0);
		sub.create_new_point_and_nor(p23, p41, nor23, nor41, M1, Mn1);

		resp = (M0 + M1) / 2;

		POINT3D newT1, newT2, newN;
		create_new_T(p12, p34, nor12, nor34, newT1);
		create_new_T(p23, p41, nor23, nor41, newT2);
		resn = (newT1%newT2).normalize();
	
	}
	else {
		POINT3D p13, p24;
		POINT3D nor13, nor24;
		subdiv sub;
		sub.create_new_point_and_nor(p1, p3, nor1, nor3, p13, nor13);
		sub.create_new_point_and_nor(p2, p4, nor2, nor4, p24, nor24);
		
		resp = (p13 + p24) / 2;
		
		POINT3D newT1, newT2, newN;
		create_new_T(p1, p3, nor1, nor3, newT1);
		create_new_T(p2, p4, nor2, nor4, newT2);
		resn = (newT1%newT2).normalize();
	}


	//sub.create_new_point_and_nor(p1, p2, nor1, nor2, p12, nor12);
	//sub.create_new_point_and_nor(p2, p3, nor2, nor3, p23, nor23);

}


void fit_mao_bezier_quad(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn) {
	POINT3D p00, p01, p02, p03,
			p10, p11, p12, p13,
			p20, p21, p22, p23,
			p30, p31, p32, p33;
	// 初始化角点
	p00 = p1;
	p30 = p2;
	p33 = p3;
	p03 = p4;

	get_bezier_ctrlpoint(p00, p30, nor1, nor2, p10, p20);
	get_bezier_ctrlpoint(p30, p33, nor2, nor3, p31, p32);
	get_bezier_ctrlpoint(p33, p03, nor3, nor4, p23, p13);
	get_bezier_ctrlpoint(p03, p00, nor4, nor1, p02, p01);

	POINT3D m12, m23, m34, m41;
	POINT3D nor12, nor23, nor34, nor41;
	subdiv sub;
	sub.create_new_point_and_nor(p1, p2, nor1, nor2, m12, nor12);
	sub.create_new_point_and_nor(p2, p3, nor2, nor3, m23, nor23);
	sub.create_new_point_and_nor(p3, p4, nor3, nor4, m34, nor34);
	sub.create_new_point_and_nor(p4, p1, nor4, nor1, m41, nor41);


	Eigen::MatrixXd A(4, 4), h, b(4, 1);
	POINT3D tp11, tp12, tp21, tp22;
	tp11 = 4.0 / 9 * p1 + 2.0 / 9 * p2 + 2.0 / 9 * p4 + 1.0 / 9 * p3;
	tp12 = 4.0 / 9 * p4 + 2.0 / 9 * p1 + 2.0 / 9 * p3 + 1.0 / 9 * p2;
	tp21 = 4.0 / 9 * p2 + 2.0 / 9 * p3 + 2.0 / 9 * p1 + 1.0 / 9 * p4;
	tp22 = 4.0 / 9 * p3 + 2.0 / 9 * p2 + 2.0 / 9 * p4 + 1.0 / 9 * p1;

	POINT3D n11, n12, n21, n22;
	n11 = ((p1 - p4) % (p2 - p1)).normalize();
	n12 = ((p4 - p3) % (p1 - p4)).normalize();
	n21 = ((p2 - p1) % (p3 - p2)).normalize();
	n22 = ((p3 - p2) % (p4 - p3)).normalize();

	A << 3 * (n11 | nor12), 3 * (n21 | nor12), 0, 0,
		0, -3 * (n21 | nor23), -3 * (n22 | nor23), 0,
		0, 0, -3 * (n22 | nor34), -3 * (n12 | nor34),
		3 * (n11 | nor41), 0, 0, 3 * (n12 | nor41);

	b << ((3 * p00 - p01 + 9 * p10 - 3 * tp11 + 9 * p20 - 3 * tp21 + 3 * p30 - p31) | nor12),
		((p20 - 3 * p30 + 3 * tp21 - 9 * p31 + 3 * tp22 - 9 * p32 + p23 - 3 * p33) | nor23),
		((p02 - 3 * p03 + 3 * tp12 - 9 * p13 + 3 * tp22 - 9 * p23 + p32 - 3 * p33) | nor34),
		((3 * p00 - p10 + 9 * p01 - 3 * tp11 + 9 * p02 - 3 * tp12 + 3 * p03 - p13) | nor41);

	cout << b;
	cout << A;
	//Eigen::HouseholderQR<Eigen::MatrixXd> qr;
	//qr.compute(A);
	//Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
	//Eigen::MatrixXd Q = qr.householderQ();

	//Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(R);
	//int rank = lu_decomp.rank();

	h = (A.inverse())*b;
	p11 = tp11 + h(1, 1)*n11;
	p12 = tp12 + h(2, 1)*n12;
	p21 = tp21 + h(3, 1)*n21;
	p22 = tp22 + h(4, 1)*n22;

	resp = 1.0 / 64 * (p00 + 3 * p01 + 3 * p02 + p03 + 3 * (p10 + 3 * p11 + 3 * p12 + p13) + 3 * (p20 + 3 * p21 + 3 * p22 + p23) + p30 + 3 * p31 + 3 * p32 + p33);

	POINT3D derivative_u = 3.0 / 32 * (-(p00 + 3 * p01 + 3 * p02 + p03) - 3 * (p10 + 3 * p11 + 3 * p12 + p13) + 3 * (p20 + 3 * p21 + 3 * p22 + p23) + p30 + 3 * p31 + 3 * p32 + p33);

	POINT3D derivative_v = 3.0 / 32 * ((-p00 - 3 * p01 + 3 * p02 + p03) + 3 * (-p10 - 3 * p11 + 3 * p12 + p13) + 3 * (-p20 - 3 * p21 + 3 * p22 + p23) + (-p30 - 3 * p31 + 3 * p32 + p33));

	resn = -(derivative_u.cross(derivative_v));
	resn /= resn.norm();
}



void cac_point_and_nor_on_quad_bezier_surface(POINT3D p1, POINT3D p2, POINT3D p3, POINT3D p4, POINT3D nor1, POINT3D nor2, POINT3D nor3, POINT3D nor4, float u, float v, POINT3D& resp, POINT3D& resn)
{
	// method 1 
	//resp = (p1 + p2 + p3 + p4) / 4;
	//resn = (nor1 + nor2 + nor3 + nor4) / 4;
	//resn.normalize();

	//method 2:bicubic
	//get_quad_bezier_suface_ctrlpoint_and_cac(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);

	//method 3:classic PN quads
	//get_quad_bezier_suface_ctrlpoint_and_cac_PN_quads(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);
	//resp = (p1 + p2 + p3 + p4) / 4;

	//method 4 两次Bezier
	//get_quad_bezier_suface_ctrlpoint_and_cac_degree2(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);

	//method: local PN
	//get_quad_bezier_ctrlpoint_local_method(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);
	//cout << "p:"<<resp << ",n:"<< resn << endl;

	//method: Mao method
	mao_method(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);

	//fit_mao_bezier_quad(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);

	//minimiaze digonal energy
	//minimize_diagonal_energy(p1, p2, p3, p4, nor1, nor2, nor3, nor4, 0.5, 0.5, resp, resn);

}


subdiv::subdiv()
{

}

subdiv::~subdiv()
{

}

int subdiv::getindex(PolygonMesh::VertexHandle vh)
{
	for (int i = 0; i < vhandle.size(); ++i)
		if (vhandle[i] == vh)
			return i;

	return -1;
}


void subdiv::create_new_point_and_nor(POINT3D P1, POINT3D P2, POINT3D N1, POINT3D N2, POINT3D &newP, POINT3D &newN)
{
	POINT3D T0, T1, V, newT, B1, B2, F0, F1, F;
	V = P2 - P1;
	T0 = V - (N1 | V)*N1;
	T1 = V - (N2 | V)*N2;
	double cos0, cos1, k0, k1, sin0, sin1;
	cos0 = (T0 | V) / (T0.norm()*V.norm());
	cos1 = (T1 | V) / (T1.norm()*V.norm());
	sin0 = sqrt(1 - cos0 * cos0);
	sin1 = sqrt(1 - cos1 * cos1);
	k0 = 2.0 / (3 * cos0*(cos0 + 1));
	k1 = 2.0 / (3 * cos1*(cos1 + 1));
	if (abs(abs(cos1) - abs(cos0)) > 0.95)
	{
		k0 = 2 * sin1 * k0 / (sin0 + sin1);
		k1 = 2 * sin0 * k1 / (sin0 + sin1);
	}

	B1 = P1 + k0 * T0;
	B2 = P2 - k1 * T1;

	newP = 1.0 / 8 * (P1 + 3 * B1 + 3 * B2 + P2);
	newT = 3.0 / 4 * (B2 + P2 - P1 - B1);

	F0 = (V % N1) / (V % N1).norm();
	F1 = (V % N2) / (V % N2).norm();
	if ((F0 + F1).norm() < 0.00001 && (F0 + F1).norm() > -0.00001)
	{
		cout << "error---------" << endl;
	}
	F = (F0 + F1) / (F0 + F1).norm();
	newN = -(newT % F) / (newT % F).norm();
}

void subdiv::subdivsion(vector<POINT3D>points, vector<POINT3D> nors, vector<POINT3D>&resp, vector<POINT3D>&resn)
{
	resp.clear();
	resn.clear();
	int n = points.size();
	int i = 0;
	vector<POINT3D> res;
	POINT3D T0, T1, V, newP, newT, newNor, B1, B2, F0, F1, F;
	double cos0, cos1, k0, k1;
	for (i = 0; i < n - 1; ++i)
	{
		create_new_point_and_nor(points[i], points[i + 1], nors[i], nors[i + 1], newP, newNor);

		cout << "newP:" << " " << newP[0] << " " << newP[1] << " " << newP[2] << " " << endl;
		cout << "newNor:" << " " << newNor[0] << " " << newNor[1] << " " << newNor[2] << " " << endl;

		resp.push_back(points[i]);
		resp.push_back(newP);
		resn.push_back(nors[i]);
		resn.push_back(newNor);
	}
	if (1)//非 闭
	{
		resp.push_back(points[i]);
		resn.push_back(nors[i]);
	}
	else
	{
		create_new_point_and_nor(points[0], points[n - 1], nors[0], nors[n - 1], newP, newNor);

		resp.push_back(points[n - 1]);
		resp.push_back(newP);
		resn.push_back(nors[n - 1]);
		resn.push_back(newNor);

	}
}


//wrong 重复加点
//void subdiv::sub_for_mesh(PolygonMesh* mesh, PolygonMesh &subvMesh)
//{
//
//	subvMesh.clear();
//
//	POINT3D from, to, from_nor, to_nor, newP, newN;
//	int id_from, id_to;
//	//vector<PolygonMesh::VertexHandle> vhandle;
//	vhandle.clear();
//	int i = 0;
//	vector<int> tmp_id;  //用于存储一个三角面的三条边的from和to index。
//	set<id_pair> id_set;
//	//map<int, int> id2newVhandle;
//	//int myMap[50][50];
//	//memset(myMap, 0, 2500);
//
//	vector<PolygonMesh::VertexHandle> face_vhandles;
//	PolygonMesh::VertexHandle tmp_vhandle;
//	PolygonMesh::Normal tmp_nor;
//
//	subvMesh.request_vertex_normals();
//	subvMesh.update_normals();
//
//	//////以为重复加点，所以需要记录已经处理过的边，根据from_id，存放到一个set中或用vector（n） flag；
//
//
//	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
//	{
//		tmp_id.clear();
//		vhandle.clear();
//		for (PolygonMesh::FaceEdgeIter fe_it = mesh->fe_iter(*f_it); fe_it.is_valid(); ++fe_it)
//		{
//			from = mesh->point(mesh->from_vertex_handle(mesh->halfedge_handle(fe_it, 0)));
//			id_from = (mesh->from_vertex_handle(mesh->halfedge_handle(fe_it, 0))).idx();		
//			from_nor = mesh->normal(mesh->from_vertex_handle(mesh->halfedge_handle(fe_it, 0)));
//
//			to = mesh->point(mesh->to_vertex_handle(mesh->halfedge_handle(fe_it, 0)));
//			id_to = (mesh->to_vertex_handle(mesh->halfedge_handle(fe_it, 0))).idx();
//			to_nor = mesh->normal(mesh->to_vertex_handle(mesh->halfedge_handle(fe_it, 0)));
//
//			//tmp_id.push_back(id_from);
//			//tmp_id.push_back(id_to);
//
//
//			create_new_point_and_nor(from, to, from_nor, to_nor, newP, newN);
//
//			id_pair tmp_id_pair;
//			tmp_id_pair.id1 = id_from;
//			tmp_id_pair.id2 = id_to;
//
//			if (id_set.count(tmp_id_pair) != 1) 
//			{
//				id_set.insert(tmp_id_pair);
//				tmp_id_pair.id1 = id_to;
//				tmp_id_pair.id2 = id_from;
//				id_set.insert(tmp_id_pair);
//
//				tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(from[0], from[1], from[2]));
//				tmp_nor[0] = from_nor[0];
//				tmp_nor[1] = from_nor[1];
//				tmp_nor[2] = from_nor[2];
//				subvMesh.set_normal(tmp_vhandle, tmp_nor);
//				vhandle.push_back(tmp_vhandle);
//
//				tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newP[0], newP[1], newP[2]));
//				tmp_nor[0] = newN[0];
//				tmp_nor[1] = newN[1];
//				tmp_nor[2] = newN[2];
//				subvMesh.set_normal(tmp_vhandle, tmp_nor);
//				vhandle.push_back(tmp_vhandle);
//			
//			}
//			else
//			{
//				tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(to[0], to[1], to[2]));//tmpMesh只用于使用add_vertex函数来创建vertexhandle;
//				tmp_nor[0] = to_nor[0];
//				tmp_nor[1] = to_nor[1];
//				tmp_nor[2] = to_nor[2];
//				subvMesh.set_normal(tmp_vhandle, tmp_nor);
//				vhandle.push_back(tmp_vhandle);
//
//				tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newP[0], newP[1], newP[2]));
//				tmp_nor[0] = newN[0];
//				tmp_nor[1] = newN[1];
//				tmp_nor[2] = newN[2];
//				subvMesh.set_normal(tmp_vhandle, tmp_nor);
//				vhandle.push_back(tmp_vhandle);
//
//			}
//
//		}
//
//		face_vhandles.clear();
//		face_vhandles.push_back(vhandle[1]);
//		face_vhandles.push_back(vhandle[3]);
//		face_vhandles.push_back(vhandle[5]);
//		subvMesh.add_face(face_vhandles);
//		face_vhandles.clear();
//		face_vhandles.push_back(vhandle[1]);
//		face_vhandles.push_back(vhandle[2]);
//		face_vhandles.push_back(vhandle[3]);
//		subvMesh.add_face(face_vhandles);
//		face_vhandles.clear();
//		face_vhandles.push_back(vhandle[3]);
//		face_vhandles.push_back(vhandle[4]);
//		face_vhandles.push_back(vhandle[5]);
//		subvMesh.add_face(face_vhandles);
//		face_vhandles.clear();
//		face_vhandles.push_back(vhandle[0]);
//		face_vhandles.push_back(vhandle[1]);
//		face_vhandles.push_back(vhandle[5]);
//		subvMesh.add_face(face_vhandles);
//	}
//
//}



void subdiv::sub_for_mesh_new(PolygonMesh* mesh, PolygonMesh &subvMesh)
{

	subvMesh.clear();

	POINT3D from, to, from_nor, to_nor, newP, newN;
	int id_from, id_to;
	vector<PolygonMesh::VertexHandle> origin_vhandle, add_vhandle;
	vhandle.clear();
	int i = 0;
	vector<int> tmp_id;  //用于存储一个三角面的三条边的from和to index。
	set<id_pair> id_set;

	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	//保存原先面的顶点
	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmpP;
		tmpP = mesh->point(*v_it);
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(tmpP[0], tmpP[1], tmpP[2]));
		tmp_nor = mesh->normal(*(v_it));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);
		origin_vhandle.push_back(tmp_vhandle);
	}

	//保存每条边新增顶点
	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++)
	{
		from = mesh->point(mesh->from_vertex_handle(mesh->halfedge_handle(e_it, 0)));
		from_nor = mesh->normal(mesh->from_vertex_handle(mesh->halfedge_handle(e_it, 0)));

		to = mesh->point(mesh->to_vertex_handle(mesh->halfedge_handle(e_it, 0)));
		to_nor = mesh->normal(mesh->to_vertex_handle(mesh->halfedge_handle(e_it, 0)));

		create_new_point_and_nor(from, to, from_nor, to_nor, newP, newN);

		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newP[0], newP[1], newP[2]));
		tmp_nor[0] = newN[0];
		tmp_nor[1] = newN[1];
		tmp_nor[2] = newN[2];
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		add_vhandle.push_back(tmp_vhandle);
	}

	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		vhandle.clear();
		for (PolygonMesh::FaceEdgeIter fe_it = mesh->fe_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			vhandle.push_back(add_vhandle[(*fe_it).idx()]);
		}

		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			vhandle.push_back(origin_vhandle[(*fv_it).idx()]);
		}

		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[2]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[3]);
		face_vhandles.push_back(vhandle[1]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[4]);
		face_vhandles.push_back(vhandle[2]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[5]);
		subvMesh.add_face(face_vhandles);
	
	
	}

}


void subdiv::sub_for_mesh_by_sqrt3(PolygonMesh* mesh, PolygonMesh &subvMesh)
{
	subvMesh.clear();

	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;
	vector<POINT3D>tmp_face_points_and_nors;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	vector<PolygonMesh::VertexHandle> midpoint_vec;
	vector<PolygonMesh::VertexHandle> origin_points_vec;
	vector<vector<int>> origin_face_points_id(faces_num);

	vector<POINT3D> my_point_nors_vec;

	for (PolygonMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_element;
		tmp_element = mesh->point(v_it);
		tmp_face_points_and_nors.push_back(tmp_element);

		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(tmp_element[0], tmp_element[1], tmp_element[2]));

		tmp_element = mesh->normal(v_it);
		tmp_face_points_and_nors.push_back(tmp_element);
		tmp_nor[0] = tmp_element[0];
		tmp_nor[1] = tmp_element[1];
		tmp_nor[2] = tmp_element[2];

		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		origin_points_vec.push_back(tmp_vhandle);
	}

	int ind = 0;
	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		
		vhandle.clear();
		tmp_face_points_and_nors.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			tmp_element = mesh->normal(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);

			origin_face_points_id[ind].push_back((*fv_it).idx());
		}

		//cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 0.1, 0.45, 0.45, newp, newNor);

		cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 1 / 3.0, 1 / 3.0, 1 / 3.0, newp, newNor);

		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		midpoint_vec.push_back(tmp_vhandle);
		ind++;
	}

	//add faces
	for (int i = 0; i < faces_num; ++i)
	{
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		subvMesh.add_face(face_vhandles);
	}
	
	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		
		if (!subvMesh.is_boundary(*e_it))
		{
				int id_from, id_to;
				id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
				id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
				if (id_from < vertices_num && id_to < vertices_num) //边不包括新生成的面中心点
				{
					subvMesh.flip(*e_it);
				}
		}
	}

	//for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	//{
	//	if (!subvMesh.is_boundary(*e_it))
	//	{
	//		if (find_extraordinary_triangle(e_it, subvMesh))
	//		{
	//			subvMesh.flip(*e_it);
	//		}
	//	}
	//}

#if 1  //翻折狭长三角形且与其相邻三角形近乎垂直的边
	for (PolygonMesh::FaceIter f_it = subvMesh.faces_begin(); f_it != subvMesh.faces_end(); ++f_it)
	{
		vector<PolygonMesh::VertexHandle> tmpPoints;
		vector<PolygonMesh::HalfedgeHandle> tmpEdges;
		for (PolygonMesh::FaceVertexIter fv_it = subvMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmpPoints.push_back(*fv_it);
			//POINT3D p;
			//p = subvMesh.point(*fv_it);
			//cout << p<<endl;
		}
		for (PolygonMesh::FaceHalfedgeIter fe_it = subvMesh.fh_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			tmpEdges.push_back(*fe_it);
			//POINT3D p;
			//p = subvMesh.point(subvMesh.to_vertex_handle(*fe_it));
			//cout << p << endl;
		}
		

		POINT3D p1, p2, p3;
		p1 = subvMesh.point(tmpPoints[0]); p2 = subvMesh.point(tmpPoints[1]); p3 = subvMesh.point(tmpPoints[2]);
		
		int face_flg = is_unregular_face(p1, p2, p3);
		switch (face_flg)
		{
		case 1:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[2]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[2]));
			}
			break;
		case 2:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[0]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[0]));
			}
			break;
		case 3:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[1]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[1]));
			}
			break;
		default:
			break;
		}

	}
#endif



	//boundary
# if 0
	subvMesh.request_face_status();
	subvMesh.request_edge_status();
	subvMesh.request_vertex_status();
	vector<vector<PolygonMesh::VertexHandle>>boundary_vec;
	vector<PolygonMesh::FaceHandle>need_to_delete_faceh;
	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		if (subvMesh.is_boundary(*e_it))
		{
			PolygonMesh::VertexHandle from, to, opposite;
			from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0)));
			to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0)));
			opposite = (subvMesh.opposite_vh(mesh->halfedge_handle(e_it, 0)));
			POINT3D p1, p2, n1, n2, p3, n3;
			p1 = subvMesh.point(from);
			p2 = subvMesh.point(to);
			n1 = subvMesh.normal(from);
			n2 = subvMesh.normal(to);
			
			create_new_point_and_nor(p1, p2, n1, n2, p3, n3);

			tmp_nor[0] = n3[0]; tmp_nor[1] = n3[1]; tmp_nor[2] = n3[2];
			tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(p3[0], p3[1], p3[2]));
			subvMesh.set_normal(tmp_vhandle, tmp_nor);


			PolygonMesh::FaceHandle fh = subvMesh.opposite_face_handle(mesh->halfedge_handle(e_it, 1));
			need_to_delete_faceh.push_back(fh);

			vector<PolygonMesh::VertexHandle> tmp;  //from to opposite new
			tmp.push_back(from);
			tmp.push_back(to);
			tmp.push_back(opposite);
			tmp.push_back(tmp_vhandle);
			boundary_vec.push_back(tmp);

			//delete face
			//subvMesh.delete_edge(e_it);
			//subvMesh.delete_face(fh, false);

		}
	}
	for (auto x : need_to_delete_faceh) {
		subvMesh.delete_face(x,false);
	}
	subvMesh.garbage_collection();

	//add boundary face;
	for (auto val : boundary_vec)
	{
		face_vhandles.clear();
		face_vhandles.push_back(val[0]);
		face_vhandles.push_back(val[3]);
		face_vhandles.push_back(val[2]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(val[3]);
		face_vhandles.push_back(val[1]);
		face_vhandles.push_back(val[2]);
		subvMesh.add_face(face_vhandles);
	}

	subvMesh.release_face_status();
	subvMesh.release_edge_status();
	subvMesh.release_vertex_status();

#endif

#if 1
	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		if (!subvMesh.is_boundary(*e_it))
		{
			cac_two_tri_bezier_surfaces_normals_continuity(e_it, subvMesh);
		}

	}
	cout << "max_two_bezier_triangle_cos: " << two_bezier_triangle_cos << endl;
	two_bezier_triangle_cos = 1.0;
	//cout << "error_cnt:" << error_normal_cnt << endl;
	//error_normal_cnt = 0;
#endif

#if 0
	subvMesh =  *mesh;
	int origin_size = subvMesh.n_vertices();
	
	subvMesh.request_face_status();
	subvMesh.request_edge_status();
	subvMesh.request_vertex_status();


	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;
	vector<POINT3D>tmp_face_points_and_nors;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	vector<PolygonMesh::VertexHandle> midpoint_vec;
	vector<PolygonMesh::VertexHandle> tmp_points_of_origin_face;

	vector<POINT3D> my_point_nors_vec;


	//先给新曲面增加所有面点并与面顶点相连
	for (PolygonMesh::FaceIter f_it = subvMesh.faces_begin(); f_it != subvMesh.faces_end(); ++f_it)
	{
		vhandle.clear();
		tmp_face_points_and_nors.clear();
		tmp_points_of_origin_face.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = subvMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			tmp_element = mesh->normal(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);

			tmp_points_of_origin_face.push_back(*fv_it);
		}

		cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 1 / 3.0, 1 / 3.0, 1 / 3.0, newp, newNor);

		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		midpoint_vec.push_back(tmp_vhandle);

		//add face
		face_vhandles.clear();
		face_vhandles.push_back(tmp_vhandle);
		face_vhandles.push_back(tmp_points_of_origin_face[1]);
		face_vhandles.push_back(tmp_points_of_origin_face[0]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(tmp_vhandle);
		face_vhandles.push_back(tmp_points_of_origin_face[2]);
		face_vhandles.push_back(tmp_points_of_origin_face[1]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(tmp_vhandle);
		face_vhandles.push_back(tmp_points_of_origin_face[0]);
		face_vhandles.push_back(tmp_points_of_origin_face[2]);
		subvMesh.add_face(face_vhandles);

		subvMesh.delete_face(*f_it);

	}

	//翻转原边（除边界边）
	//for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	//{
	//	if (!subvMesh.is_boundary(*e_it))
	//	{
	//		int id_from, id_to;
	//		id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
	//		id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
	//		if (id_from < origin_size && id_to < origin_size)  //边不包括新生成的面中心点
	//		{
	//			subvMesh.flip(*e_it);
	//		}
	//	}
	//}

#endif

# if 0
	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;

	set<int> midpoint_id_set;
	set<POINT3D> newMesh_point_set;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();
	vector<POINT3D>tmp_face_points_and_nors;

	
	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		tmp_face_points_and_nors.clear();
		vhandle.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			POINT3D tmp_point;
			tmp_point = mesh->point(fv_it);
			if (0 == newMesh_point_set.count(tmp_point))
			{
				newMesh_point_set.insert(tmp_point);
				tmp_vhandle = subvMesh.add_vertex(mesh->point(fv_it));
			}
			else
			{
				tmp_vhandle = *fv_it;
			}

			//tmp_vhandle = subvMesh.add_vertex(mesh->point(fv_it));
			
			vhandle.push_back(tmp_vhandle);
			
			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			tmp_element = mesh->normal(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			
			tmp_nor[0] = tmp_element[0]; tmp_nor[1] = tmp_element[1]; tmp_nor[2] = tmp_element[2];
			subvMesh.set_normal(tmp_vhandle, tmp_nor);
		}
		
		cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 1 / 3.0, 1 / 3.0, 1 / 3.0, newp, newNor);


		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);
		vhandle.push_back(tmp_vhandle);

		midpoint_id_set.insert(tmp_vhandle.idx());

		//添加面
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[3]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[3]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[3]);
		subvMesh.add_face(face_vhandles);

	}

	

	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		if (!subvMesh.is_boundary(*e_it))
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			if (midpoint_id_set.count(id_from) != 1 && midpoint_id_set.count(id_to) != 1)  //边不包括新生成的面中心点
			{
				//subvMesh.flip(*e_it);
			}
		}
		else
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			//printf("-----------id_from:%d,id_to:%d------------\n", id_from,id_to);
		}
	}

#endif

}

void subdiv::sub_for_mesh_by_sqrt3_PN_method(PolygonMesh* mesh, PolygonMesh &subvMesh)
{
	subvMesh.clear();

	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;
	vector<POINT3D>tmp_face_points_and_nors;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	vector<PolygonMesh::VertexHandle> midpoint_vec;
	vector<PolygonMesh::VertexHandle> origin_points_vec;
	vector<vector<int>> origin_face_points_id(faces_num);

	vector<POINT3D> my_point_nors_vec;

	for (PolygonMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_element;
		tmp_element = mesh->point(v_it);
		tmp_face_points_and_nors.push_back(tmp_element);

		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(tmp_element[0], tmp_element[1], tmp_element[2]));

		tmp_element = mesh->normal(v_it);
		tmp_face_points_and_nors.push_back(tmp_element);
		tmp_nor[0] = tmp_element[0];
		tmp_nor[1] = tmp_element[1];
		tmp_nor[2] = tmp_element[2];

		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		origin_points_vec.push_back(tmp_vhandle);
	}

	int ind = 0;
	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{

		vhandle.clear();
		tmp_face_points_and_nors.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			tmp_element = mesh->normal(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);

			origin_face_points_id[ind].push_back((*fv_it).idx());
		}

		//cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 0.1, 0.45, 0.45, newp, newNor);

		cac_point_and_nor_on_tri_bezier_surface_PN(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 1 / 3.0, 1 / 3.0, 1 / 3.0, newp, newNor);

		//newp = 1.0 / 3 * (tmp_face_points_and_nors[0] + tmp_face_points_and_nors[1] + tmp_face_points_and_nors[2]);
		//newNor = 1.0 / 3 * (tmp_face_points_and_nors[3] + tmp_face_points_and_nors[4] + tmp_face_points_and_nors[5]);



		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		midpoint_vec.push_back(tmp_vhandle);
		ind++;
	}

	//add faces
	for (int i = 0; i < faces_num; ++i)
	{
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		subvMesh.add_face(face_vhandles);
	}

	/***************************/
	PolygonMesh::Point final_point;
	for (int i = 0; i < midpoint_vec.size(); ++i)
	{
		POINT3D tmp_sum_nor(0, 0, 0);
		double sum_height = 0.0;
		POINT3D new_face_nor, new_face_point;
		new_face_point = subvMesh.point(midpoint_vec[i]);
		new_face_nor = subvMesh.normal(midpoint_vec[i]);
		for (auto vv_it = subvMesh.vv_iter(midpoint_vec[i]); vv_it.is_valid(); ++vv_it)
		{
			POINT3D new_origin_ver_nor, new_origin_ver;
			new_origin_ver = subvMesh.point(vv_it);
			new_origin_ver_nor = subvMesh.normal(vv_it);
			sum_height += ((new_face_nor + new_origin_ver_nor) | (new_origin_ver - new_face_point)) / ((new_face_nor + new_origin_ver_nor) | new_face_nor);
		}
		tmp_sum_nor = 1 / 6.0*sum_height*new_face_nor;
		new_face_point += tmp_sum_nor;
		final_point[0] = new_face_point[0]; final_point[1] = new_face_point[1]; final_point[2] = new_face_point[2];
		subvMesh.set_point(midpoint_vec[i], final_point);
	}



	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{

		if (!subvMesh.is_boundary(*e_it))
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			if (id_from < vertices_num && id_to < vertices_num) //边不包括新生成的面中心点
			{
				subvMesh.flip(*e_it);
			}
		}
	}

	//for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	//{
	//	if (!subvMesh.is_boundary(*e_it))
	//	{
	//		if (find_extraordinary_triangle(e_it, subvMesh))
	//		{
	//			subvMesh.flip(*e_it);
	//		}
	//	}
	//}

#if 1  //翻折狭长三角形且与其相邻三角形近乎垂直的边
	for (PolygonMesh::FaceIter f_it = subvMesh.faces_begin(); f_it != subvMesh.faces_end(); ++f_it)
	{
		vector<PolygonMesh::VertexHandle> tmpPoints;
		vector<PolygonMesh::HalfedgeHandle> tmpEdges;
		for (PolygonMesh::FaceVertexIter fv_it = subvMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmpPoints.push_back(*fv_it);
			//POINT3D p;
			//p = subvMesh.point(*fv_it);
			//cout << p<<endl;
		}
		for (PolygonMesh::FaceHalfedgeIter fe_it = subvMesh.fh_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			tmpEdges.push_back(*fe_it);
			//POINT3D p;
			//p = subvMesh.point(subvMesh.to_vertex_handle(*fe_it));
			//cout << p << endl;
		}


		POINT3D p1, p2, p3;
		p1 = subvMesh.point(tmpPoints[0]); p2 = subvMesh.point(tmpPoints[1]); p3 = subvMesh.point(tmpPoints[2]);

		int face_flg = is_unregular_face(p1, p2, p3);
		switch (face_flg)
		{
		case 1:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[2]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[2]));
			}
			break;
		case 2:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[0]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[0]));
			}
			break;
		case 3:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[1]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[1]));
			}
			break;
		default:
			break;
		}

	}
#endif

#if 0
	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		if (!subvMesh.is_boundary(*e_it))
		{
			cac_two_tri_bezier_surfaces_normals_continuity(e_it, subvMesh);
		}

	}
	//cout << "max_two_bezier_triangle_cos: " << two_bezier_triangle_cos << endl;
	//cout << "error_cnt:" << error_normal_cnt << endl;
	error_normal_cnt = 0;
#endif

	//boundary
# if 0
	subvMesh.request_face_status();
	subvMesh.request_edge_status();
	subvMesh.request_vertex_status();
	vector<vector<PolygonMesh::VertexHandle>>boundary_vec;
	vector<PolygonMesh::FaceHandle>need_to_delete_faceh;
	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		if (subvMesh.is_boundary(*e_it))
		{
			cout << "------------" << endl;
			PolygonMesh::VertexHandle from, to, opposite;
			from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0)));
			to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0)));
			opposite = (subvMesh.opposite_vh(mesh->halfedge_handle(e_it, 0)));
			POINT3D p1, p2, n1, n2, p3, n3;
			p1 = subvMesh.point(from);
			p2 = subvMesh.point(to);
			n1 = subvMesh.normal(from);
			n2 = subvMesh.normal(to);

			create_new_point_and_nor(p1, p2, n1, n2, p3, n3);

			tmp_nor[0] = n3[0]; tmp_nor[1] = n3[1]; tmp_nor[2] = n3[2];
			tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(p3[0], p3[1], p3[2]));
			subvMesh.set_normal(tmp_vhandle, tmp_nor);


			PolygonMesh::FaceHandle fh = subvMesh.opposite_face_handle(mesh->halfedge_handle(e_it, 1));
			need_to_delete_faceh.push_back(fh);

			vector<PolygonMesh::VertexHandle> tmp;  //from to opposite new
			tmp.push_back(from);
			tmp.push_back(to);
			tmp.push_back(opposite);
			tmp.push_back(tmp_vhandle);
			boundary_vec.push_back(tmp);

			//delete face
			//subvMesh.delete_edge(e_it);
			//subvMesh.delete_face(fh, false);

		}
	}
	for (auto x : need_to_delete_faceh) {
		subvMesh.delete_face(x, false);
	}
	subvMesh.garbage_collection();

	//add boundary face;
	for (auto val : boundary_vec)
	{
		face_vhandles.clear();
		face_vhandles.push_back(val[0]);
		face_vhandles.push_back(val[3]);
		face_vhandles.push_back(val[2]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(val[3]);
		face_vhandles.push_back(val[1]);
		face_vhandles.push_back(val[2]);
		subvMesh.add_face(face_vhandles);
	}

	subvMesh.release_face_status();
	subvMesh.release_edge_status();
	subvMesh.release_vertex_status();

#endif



#if 0
	subvMesh = *mesh;
	int origin_size = subvMesh.n_vertices();

	subvMesh.request_face_status();
	subvMesh.request_edge_status();
	subvMesh.request_vertex_status();


	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;
	vector<POINT3D>tmp_face_points_and_nors;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	vector<PolygonMesh::VertexHandle> midpoint_vec;
	vector<PolygonMesh::VertexHandle> tmp_points_of_origin_face;

	vector<POINT3D> my_point_nors_vec;


	//先给新曲面增加所有面点并与面顶点相连
	for (PolygonMesh::FaceIter f_it = subvMesh.faces_begin(); f_it != subvMesh.faces_end(); ++f_it)
	{
		vhandle.clear();
		tmp_face_points_and_nors.clear();
		tmp_points_of_origin_face.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = subvMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			tmp_element = mesh->normal(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);

			tmp_points_of_origin_face.push_back(*fv_it);
		}

		cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 1 / 3.0, 1 / 3.0, 1 / 3.0, newp, newNor);

		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		midpoint_vec.push_back(tmp_vhandle);

		//add face
		face_vhandles.clear();
		face_vhandles.push_back(tmp_vhandle);
		face_vhandles.push_back(tmp_points_of_origin_face[1]);
		face_vhandles.push_back(tmp_points_of_origin_face[0]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(tmp_vhandle);
		face_vhandles.push_back(tmp_points_of_origin_face[2]);
		face_vhandles.push_back(tmp_points_of_origin_face[1]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(tmp_vhandle);
		face_vhandles.push_back(tmp_points_of_origin_face[0]);
		face_vhandles.push_back(tmp_points_of_origin_face[2]);
		subvMesh.add_face(face_vhandles);

		subvMesh.delete_face(*f_it);

	}

	//翻转原边（除边界边）
	//for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	//{
	//	if (!subvMesh.is_boundary(*e_it))
	//	{
	//		int id_from, id_to;
	//		id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
	//		id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
	//		if (id_from < origin_size && id_to < origin_size)  //边不包括新生成的面中心点
	//		{
	//			subvMesh.flip(*e_it);
	//		}
	//	}
	//}

#endif

# if 0
	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;

	set<int> midpoint_id_set;
	set<POINT3D> newMesh_point_set;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();
	vector<POINT3D>tmp_face_points_and_nors;


	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		tmp_face_points_and_nors.clear();
		vhandle.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			POINT3D tmp_point;
			tmp_point = mesh->point(fv_it);
			if (0 == newMesh_point_set.count(tmp_point))
			{
				newMesh_point_set.insert(tmp_point);
				tmp_vhandle = subvMesh.add_vertex(mesh->point(fv_it));
			}
			else
			{
				tmp_vhandle = *fv_it;
			}

			//tmp_vhandle = subvMesh.add_vertex(mesh->point(fv_it));

			vhandle.push_back(tmp_vhandle);

			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);
			tmp_element = mesh->normal(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);

			tmp_nor[0] = tmp_element[0]; tmp_nor[1] = tmp_element[1]; tmp_nor[2] = tmp_element[2];
			subvMesh.set_normal(tmp_vhandle, tmp_nor);
		}

		cac_point_and_nor_on_tri_bezier_surface(tmp_face_points_and_nors[0], tmp_face_points_and_nors[2], tmp_face_points_and_nors[4], tmp_face_points_and_nors[1], tmp_face_points_and_nors[3], tmp_face_points_and_nors[5], 1 / 3.0, 1 / 3.0, 1 / 3.0, newp, newNor);


		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);
		vhandle.push_back(tmp_vhandle);

		midpoint_id_set.insert(tmp_vhandle.idx());

		//添加面
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[3]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[3]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[3]);
		subvMesh.add_face(face_vhandles);

	}



	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{
		if (!subvMesh.is_boundary(*e_it))
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			if (midpoint_id_set.count(id_from) != 1 && midpoint_id_set.count(id_to) != 1)  //边不包括新生成的面中心点
			{
				//subvMesh.flip(*e_it);
			}
		}
		else
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			//printf("-----------id_from:%d,id_to:%d------------\n", id_from,id_to);
		}
	}

#endif

}

//don't have boundary
void subdiv::sub_for_mesh_by_sqrt2(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh)
{
	subvMesh.clear();
	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<T_PolygonMesh::VertexHandle>old_vertices(vertices_num);
	vector<T_PolygonMesh::VertexHandle>new_face_vertices(faces_num);
	T_PolygonMesh::VertexHandle tmp_add_vh;
	T_PolygonMesh::Normal tmp_nor;
	vector<T_PolygonMesh::VertexHandle>faces_handle;

	if (0)
	{
		preNormalDisjustion(mesh);
		//cnttt++;
		//char t[1];
		//t[0] = 'a';
		//if (!OpenMesh::IO::write_mesh(mesh, t, OpenMesh::IO::Options::VertexNormal))
		//{
		//	printf("write mesh file is error!\n");
		//}
	}
	

	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_point,tmp_norr;
		tmp_point = mesh->point(v_it);
		tmp_add_vh = subvMesh.add_vertex(T_PolygonMesh::Point(tmp_point[0], tmp_point[1], tmp_point[2]));
		tmp_norr = mesh->normal(v_it);
		tmp_nor[0] = tmp_norr[0]; tmp_nor[1] = tmp_norr[1]; tmp_nor[2] = tmp_norr[2];
		subvMesh.set_normal(tmp_add_vh, tmp_nor);
		old_vertices[(*v_it).idx()] = tmp_add_vh;
	}


	//cac_new_face_points
	for (auto f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		vector<POINT3D>face_points;
		vector<POINT3D>face_nors;
		for (auto fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			POINT3D tmp_point,tmp_norr;
			tmp_point = mesh->point(fv_it);
			face_points.push_back(tmp_point);
			tmp_norr = mesh->normal(fv_it);
			face_nors.push_back(tmp_norr);
		}

		POINT3D new_p, new_nor;
		cac_point_and_nor_on_quad_bezier_surface(face_points[0], face_points[1], face_points[2], face_points[3], face_nors[0], face_nors[1], face_nors[2], face_nors[3], 0.5, 0.5, new_p, new_nor);

		tmp_add_vh = subvMesh.add_vertex(T_PolygonMesh::Point(new_p[0], new_p[1], new_p[2]));
		tmp_nor[0] = new_nor[0]; tmp_nor[1] = new_nor[1]; tmp_nor[2] = new_nor[2];
		subvMesh.set_normal(tmp_add_vh, tmp_nor);

		new_face_vertices[(*f_it).idx()] = tmp_add_vh;
	}

	//POINT3D from, to, from_nor, to_nor;
	T_PolygonMesh::HalfedgeHandle tmp_heh;
	T_PolygonMesh::VertexHandle tmp_vh;
	T_PolygonMesh::FaceHandle tmp_fh;
	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		tmp_heh = mesh->halfedge_handle(e_it, 0);
		tmp_vh = mesh->from_vertex_handle(tmp_heh);
		tmp_fh = mesh->opposite_face_handle(tmp_heh);

		faces_handle.push_back(old_vertices[tmp_vh.idx()]);
		faces_handle.push_back(new_face_vertices[tmp_fh.idx()]);
			
		tmp_heh = mesh->halfedge_handle(e_it, 1);
		tmp_vh = mesh->from_vertex_handle(tmp_heh);
		tmp_fh = mesh->opposite_face_handle(tmp_heh);

		faces_handle.push_back(old_vertices[tmp_vh.idx()]);
		faces_handle.push_back(new_face_vertices[tmp_fh.idx()]);

		subvMesh.add_face(faces_handle);
		faces_handle.clear();
	}


#if 0

	//更新原顶点法向 
	//for (int i = 0; i < vertices_num; ++i)
	//{
	//	POINT3D tmp_sum_nor(0,0,0);
	//	for (auto vv_it = subvMesh.vv_iter(old_vertices[i]); vv_it.is_valid(); ++vv_it)
	//	{
	//		POINT3D vv_nor;
	//		vv_nor = subvMesh.normal(vv_it);
	//		tmp_sum_nor += vv_nor;
	//	}
	//	tmp_sum_nor = tmp_sum_nor / tmp_sum_nor.norm();
	//	tmp_nor[0] = tmp_sum_nor[0]; tmp_nor[1] = tmp_sum_nor[1]; tmp_nor[2] = tmp_sum_nor[2];
	//	subvMesh.set_normal(old_vertices[i], tmp_nor);
	//}

	//计算新面点在法向上的衍生新点
	T_PolygonMesh::Point final_point;
	for (int i = 0; i < new_face_vertices.size(); ++i)
	{
		POINT3D tmp_sum_nor(0, 0, 0);
		double sum_height = 0.0;
		POINT3D new_face_nor, new_face_point;
		new_face_point = subvMesh.point(new_face_vertices[i]);
		new_face_nor = subvMesh.normal(new_face_vertices[i]);
		int test_cnt = 0;
		for (auto vv_it = subvMesh.vv_iter(new_face_vertices[i]); vv_it.is_valid(); ++vv_it)
		{
			POINT3D new_origin_ver_nor, new_origin_ver;
			new_origin_ver = subvMesh.point(vv_it);
			new_origin_ver_nor = subvMesh.normal(vv_it);
			sum_height += ((new_face_nor + new_origin_ver_nor) | (new_origin_ver - new_face_point)) / ((new_face_nor + new_origin_ver_nor) | new_face_nor);
			test_cnt++;
		}
		if (test_cnt > 4) {
			cout << "cnt>4" << endl;
		}

		//最基础
		tmp_sum_nor = 1 / 4.0*sum_height*new_face_nor;


		new_face_point += tmp_sum_nor;
		final_point[0] = new_face_point[0]; final_point[1] = new_face_point[1]; final_point[2] = new_face_point[2];
		subvMesh.set_point(new_face_vertices[i], final_point);
	}
#endif
}


void subdiv::sqrt3(PolygonMesh* mesh, PolygonMesh &subvMesh) {
	subvMesh.clear();

	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;
	vector<POINT3D>tmp_face_points_and_nors;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	vector<PolygonMesh::VertexHandle> midpoint_vec;
	vector<PolygonMesh::VertexHandle> origin_points_vec;
	vector<vector<int>> origin_face_points_id(faces_num);

	vector<POINT3D> my_point_nors_vec;

	for (PolygonMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_element;
		tmp_element = mesh->point(v_it);
		tmp_face_points_and_nors.push_back(tmp_element);

		//updata v
		int deg = 0;
		POINT3D surround_sum(0, 0, 0);
		for (PolygonMesh::VertexVertexIter vv_it = mesh->vv_iter(v_it); vv_it.is_valid(); vv_it++) {
			deg++;
			surround_sum += mesh->point(vv_it);
		}
		POINT3D newp;
		double alpha = (4 - 2 * cos(2 * PI / deg)) / 9;
		newp = (1 - alpha)*tmp_element + alpha / deg * surround_sum;

		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		origin_points_vec.push_back(tmp_vhandle);
	}

	int ind = 0;
	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		vhandle.clear();
		tmp_face_points_and_nors.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmp_element = mesh->point(fv_it);
			tmp_face_points_and_nors.push_back(tmp_element);

			origin_face_points_id[ind].push_back((*fv_it).idx());
		}

		newp = (tmp_face_points_and_nors[0] + tmp_face_points_and_nors[1] + tmp_face_points_and_nors[2]) / 3.0;
		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		midpoint_vec.push_back(tmp_vhandle);
		ind++;
	}

	//add faces
	for (int i = 0; i < faces_num; ++i)
	{
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		subvMesh.add_face(face_vhandles);
	}

	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{

		if (!subvMesh.is_boundary(*e_it))
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			if (id_from < vertices_num && id_to < vertices_num) //边不包括新生成的面中心点
			{
				subvMesh.flip(*e_it);
			}
		}
	}


#if 0  //翻折狭长三角形且与其相邻三角形近乎垂直的边
	for (PolygonMesh::FaceIter f_it = subvMesh.faces_begin(); f_it != subvMesh.faces_end(); ++f_it)
	{
		vector<PolygonMesh::VertexHandle> tmpPoints;
		vector<PolygonMesh::HalfedgeHandle> tmpEdges;
		for (PolygonMesh::FaceVertexIter fv_it = subvMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmpPoints.push_back(*fv_it);
			//POINT3D p;
			//p = subvMesh.point(*fv_it);
			//cout << p<<endl;
		}
		for (PolygonMesh::FaceHalfedgeIter fe_it = subvMesh.fh_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			tmpEdges.push_back(*fe_it);
			//POINT3D p;
			//p = subvMesh.point(subvMesh.to_vertex_handle(*fe_it));
			//cout << p << endl;
		}


		POINT3D p1, p2, p3;
		p1 = subvMesh.point(tmpPoints[0]); p2 = subvMesh.point(tmpPoints[1]); p3 = subvMesh.point(tmpPoints[2]);

		int face_flg = is_unregular_face(p1, p2, p3);
		switch (face_flg)
		{
		case 1:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[2]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[2]));
			}
			break;
		case 2:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[0]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[0]));
			}
			break;
		case 3:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[1]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[1]));
			}
			break;
		default:
			break;
		}

	}
#endif
}


void subdiv::sqrt2(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh) {
	subvMesh.clear();
	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<T_PolygonMesh::VertexHandle>old_vertices(vertices_num);
	vector<T_PolygonMesh::VertexHandle>new_face_vertices(faces_num);
	T_PolygonMesh::VertexHandle tmp_add_vh;
	T_PolygonMesh::Normal tmp_nor;
	vector<T_PolygonMesh::VertexHandle>faces_handle;

	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_point, tmp_norr;
		tmp_point = mesh->point(v_it);
	
		//updata v
		int deg = 0;
		POINT3D surround_sum(0, 0, 0);
		for (PolygonMesh::VertexVertexIter vv_it = mesh->vv_iter(v_it); vv_it.is_valid(); vv_it++) {
			deg++;
			surround_sum += mesh->point(vv_it);
		}
		POINT3D newp;
		double alpha = (1 - cos(2 * PI / deg)) / 2;
		newp = (1 - alpha)*tmp_point + alpha / deg * surround_sum;

		tmp_add_vh = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_add_vh, tmp_nor);
		old_vertices[(*v_it).idx()] = tmp_add_vh;
	}


	//cac_new_face_points
	for (auto f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		vector<POINT3D>face_points;
		vector<POINT3D>face_nors;
		for (auto fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			POINT3D tmp_point, tmp_norr;
			tmp_point = mesh->point(fv_it);
			face_points.push_back(tmp_point);
			tmp_norr = mesh->normal(fv_it);
			face_nors.push_back(tmp_norr);
		}

		POINT3D new_p, new_nor;
		new_p = (face_points[0] + face_points[1] + face_points[2] + face_points[3]) / 4;

		tmp_add_vh = subvMesh.add_vertex(T_PolygonMesh::Point(new_p[0], new_p[1], new_p[2]));
		tmp_nor[0] = new_nor[0]; tmp_nor[1] = new_nor[1]; tmp_nor[2] = new_nor[2];
		subvMesh.set_normal(tmp_add_vh, tmp_nor);

		new_face_vertices[(*f_it).idx()] = tmp_add_vh;
	}

	T_PolygonMesh::HalfedgeHandle tmp_heh;
	T_PolygonMesh::VertexHandle tmp_vh;
	T_PolygonMesh::FaceHandle tmp_fh;
	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		tmp_heh = mesh->halfedge_handle(e_it, 0);
		tmp_vh = mesh->from_vertex_handle(tmp_heh);
		tmp_fh = mesh->opposite_face_handle(tmp_heh);

		faces_handle.push_back(old_vertices[tmp_vh.idx()]);
		faces_handle.push_back(new_face_vertices[tmp_fh.idx()]);

		tmp_heh = mesh->halfedge_handle(e_it, 1);
		tmp_vh = mesh->from_vertex_handle(tmp_heh);
		tmp_fh = mesh->opposite_face_handle(tmp_heh);

		faces_handle.push_back(old_vertices[tmp_vh.idx()]);
		faces_handle.push_back(new_face_vertices[tmp_fh.idx()]);

		subvMesh.add_face(faces_handle);
		faces_handle.clear();
	}

}

void subdiv::sub_for_quad_mesh_new(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh) {
	subvMesh.clear();

	POINT3D from, to, from_nor, to_nor, newP, newN;
	int id_from, id_to;
	vector<T_PolygonMesh::VertexHandle> origin_vhandle, add_vhandle, facePoint_vhandle;
	Tvhandle.clear();
	int i = 0;
	vector<int> tmp_id;  
	set<id_pair> id_set;

	vector<T_PolygonMesh::VertexHandle> face_vhandles;
	T_PolygonMesh::VertexHandle tmp_vhandle;
	T_PolygonMesh::Normal tmp_nor;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	//保存原先面的顶点
	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmpP;
		tmpP = mesh->point(*v_it);
		tmp_vhandle = subvMesh.add_vertex(T_PolygonMesh::Point(tmpP[0], tmpP[1], tmpP[2]));
		tmp_nor = mesh->normal(*(v_it));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);
		origin_vhandle.push_back(tmp_vhandle);
	}

	//保存每条边新增顶点
	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++)
	{
		from = mesh->point(mesh->from_vertex_handle(mesh->halfedge_handle(e_it, 0)));
		from_nor = mesh->normal(mesh->from_vertex_handle(mesh->halfedge_handle(e_it, 0)));

		to = mesh->point(mesh->to_vertex_handle(mesh->halfedge_handle(e_it, 0)));
		to_nor = mesh->normal(mesh->to_vertex_handle(mesh->halfedge_handle(e_it, 0)));

		create_new_point_and_nor(from, to, from_nor, to_nor, newP, newN);

		tmp_vhandle = subvMesh.add_vertex(T_PolygonMesh::Point(newP[0], newP[1], newP[2]));
		tmp_nor[0] = newN[0];
		tmp_nor[1] = newN[1];
		tmp_nor[2] = newN[2];
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		add_vhandle.push_back(tmp_vhandle);
	}

	//保存每个面顶点
	for (T_PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		Tvhandle.clear();
		for (T_PolygonMesh::FaceEdgeIter fe_it = mesh->fe_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			Tvhandle.push_back(add_vhandle[(*fe_it).idx()]);
		}

		POINT3D p1, p2, p3, p4;
		POINT3D n1, n2, n3, n4;

		p1 = subvMesh.point(Tvhandle[0]);
		p2 = subvMesh.point(Tvhandle[1]); 
		p3 = subvMesh.point(Tvhandle[2]);
		p4 = subvMesh.point(Tvhandle[3]);
		n1 = subvMesh.normal(Tvhandle[0]); n2 = subvMesh.normal(Tvhandle[1]); n3 = subvMesh.normal(Tvhandle[2]); n4 = subvMesh.normal(Tvhandle[3]);
		POINT3D pp1, pp2, pp3, pp4;
		POINT3D nn1, nn2, nn3, nn4;
		create_new_point_and_nor(p1, p2, n1, n2, pp1, nn1);
		create_new_point_and_nor(p2, p3, n2, n3, pp2, nn2);
		create_new_point_and_nor(p3, p4, n3, n4, pp3, nn3);
		create_new_point_and_nor(p4, p1, n4, n1, pp4, nn4);
		
		POINT3D M0, M1, Mn0, Mn1;
		create_new_point_and_nor(pp1, pp3, nn1, nn3, M0, Mn0);
		create_new_point_and_nor(pp2, pp4, nn2, nn4, M1, Mn1);

		POINT3D facePoint = (M0 + M1) / 2;
		
		POINT3D newT1, newT2,newN;
		create_new_T(pp1, pp3, nn1, nn3, newT1);
		create_new_T(pp2, pp4, nn2, nn4, newT2);
		newN = (newT1%newT2).normalize();


		tmp_vhandle = subvMesh.add_vertex(T_PolygonMesh::Point(facePoint[0], facePoint[1], facePoint[2]));
		tmp_nor[0] = newN[0];
		tmp_nor[1] = newN[1];
		tmp_nor[2] = newN[2];
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		facePoint_vhandle.push_back(tmp_vhandle);
	}

	for (T_PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		Tvhandle.clear();
		for (T_PolygonMesh::FaceEdgeIter fe_it = mesh->fe_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			Tvhandle.push_back(add_vhandle[(*fe_it).idx()]);
		}

		for (T_PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			Tvhandle.push_back(origin_vhandle[(*fv_it).idx()]);
		}

		Tvhandle.push_back(facePoint_vhandle[(*f_it).idx()]);


		face_vhandles.clear();
		face_vhandles.push_back(Tvhandle[4]);
		face_vhandles.push_back(Tvhandle[1]);
		face_vhandles.push_back(Tvhandle[8]);
		face_vhandles.push_back(Tvhandle[0]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(Tvhandle[1]);
		face_vhandles.push_back(Tvhandle[5]);
		face_vhandles.push_back(Tvhandle[2]);
		face_vhandles.push_back(Tvhandle[8]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(Tvhandle[8]);
		face_vhandles.push_back(Tvhandle[2]);
		face_vhandles.push_back(Tvhandle[6]);
		face_vhandles.push_back(Tvhandle[3]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(Tvhandle[0]);
		face_vhandles.push_back(Tvhandle[8]);
		face_vhandles.push_back(Tvhandle[3]);
		face_vhandles.push_back(Tvhandle[7]);
		subvMesh.add_face(face_vhandles);
	}
}

void subdiv::Isqrt3(PolygonMesh* mesh, PolygonMesh &subvMesh) {
	subvMesh.clear();

	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<PolygonMesh::VertexHandle> face_vhandles;
	PolygonMesh::VertexHandle tmp_vhandle;
	PolygonMesh::Normal tmp_nor;
	vector<POINT3D>tmp_face_points_and_nors;

	subvMesh.request_vertex_normals();
	subvMesh.update_normals();

	vector<PolygonMesh::VertexHandle> midpoint_vec;
	vector<PolygonMesh::VertexHandle> origin_points_vec;
	vector<vector<int>> origin_face_points_id(faces_num);

	vector<POINT3D> my_point_nors_vec;

	for (PolygonMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_element;
		tmp_element = mesh->point(v_it);
		tmp_face_points_and_nors.push_back(tmp_element);

		//updata v
		//int deg = 0;
		//POINT3D surround_sum(0, 0, 0);
		//for (PolygonMesh::VertexVertexIter vv_it = mesh->vv_iter(v_it); vv_it.is_valid(); vv_it++) {
		//	deg++;
		//	surround_sum += mesh->point(vv_it);
		//}
		//POINT3D newp;
		//double alpha = (4 - 2 * cos(2 * PI / deg)) / 9;
		//newp = (1 - alpha)*tmp_element + alpha / deg * surround_sum;

		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(tmp_element[0], tmp_element[1], tmp_element[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		origin_points_vec.push_back(tmp_vhandle);
	}

	int ind = 0;
	for (PolygonMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		vhandle.clear();
		tmp_face_points_and_nors.clear();
		POINT3D newp(0, 0, 0);
		POINT3D newNor(0, 0, 0);
		POINT3D tmp_element;
		POINT3D point3AreaSum(0, 0, 0);
		POINT3D pointsOnTriangleSum(0, 0, 0);
		POINT3D pointsBoundaryMidSum(0, 0, 0);
		for (auto fh_it = mesh->fh_iter(*f_it); fh_it.is_valid(); ++fh_it) {
			auto opHalfEdge = mesh->opposite_halfedge_handle(fh_it);
			auto boundaryMidVhandle = mesh->opposite_he_opposite_vh(opHalfEdge);
			POINT3D edgeMidPoint;
			edgeMidPoint = mesh->point(boundaryMidVhandle);
			pointsBoundaryMidSum += edgeMidPoint;
		}

		for (PolygonMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			pointsOnTriangleSum += mesh->point(fv_it);

			for (auto fvv_it = mesh->vv_iter(*fv_it); fvv_it.is_valid(); ++fvv_it) {
				point3AreaSum += mesh->point(fvv_it);
			}

			//tmp_element = mesh->point(fv_it);
			//tmp_face_points_and_nors.push_back(tmp_element);

			origin_face_points_id[ind].push_back((*fv_it).idx());
		}

		newp = -2.0 / 81 * (point3AreaSum - 2*pointsBoundaryMidSum - 2*pointsOnTriangleSum) + 32.0 / 81 * pointsOnTriangleSum - 1.0 / 81 * pointsBoundaryMidSum;
		tmp_nor[0] = newNor[0]; tmp_nor[1] = newNor[1]; tmp_nor[2] = newNor[2];
		tmp_vhandle = subvMesh.add_vertex(PolygonMesh::Point(newp[0], newp[1], newp[2]));
		subvMesh.set_normal(tmp_vhandle, tmp_nor);

		midpoint_vec.push_back(tmp_vhandle);
		ind++;
	}

	//add faces
	for (int i = 0; i < faces_num; ++i)
	{
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][1]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		subvMesh.add_face(face_vhandles);
		face_vhandles.clear();
		face_vhandles.push_back(midpoint_vec[i]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][2]]);
		face_vhandles.push_back(origin_points_vec[origin_face_points_id[i][0]]);
		subvMesh.add_face(face_vhandles);
	}

	for (PolygonMesh::EdgeIter e_it = subvMesh.edges_begin(); e_it != subvMesh.edges_end(); ++e_it)
	{

		if (!subvMesh.is_boundary(*e_it))
		{
			int id_from, id_to;
			id_from = (subvMesh.from_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			id_to = (subvMesh.to_vertex_handle(mesh->halfedge_handle(e_it, 0))).idx();
			if (id_from < vertices_num && id_to < vertices_num) //边不包括新生成的面中心点
			{
				subvMesh.flip(*e_it);
			}
		}
	}


#if 0  //翻折狭长三角形且与其相邻三角形近乎垂直的边
	for (PolygonMesh::FaceIter f_it = subvMesh.faces_begin(); f_it != subvMesh.faces_end(); ++f_it)
	{
		vector<PolygonMesh::VertexHandle> tmpPoints;
		vector<PolygonMesh::HalfedgeHandle> tmpEdges;
		for (PolygonMesh::FaceVertexIter fv_it = subvMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			tmpPoints.push_back(*fv_it);
			//POINT3D p;
			//p = subvMesh.point(*fv_it);
			//cout << p<<endl;
		}
		for (PolygonMesh::FaceHalfedgeIter fe_it = subvMesh.fh_iter(*f_it); fe_it.is_valid(); ++fe_it)
		{
			tmpEdges.push_back(*fe_it);
			//POINT3D p;
			//p = subvMesh.point(subvMesh.to_vertex_handle(*fe_it));
			//cout << p << endl;
		}


		POINT3D p1, p2, p3;
		p1 = subvMesh.point(tmpPoints[0]); p2 = subvMesh.point(tmpPoints[1]); p3 = subvMesh.point(tmpPoints[2]);

		int face_flg = is_unregular_face(p1, p2, p3);
		switch (face_flg)
		{
		case 1:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[2]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[2]));
			}
			break;
		case 2:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[0]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[0]));
			}
			break;
		case 3:
			if (subvMesh.is_flip_ok(subvMesh.edge_handle(tmpEdges[1]))) {
				subvMesh.flip(subvMesh.edge_handle(tmpEdges[1]));
			}
			break;
		default:
			break;
		}

	}
#endif
}



bool hasIregularPointsQuad(T_PolygonMesh::FaceHandle fh,T_PolygonMesh *mesh) {
	int deg = 0;
	for (auto vh = mesh->fv_iter(fh); vh.is_valid(); vh++) {
		for (auto it = mesh->vv_iter(vh); it.is_valid(); it++) {
			deg++;
		}
		if (deg != 4) return true;
		deg = 0;
	}
	cout << "has Iregular" << endl;
	return false;
}

void printV(vector<POINT3D>&v) {
	cout << "====: "<<endl;
	for (auto &x : v) {
		cout << x << endl;
	}
	cout << "====: " << endl;
}


void subdiv::Isqrt2(T_PolygonMesh* mesh, T_PolygonMesh &subvMesh) {
	subvMesh.clear();
	int faces_num = mesh->n_faces();
	int vertices_num = mesh->n_vertices();

	vector<T_PolygonMesh::VertexHandle>old_vertices(vertices_num);
	vector<T_PolygonMesh::VertexHandle>new_face_vertices(faces_num);
	T_PolygonMesh::VertexHandle tmp_add_vh;
	T_PolygonMesh::Normal tmp_nor;
	vector<T_PolygonMesh::VertexHandle>faces_handle;

	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
	{
		POINT3D tmp_point, tmp_norr;
		tmp_point = mesh->point(v_it);
		//tmp_nor = mesh->normal(v_it);

		tmp_add_vh = subvMesh.add_vertex(PolygonMesh::Point(tmp_point[0], tmp_point[1], tmp_point[2]));
		subvMesh.set_normal(tmp_add_vh, tmp_nor);
		old_vertices[(*v_it).idx()] = tmp_add_vh;
	}


	//cac_new_face_points
	for (auto f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		vector<POINT3D>face_points;
		vector<POINT3D>face_nors;
		vector<POINT3D>originFacePoints;
		for (auto fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
			POINT3D tmp_p;
			tmp_p = mesh->point(fv_it);
			originFacePoints.push_back(tmp_p);
		}
		//printV(originFacePoints);

		int test_n = 0;
		if (hasIregularPointsQuad(*f_it,mesh)) {
			POINT3D new_p(0, 0, 0),new_nor;
			for (auto fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
				POINT3D quadV;
				quadV = mesh->point(fv_it);
				new_p += 0.5 * quadV;

				int quadVecIndex = 0;
				for (int i = 0; i < 4; ++i) {
					if (quadV == originFacePoints[i]) {
						quadVecIndex = i;
						break;
					}
				}
				POINT3D quadVNext = originFacePoints[(quadVecIndex + 1) % 4]; //获得1邻域内的起始第一个顶点

				for (auto fvv_it = mesh->vv_iter(fv_it); fvv_it.is_valid(); ++fvv_it) {
					POINT3D tmpP;
					tmpP = mesh->point(fvv_it);
					face_points.push_back(tmpP);
				}

				int firstIndex = 0;
				for (int i = 0; i < face_points.size(); ++i) {
					if (face_points[i] == quadVNext) {
						firstIndex = i;
						break;
					}
				}


				int n = face_points.size();
				for (int i = 0; i < n; ++i) {
					float alpha = 1.0 / 2 / n * (1+cos(2 * i*PI / n) + cos(2 * (i - 1)*PI / n) - sin(2 * i*PI / n) - sin(2 * (i - 1)*PI / n));
					new_p += alpha * face_points[(firstIndex - i +  n) % n]; //逆时针
					cout << test_n << ":" << face_points[(firstIndex - i + n) % n] << endl;
					//cout << test_n << " pz:" << new_p[2] << endl;
				}

				test_n++;
				//cout << "face_points" << endl;
				//printV(face_points);
				//cout << "id:" << firstIndex <<endl;
				face_points.clear();
			}

			new_p /= 4.0;
			tmp_add_vh = subvMesh.add_vertex(T_PolygonMesh::Point(new_p[0], new_p[1], new_p[2]));
			tmp_nor[0] = new_nor[0]; tmp_nor[1] = new_nor[1]; tmp_nor[2] = new_nor[2];
			subvMesh.set_normal(tmp_add_vh, tmp_nor);

			new_face_vertices[(*f_it).idx()] = tmp_add_vh;
		}
		else
		{
			for (auto fv_it = mesh->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
			{
				POINT3D tmp_point, tmp_norr;
				tmp_point = mesh->point(fv_it);
				face_points.push_back(tmp_point);
				tmp_norr = mesh->normal(fv_it);
				face_nors.push_back(tmp_norr);
			}

			for (auto fh_it = mesh->fh_iter(*f_it); fh_it.is_valid(); ++fh_it) {
				auto ofh_it = mesh->opposite_halfedge_handle(fh_it);
				auto nextOfh_it = mesh->next_halfedge_handle(ofh_it);
				auto nNextOfh_it = mesh->next_halfedge_handle(nextOfh_it);
				POINT3D t1, t2;
				t1 = mesh->point(mesh->from_vertex_handle(nNextOfh_it));
				t2 = mesh->point(mesh->to_vertex_handle(nNextOfh_it));
				face_points.push_back(t1);
				face_points.push_back(t2);
			}

			POINT3D new_p(0, 0, 0), new_nor;
			for (int i = 0; i < 12; ++i) {
				if (i < 4) {
					new_p += 5.0 / 16 * face_points[i];
				}
				else
				{
					new_p -= 1.0 / 32 * face_points[i];
				}
			}

			tmp_add_vh = subvMesh.add_vertex(T_PolygonMesh::Point(new_p[0], new_p[1], new_p[2]));
			tmp_nor[0] = new_nor[0]; tmp_nor[1] = new_nor[1]; tmp_nor[2] = new_nor[2];
			subvMesh.set_normal(tmp_add_vh, tmp_nor);

			new_face_vertices[(*f_it).idx()] = tmp_add_vh;

		}
	}

	T_PolygonMesh::HalfedgeHandle tmp_heh;
	T_PolygonMesh::VertexHandle tmp_vh;
	T_PolygonMesh::FaceHandle tmp_fh;
	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
	{
		tmp_heh = mesh->halfedge_handle(e_it, 0);
		tmp_vh = mesh->from_vertex_handle(tmp_heh);
		tmp_fh = mesh->opposite_face_handle(tmp_heh);

		faces_handle.push_back(old_vertices[tmp_vh.idx()]);
		faces_handle.push_back(new_face_vertices[tmp_fh.idx()]);

		tmp_heh = mesh->halfedge_handle(e_it, 1);
		tmp_vh = mesh->from_vertex_handle(tmp_heh);
		tmp_fh = mesh->opposite_face_handle(tmp_heh);

		faces_handle.push_back(old_vertices[tmp_vh.idx()]);
		faces_handle.push_back(new_face_vertices[tmp_fh.idx()]);

		subvMesh.add_face(faces_handle);
		faces_handle.clear();
	}

}