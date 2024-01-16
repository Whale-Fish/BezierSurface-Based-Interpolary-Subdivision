#include "NorTest2.h"
#include<string>
#include <QFileDialog.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<stack>
#include"subdiv.h"
#include <sys/timeb.h>
#include<sys/utime.h>
#include "eigenSolution.h"

using namespace std;

#define saved_subdiv_num 4

int SHOW_NORMALS = 0;
int SHOW_SURFACE = 1;
int SHOW_TRANGLES = 0;
int is_Triangle = 1;
int subdiv_cnt = 0;
int closed = 1;

int show_flg = 0;

//全局变量
vector<POINT3D> points, nors;
PolygonMesh Mymesh;
T_PolygonMesh TMymesh;

string cur_file_name;


void OpenFile(string filename, vector<POINT3D>&points, vector<POINT3D>&nors)
{
	int n;
	ifstream bfile(filename);
	string sline;

	getline(bfile, sline);
	istringstream sin(sline);
	sin >> n;
	for (int i = 0; i < n; ++i)
	{
		getline(bfile, sline);
		sin.clear();
		sin.str(sline);
		POINT3D p, nor;
		sin >> p[0] >> p[1] >> p[2] >> nor[0] >> nor[1] >> nor[2];
		points.push_back(p);
		nors.push_back(nor);
	}
	getline(bfile, sline);
	sin.clear();
	sin.str(sline);
	sin >> closed;
	bfile.close();
}

void setMinMax(POINT3D &min, POINT3D &max, const vector<POINT3D>tmp)
{
	min[0] = min[1] = min[2] = 100000000;
	max[0] = max[1] = max[2] = -100000000;
	for (int i = 0; i < tmp.size(); i++)
	{
		if (tmp[i][0] > max[0])max[0] = tmp[i][0];
		if (tmp[i][1] > max[1])max[1] = tmp[i][1];
		if (tmp[i][2] > max[2])max[2] = tmp[i][2];
		if (tmp[i][0] < min[0])min[0] = tmp[i][0];
		if (tmp[i][1] < min[1])min[1] = tmp[i][1];
		if (tmp[i][2] < min[2])min[2] = tmp[i][2];
	}
	
}

void setMinMax(POINT3D &min, POINT3D &max, const POINT3D verpo)
{
	if (verpo[0] > max[0])max[0] = verpo[0];
	if (verpo[1] > max[1])max[1] = verpo[1];
	if (verpo[2] > max[2])max[2] = verpo[2];
	if (verpo[0] < min[0])min[0] = verpo[0];
	if (verpo[1] < min[1])min[1] = verpo[1];
	if (verpo[2] < min[2])min[2] = verpo[2];

}

float Normalize(vector<POINT3D>&tmp)//return scale for the transformation of nors.
{
	POINT3D min, max, dimv, transv;
	float scale;
	setMinMax(min, max, tmp);
	dimv[0] = max[0] - min[0];
	dimv[1] = max[1] - min[1];
	dimv[2] = max[2] - min[2];
	if (dimv[0] >= dimv[1] && dimv[0] > dimv[2])
	{
		scale = 2.0f / dimv[0];
	}
	else if (dimv[1] >= dimv[0] && dimv[1] > dimv[2])
	{
		scale = 2.0f / dimv[1];
	}
	else 
	{
		scale = 2.0f / dimv[2];
	}
	transv[0] = 0.5*(min[0] + max[0]);
	transv[1] = 0.5*(min[1] + max[1]);
	transv[2] = 0.5*(min[2] + max[2]);

	for (int i = 0; i < tmp.size(); i++)
	{
		tmp[i][0] = (tmp[i][0] - transv[0])*scale;
		tmp[i][1] = (tmp[i][1] - transv[1])*scale;
		tmp[i][2] = (tmp[i][2] - transv[2])*scale;
	}

	return scale;
}


POINT3D cac_normal_of_vertex(PolygonMesh::VertexIter v_it, PolygonMesh& mesh)
{
	vector<POINT3D>vf_vertices;
	POINT3D nor(0,0,0);
	for (PolygonMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
	{
		POINT3D tmp;
		tmp = mesh.point((*vv_it));
		vf_vertices.push_back(tmp);
	}
	int sz = vf_vertices.size();
	for (int i = 0; i < sz - 1; ++i)
	{
		nor += vf_vertices[i]%vf_vertices[i + 1];
	}
	nor += vf_vertices[sz - 1]%vf_vertices[0];
	
	return (nor * 1.0 / sz)/(nor * 1.0 / sz).norm();
}

POINT3D cac_normal_of_vertex(PolygonMesh::VertexIter v_it, T_PolygonMesh& mesh)
{
	vector<POINT3D>vf_vertices;
	POINT3D nor(0, 0, 0);
	for (T_PolygonMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
	{
		POINT3D tmp;
		tmp = mesh.point((*vv_it));
		vf_vertices.push_back(tmp);
	}
	int sz = vf_vertices.size();
	for (int i = 0; i < sz - 1; ++i)
	{
		nor += vf_vertices[i] % vf_vertices[i + 1];
	}
	nor += vf_vertices[sz - 1] % vf_vertices[0];

	return (nor * 1.0 / sz) / (nor * 1.0 / sz).norm();
}

void NorTest2::cac_normal_and_write_quad_mesh()
{
	T_PolygonMesh::Normal tnor;
	for (auto v_it = TMymesh.vertices_begin(); v_it != TMymesh.vertices_end(); ++v_it)
	{
		POINT3D cac_nor = -cac_normal_of_vertex(v_it, TMymesh);
		tnor[0] = cac_nor[0]; tnor[1] = cac_nor[1]; tnor[2] = cac_nor[2];
		TMymesh.set_normal(v_it, tnor);
	}

	char mesh_name[80];
	const char* s = cur_file_name.data();
	sprintf(mesh_name, "outputFile//%s_with_nors.obj", s);
	if (!OpenMesh::IO::write_mesh(TMymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
	{
		printf("write mesh file is error!\n");
	}
}

void NorTest2::cac_normal_and_write_tri_mesh()
{
	T_PolygonMesh::Normal tnor;
	for (auto v_it = Mymesh.vertices_begin(); v_it != Mymesh.vertices_end(); ++v_it)
	{
		POINT3D cac_nor = -cac_normal_of_vertex(v_it, Mymesh);
		tnor[0] = cac_nor[0]; tnor[1] = cac_nor[1]; tnor[2] = cac_nor[2];
		Mymesh.set_normal(v_it, tnor);
	}

	char mesh_name[80];
	const char* s = cur_file_name.data();
	sprintf(mesh_name, "outputFile//%s_with_nors.obj", s);
	if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
	{
		printf("write mesh file is error!\n");
	}
}

void normals_comparation_of_limit_surface(PolygonMesh& mesh)
{
	double max_cos = 1.0;
	double mid_cos = 0.0;
	int n = mesh.n_vertices();
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		POINT3D cac_nor = -cac_normal_of_vertex(v_it, mesh);
		POINT3D v_nor;
		v_nor = mesh.normal(*v_it);

		double costheta = cac_nor | v_nor;
		if (isnan(costheta))
		{
			continue;
		}
		if (abs(costheta) < max_cos) { max_cos = abs(costheta); }
		mid_cos += abs(costheta);
	}

	cout << "the mid cos theta of normals:" << mid_cos/n << endl;
	cout << "max_vertices_normal_theta:" << max_cos << endl;
}



void NorTest2::CreateMenus()
{
	file = menuBar()->addMenu(QStringLiteral("Open File"));
	subdivision = menuBar()->addMenu(QStringLiteral("Subdivsion"));
	mode = menuBar()->addMenu(QStringLiteral("Mode"));
	show_or_hide = menuBar()->addMenu(QStringLiteral("Show or hide"));
}


void NorTest2::CreateActions()
{
	open_for_line = new QAction(QStringLiteral("&Line"), this);
	file->addAction(open_for_line);

	open_for_face = new QAction(QStringLiteral("&Tri Surface"), this);
	file->addAction(open_for_face);

	open_for_quad_face = new QAction(QStringLiteral("&Quad Surface"), this);
	file->addAction(open_for_quad_face);

	subdivide_line = new QAction(QStringLiteral("&sub_line"), this);
	subdivision->addAction(subdivide_line);

	subdivide_mesh = new QAction(QStringLiteral("&sub_mesh_nor"), this);
	subdivision->addAction(subdivide_mesh);

	subdivide_mesh_sqrt3 = new QAction(QStringLiteral("&sub_mesh_sqrt3"), this);
	subdivision->addAction(subdivide_mesh_sqrt3);

	subdivide_mesh_sqrt3_PN_method = new QAction(QStringLiteral("&sub_mesh_sqrt3_PN"), this);
	subdivision->addAction(subdivide_mesh_sqrt3_PN_method);

	subdivide_mesh_sqrt2 = new QAction(QStringLiteral("&sub_mesh_sqrt2"), this);
	subdivision->addAction(subdivide_mesh_sqrt2);

	sqrt3 = new QAction(QStringLiteral("&sqrt3"), this);
	subdivision->addAction(sqrt3);

	sqrt2 = new QAction(QStringLiteral("&sqrt2"), this);
	subdivision->addAction(sqrt2);

	Isqrt3 = new QAction(QStringLiteral("&interpolatory sqrt3"), this);
	subdivision->addAction(Isqrt3);

	Isqrt2 = new QAction(QStringLiteral("&interpolatory sqrt2"), this);
	subdivision->addAction(Isqrt2);

	subdivide_quad_mesh = new QAction(QStringLiteral("&sub_quad_mesh_nor"), this);
	subdivision->addAction(subdivide_quad_mesh);

	change_mode = new QAction(QStringLiteral("&mode_change"), this);
	mode->addAction(change_mode);

	//surface_mode = new QAction(QStringLiteral("&surface_mode"), this);
	//mode->addAction(surface_mode);

	nor_show = new QAction(QStringLiteral("nors_show"), this);
	show_or_hide->addAction(nor_show);

	nor_hide = new QAction(QStringLiteral("nors_hide"), this);
	show_or_hide->addAction(nor_hide);

	trangles_show = new QAction(QStringLiteral("trangles_show"), this);
	show_or_hide->addAction(trangles_show);

	trangles_hide = new QAction(QStringLiteral("trangles_hide"), this);
	show_or_hide->addAction(trangles_hide);

	cac_quad_normals = new QAction(QStringLiteral("cac_quad_normals"), this);
	mode->addAction(cac_quad_normals);

	cac_tri_normals = new QAction(QStringLiteral("cac_tri_normals"), this);
	mode->addAction(cac_tri_normals);
}

void NorTest2::file_open_txt()
{
	SHOW_SURFACE = 0;

	open_for_face->setEnabled(false);
	open_for_quad_face->setEnabled(false);
	subdivide_mesh->setEnabled(false);
	subdivide_mesh_sqrt3->setEnabled(false);
	subdivide_mesh_sqrt3_PN_method->setEnabled(false);
	subdivide_mesh_sqrt2->setEnabled(false);


	QString path = QFileDialog::getOpenFileName(this, QStringLiteral("open"), "./square.txt", QStringLiteral("file type (*.txt)"));
	std::string path1 = path.toStdString();	

	if (path1.size() == 0)
	{
		printf("Don't choose file!\n");
		return;
	}

	if (!points.empty() || !nors.empty())
	{
		vector<POINT3D> tmp1, tmp2;
		points.swap(tmp1);
		nors.swap(tmp2);
	}

	OpenFile(path1, points, nors);

	show_flg = 1;
}

void NorTest2::file_open_obj()
{
	SHOW_SURFACE = 1;
	is_Triangle = 1;
	subdiv_cnt = 0;

	open_for_line->setEnabled(false);
	subdivide_line->setEnabled(false);
	subdivide_mesh_sqrt2->setEnabled(false);
	subdivide_mesh_sqrt3->setEnabled(true);
	subdivide_mesh_sqrt3_PN_method->setEnabled(true);
	subdivide_mesh->setEnabled(true);

	OpenMesh::IO::Options opt = OpenMesh::IO::Options::VertexNormal;

	QString path = QFileDialog::getOpenFileName(this, QStringLiteral("open"), "./tri_mesh/cube_nor.obj", QStringLiteral("file type (*.obj);;file type (*.off)"));
	std::string path1 = path.toStdString();
	if (path1.size() == 0)
	{
		printf("Don't choose file!\n");
		return;
	}

	//get cur_file_name
	{
		cur_file_name = "";
		stack<char>stk;
		string tmp_s;
		int s_size = path1.size();
		char c = path1[s_size - 5];
		while (c != '/')
		{
			stk.push(c);
			s_size--;//嫌加变量烦
			c = path1[s_size - 5];
		}
		
		while (!stk.empty())
		{
			char c1 = stk.top();
			cur_file_name += c1;
			stk.pop();
		}

		//cur_file_name = tmp_s.data();
	}



	(&Mymesh)->request_vertex_normals();
	if (!OpenMesh::IO::read_mesh(Mymesh, path1, opt))
	{
		printf("mesh file open error!");
	}

	show_flg = 1;

	//if (true)
	//{
	//	char mesh_name[80];
	//	const char* s = cur_file_name.data();
	//	sprintf(mesh_name, "outputFile//tri_%s.obj", s);
	//	if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
	//	{
	//		printf("write mesh file is error!\n");
	//	}
	//}

	

}

void test1(vector<POINT3D>&p) {
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

	cout << H << endl;

	double p1x = 0, p2x = 1, n1x = 0, n2x = 0;
	double p1y = 0, p2y = 1, n1y = 1, n2y = 1;
	double p1z = 0, p2z = 0, n1z = 0, n2z = 0;

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
	Eigen::MatrixXd x, lambda; //x1,x2,x3,x4,...z3,z4;
	s.qp_lagrange(H, c, A, b, x, lambda, 12, 4);


	p.push_back(POINT3D(p1x, p1y, p1z));
	p.push_back(POINT3D(x(0, 0), x(0, 4), x(0, 8)));
	p.push_back(POINT3D(x(0, 1), x(0, 5), x(0, 9)));
	p.push_back(POINT3D(x(0, 2), x(0, 6), x(0, 10)));
	p.push_back(POINT3D(x(0, 3), x(0, 7), x(0, 11)));
	p.push_back(POINT3D(p2x, p2y, p2z));
}

void test2(vector<POINT3D>&p) {
	Eigen::MatrixXd H = Eigen::MatrixXd::Identity(7, 7);
	H = H * 2;

	double p1x = 0, p2x = 1, n1x = 0, n2x = 0;
	double p1y = 0, p2y = 1, n1y = 1, n2y = 1;
	double p1z = 0, p2z = 0, n1z = 0, n2z = 0;

	double tx = 0.5*(n1x + n2x);
	double ty = 0.5*(n1y + n2y);
	double tz = 0.5*(n1z + n2z);

	POINT3D n1(n1x, n1y, n1z);
	POINT3D n2(n2x, n2y, n2z);

	Eigen::MatrixXd A(5, 7);
	A << (p[1][0] - p[0][0]), (p[1][1] - p[0][1]), (p[1][2] - p[0][2]), (p[3] - p[2]) | n1, 0, 0, 0,
		6 * (p[2][0] - p[1][0]), 6 * (p[2][1] - p[1][1]), 6 * (p[2][2] - p[1][2]), (p[4] - p[3]) | n1, 6 * (p[3] - p[2]) | n1, (p[1] - p[0]) | n2, 0,
		36 * (p[3][0] - p[2][0]), 36 * (p[3][1] - p[2][1]), 36 * (p[3][2] - p[2][2]), (p[5] - p[4]) | n1, 16 * (p[4] - p[3]) | n1, 16 * (p[2] - p[1]) | n2, (p[1] - p[0]) | n2,
		6 * (p[4][0] - p[3][0]), 6 * (p[4][1] - p[3][1]), 6 * (p[4][2] - p[3][2]), 0, (p[5] - p[4]) | n1, 6 * (p[3] - p[2]) | n2, (p[2] - p[1]) | n2,
		(p[5][0] - p[4][0]), (p[5][1] - p[4][1]), (p[5][2] - p[4][2]), 0, 0, 0, (p[3] - p[2]) | n2;


	Eigen::MatrixXd c(7, 1);
	Eigen::MatrixXd b(5, 1);

	c << -2 * tx, -2 * ty, -2 * tz, 0, 0, 0, 0;
	b << 0, 0, 0, 0, 0;

	solution s;
	Eigen::MatrixXd x, lambda; //x,y,z,s0,s1,s3,s4;
	s.qp_lagrange(H, c, A, b, x, lambda, 7, 5);

	cout << "newn:" << x;


}

void NorTest2::file_open_quad_obj()
{
	is_Triangle = 0;
	SHOW_SURFACE = 1;
	subdiv_cnt = 0;

	subdivide_line->setEnabled(false);
	open_for_line->setEnabled(false);
	subdivide_mesh_sqrt3->setEnabled(false);
	subdivide_mesh_sqrt3_PN_method->setEnabled(false);
	subdivide_mesh->setEnabled(false);
	subdivide_mesh_sqrt2->setEnabled(true);

	OpenMesh::IO::Options opt = OpenMesh::IO::Options::VertexNormal;

	QString path = QFileDialog::getOpenFileName(this, QStringLiteral("open"), "./quad_mesh/cube_nor.obj", QStringLiteral("file type (*.obj);;file type (*.off)"));
	std::string path1 = path.toStdString();
	if (path1.size() == 0)
	{
		printf("Don't choose file!\n");
		return;
	}

	//get cur_file_name
	{
		cur_file_name = "";
		stack<char>stk;
		string tmp_s;
		int s_size = path1.size();
		char c = path1[s_size - 5];
		while (c != '/')
		{
			stk.push(c);
			s_size--;//嫌加变量烦
			c = path1[s_size - 5];
		}

		while (!stk.empty())
		{
			char c1 = stk.top();
			cur_file_name += c1;
			stk.pop();
		}

		//cur_file_name = tmp_s.data();
	}



	(&TMymesh)->request_vertex_normals();
	if (!OpenMesh::IO::read_mesh(TMymesh, path1, opt))
	{
		printf("mesh file open error!");
	}

	show_flg = 1;


	//test
	vector<POINT3D>p;
	//test1(p);
	//test2(p);

}

void NorTest2::sub_for_1_time()
{
	subdiv_cnt = 1;

	subdiv sub;

	vector<POINT3D> newpoints(points);
	vector<POINT3D> newnors(nors);
	for (int i = 0; i < subdiv_cnt; i++)
	{
		if (i == subdiv_cnt - 1)
		{
			sub.subdivsion(newpoints, newnors, points, nors);
		}
		else
		{
			sub.subdivsion(newpoints, newnors, points, nors);
			newpoints.swap(points);
			newnors.swap(nors);
		}
	}
	op->update();
}

void NorTest2::sub_for_1_time_mesh()
{
	PolygonMesh newMesh(Mymesh);
    subdiv sub;
	long t1, t2;

	t1 = GetTickCount();


	sub.sub_for_mesh_new(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//nor_outputfile_%s_with_subdivtime_%d.obj",s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}
}

void NorTest2::sub_for_1_time_mesh_sqrt3()
{
	PolygonMesh newMesh(Mymesh);
	subdiv sub;
	//sub.sub_for_mesh(&newMesh,Mymesh);

	//OpenMesh::IO::write_mesh(Mymesh, "123.obj", OpenMesh::IO::Options::VertexNormal);

	long t1, t2;

	t1 = GetTickCount();
	sub.sub_for_mesh_by_sqrt3(&newMesh, Mymesh);
	//sub.sub_for_mesh_by_sqrt3_PN_method(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//sqrt3_outputfile_%s_with_subdivtime_%d.obj", s,subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}

	normals_comparation_of_limit_surface(Mymesh);
}

void NorTest2::sub_for_1_time_mesh_sqrt3_PN_method()
{
	PolygonMesh newMesh(Mymesh);
	subdiv sub;
	//sub.sub_for_mesh(&newMesh,Mymesh);

	//OpenMesh::IO::write_mesh(Mymesh, "123.obj", OpenMesh::IO::Options::VertexNormal);

	long t1, t2;

	t1 = GetTickCount();
	//sub.sub_for_mesh_by_sqrt3(&newMesh, Mymesh);
	sub.sub_for_mesh_by_sqrt3_PN_method(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//sqrt3_PN_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}

	normals_comparation_of_limit_surface(Mymesh);
}

void NorTest2::sub_for_1_time_mesh_sqrt2()
{
	T_PolygonMesh newMesh(TMymesh);
	subdiv sub;
	//sub.sub_for_mesh(&newMesh,Mymesh);

	//OpenMesh::IO::write_mesh(Mymesh, "123.obj", OpenMesh::IO::Options::VertexNormal);





	long t1, t2;

	t1 = GetTickCount();
	sub.sub_for_mesh_by_sqrt2(&newMesh, TMymesh);
	//sub.sub_for_mesh_by_sqrt3_test(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//sqrt2_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(TMymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}

	//normals_comparation_of_limit_surface(TMymesh);
}

void NorTest2::sub_local_quad() {
	T_PolygonMesh newMesh(TMymesh);
	subdiv sub;

	long t1, t2;

	t1 = GetTickCount();
	sub.sub_for_quad_mesh_new(&newMesh, TMymesh);
	//sub.sub_for_mesh_by_sqrt3_test(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//norQuad_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(TMymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}

}

void NorTest2::sub_sqrt3() {
	PolygonMesh newMesh(Mymesh);
	subdiv sub;

	long t1, t2;

	t1 = GetTickCount();
	sub.sqrt3(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("sqrt3 subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//origin_sqrt3_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}
}

void NorTest2::sub_sqrt2() {
	T_PolygonMesh newMesh(TMymesh);
	subdiv sub;

	long t1, t2;

	t1 = GetTickCount();
	sub.sqrt2(&newMesh, TMymesh);
	//sub.sub_for_mesh_by_sqrt3_test(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//sqrt2_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(TMymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}


}

void NorTest2::sub_Isqrt3() {
	PolygonMesh newMesh(Mymesh);
	subdiv sub;

	long t1, t2;

	t1 = GetTickCount();
	sub.Isqrt3(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("Isqrt3 subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//origin_Isqrt3_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(Mymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}
}

void NorTest2::sub_Isqrt2() {
	T_PolygonMesh newMesh(TMymesh);
	subdiv sub;

	long t1, t2;

	t1 = GetTickCount();
	sub.Isqrt2(&newMesh, TMymesh);
	//sub.sub_for_mesh_by_sqrt3_test(&newMesh, Mymesh);
	op->update();

	t2 = GetTickCount();
	printf("this subdivsion cost time: %d\n", t2 - t1);

	subdiv_cnt++;
	if (subdiv_cnt >= saved_subdiv_num)
	{
		char mesh_name[80];
		const char* s = cur_file_name.data();
		sprintf(mesh_name, "outputFile//Isqrt2_outputfile_%s_with_subdivtime_%d.obj", s, subdiv_cnt);
		if (!OpenMesh::IO::write_mesh(TMymesh, mesh_name, OpenMesh::IO::Options::VertexNormal))
		{
			printf("write mesh file is error!\n");
		}
	}


}

void NorTest2::line_mode_change()
{
	SHOW_SURFACE = 0;
	file_open_txt();
	op->update();
}

void NorTest2::mode_change()
{
	if (0 == SHOW_SURFACE && (!points.empty() || !nors.empty()))
	{
		vector<POINT3D> tmp1, tmp2;
		points.swap(tmp1);
		nors.swap(tmp2);
	}
	else if (1 == SHOW_SURFACE)
	{
		Mymesh.clear();
	}

	subdivide_line->setEnabled(true);
	open_for_line->setEnabled(true);
	subdivide_mesh->setEnabled(true);
	open_for_face->setEnabled(true);
	open_for_quad_face->setEnabled(true);
	show_flg = 0;

	op->update();
}

void NorTest2::surface_mode_change()
{
	SHOW_SURFACE = 1;
	file_open_obj();
	op->update();
}

void NorTest2::show_normals()
{
	SHOW_NORMALS = 1;
	op->update();
}

void NorTest2::hide_normals()
{
	SHOW_NORMALS = 0;
	op->update();
}

void NorTest2::show_trangles()
{
	SHOW_TRANGLES = 1;
	op->update();
}

void NorTest2::hide_trangles()
{
	SHOW_TRANGLES = 0;
	op->update();
}


NorTest2::NorTest2(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);
	op = new viewwidget(this);
	
	CreateMenus();
	CreateActions();
	this->setCentralWidget(op);

	connect(open_for_line, SIGNAL(triggered()), this, SLOT(file_open_txt()));
	connect(open_for_face, SIGNAL(triggered()), this, SLOT(file_open_obj()));
	connect(open_for_quad_face, SIGNAL(triggered()), this, SLOT(file_open_quad_obj()));
	connect(subdivide_line, SIGNAL(triggered()), this, SLOT(sub_for_1_time()));
	connect(subdivide_mesh, SIGNAL(triggered()), this, SLOT(sub_for_1_time_mesh()));
	connect(subdivide_mesh_sqrt3, SIGNAL(triggered()), this, SLOT(sub_for_1_time_mesh_sqrt3()));
	connect(subdivide_mesh_sqrt3_PN_method, SIGNAL(triggered()), this, SLOT(sub_for_1_time_mesh_sqrt3_PN_method()));
	connect(subdivide_mesh_sqrt2, SIGNAL(triggered()), this, SLOT(sub_for_1_time_mesh_sqrt2()));
	connect(sqrt3, SIGNAL(triggered()), this, SLOT(sub_sqrt3()));
	connect(sqrt2, SIGNAL(triggered()), this, SLOT(sub_sqrt2()));
	connect(Isqrt3, SIGNAL(triggered()), this, SLOT(sub_Isqrt3()));
	connect(Isqrt2, SIGNAL(triggered()), this, SLOT(sub_Isqrt2()));
	connect(subdivide_quad_mesh, SIGNAL(triggered()), this, SLOT(sub_local_quad()));

	connect(nor_show, SIGNAL(triggered()), this, SLOT(show_normals()));
	connect(nor_hide, SIGNAL(triggered()), this, SLOT(hide_normals()));
	connect(trangles_show, SIGNAL(triggered()), this, SLOT(show_trangles()));
	connect(trangles_hide, SIGNAL(triggered()), this, SLOT(hide_trangles()));

	connect(change_mode, SIGNAL(triggered()), this, SLOT(mode_change()));
	connect(cac_quad_normals, SIGNAL(triggered()), this, SLOT(cac_normal_and_write_quad_mesh()));
	connect(cac_tri_normals, SIGNAL(triggered()), this, SLOT(cac_normal_and_write_tri_mesh()));
	//connect(surface_mode, SIGNAL(triggered()), this, SLOT(surface_mode_change()));

}

NorTest2::~NorTest2()
{



}
