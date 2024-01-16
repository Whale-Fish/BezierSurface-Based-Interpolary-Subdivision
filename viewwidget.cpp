#include "viewwidget.h"
#include "glut.h"
#include<vector>

using namespace std;


int cnt = 0;
bool d_grid = true;
bool d_fill = false;
bool d_subsf = false;
bool d_bsf = true;

extern int SHOW_NORMALS;
extern int SHOW_SURFACE;
extern int SHOW_TRANGLES;
extern int is_Triangle;
extern int closed;
extern vector<POINT3D> points, nors;
extern PolygonMesh Mymesh;
extern T_PolygonMesh TMymesh;
extern int show_flg;




viewwidget::viewwidget(QWidget *parent) :QOpenGLWidget(parent)
{

}

viewwidget::~viewwidget()
{

}

void viewwidget::initializeGL()                         //此处开始对OpenGL进行所以设置
{
	this->initializeOpenGLFunctions();
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);               //黑色背景
	glShadeModel(GL_SMOOTH);                            //启用阴影平滑
	glClearDepth(1.0);                                  //设置深度缓存
	glEnable(GL_DEPTH_TEST);                            //启用深度测试
	glDepthFunc(GL_LEQUAL);                             //所作深度测试的类型
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  //告诉系统对透视进行修正
}
	

void viewwidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(0.0, 1.0, 0.0);

	glLoadIdentity();
	gluLookAt(0, 0, dist_, 0, 0, 0, 0, 1, 0);
	glRotatef(elevation_, 1, 0, 0);
	glRotatef(azimuth_, 0, 1, 0);

	//代码开关3：设置材质与光源
	GLfloat ambient[] = { 0.4, 0.6, 0.2, 1.0 };
	GLfloat position[] = { 1.0, 1.0, 1.0, 1.0 };

	GLfloat ambient1[] = { 0.4, 0.6, 0.2, 1.0 };
	GLfloat position1[] = { -1.0, 1.0, -1.0, 1.0 };

	GLfloat ambient2[] = { 0.4, 0.6, 0.2, 1.0 };
	GLfloat position2[] = { 1.0, -1.0, -1.0, 1.0 };

	glPushMatrix();
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_POSITION, position);

	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient1);
	glLightfv(GL_LIGHT1, GL_POSITION, position1);

	glEnable(GL_LIGHT2);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient2);
	glLightfv(GL_LIGHT1, GL_POSITION, position2);

	GLfloat mat_ambient1[] = { 0.247250, 0.199500, 0.074500, 1.000000, };
	GLfloat mat_diffuse1[] = { 0.751640, 0.606480, 0.226480, 1.000000, };
	GLfloat mat_specular1[] = { 0.628281, 0.555802, 0.366065, 1.000000, };
	GLfloat mat_shininess1[] = { 51.200003, }; //黄铜
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient1);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular1);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess1);

	if (1 == show_flg)
	{
		if (0 == SHOW_SURFACE)
		{
			glBegin(GL_LINE_STRIP);
			glColor3f(1.0, 0.0, 0.0);
			for (int i = 0; i < points.size(); ++i)
			{
				glVertex3d(points[i][0], points[i][1], points[i][2]);
			}

			//是否为闭合曲线
			if (1 == closed) {
				glVertex3d(points[0][0], points[0][1], points[0][2]);
			}
			glEnd();

			if (1 == SHOW_NORMALS)
			{
				glBegin(GL_LINES);
				glColor3f(0.0, 1.0, 0.0);
				for (int i = 0; i < points.size(); ++i)
				{
					glVertex3d(points[i][0], points[i][1], points[i][2]);
					glVertex3d(points[i][0] + nors[i][0], points[i][1] + nors[i][1], points[i][2] + nors[i][2]);
				}
				glEnd();
			}

		}
		else
		{
			glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);

			PolygonMesh::Point from, to;
			int tmp_cnt = 0;
			if (1 == is_Triangle) 
			{
				for (PolygonMesh::EdgeIter e_it = (&Mymesh)->edges_begin(); e_it != (&Mymesh)->edges_end(); e_it++)
				{
					from = (&Mymesh)->point((&Mymesh)->from_vertex_handle((&Mymesh)->halfedge_handle(e_it, 0)));
					to = (&Mymesh)->point((&Mymesh)->to_vertex_handle((&Mymesh)->halfedge_handle(e_it, 0)));
					glVertex3d(from[0], from[1], from[2]);
					glVertex3d(to[0], to[1], to[2]);
				}
			}
			else
			{
				for (T_PolygonMesh::EdgeIter e_it = (&TMymesh)->edges_begin(); e_it != (&TMymesh)->edges_end(); e_it++)
				{
					from = (&TMymesh)->point((&TMymesh)->from_vertex_handle((&TMymesh)->halfedge_handle(e_it, 0)));
					to = (&TMymesh)->point((&TMymesh)->to_vertex_handle((&TMymesh)->halfedge_handle(e_it, 0)));
					glVertex3d(from[0], from[1], from[2]);
					glVertex3d(to[0], to[1], to[2]);
				}
			}
			glEnd();

			//painting normal lines
			if (1 == SHOW_NORMALS)
			{
				glBegin(GL_LINES);
				glColor3f(0.0, 1.0, 0.0);

				PolygonMesh::Point ver, nor;
				if (1 == is_Triangle) 
				{
					for (auto v_it = (&Mymesh)->vertices_begin(); v_it != (&Mymesh)->vertices_end(); v_it++)
					{
						ver = (&Mymesh)->point(*v_it);
						nor = (&Mymesh)->normal(*v_it);
						glVertex3d(ver[0], ver[1], ver[2]);
						glVertex3d(ver[0] + nor[0], ver[1] + nor[1], ver[2] + nor[2]);
					}
				}
				else
				{
					for (auto v_it = (&TMymesh)->vertices_begin(); v_it != (&TMymesh)->vertices_end(); v_it++)
					{
						ver = (&TMymesh)->point(*v_it);
						nor = (&TMymesh)->normal(*v_it);
						glVertex3d(ver[0], ver[1], ver[2]);
						glVertex3d(ver[0] + nor[0], ver[1] + nor[1], ver[2] + nor[2]);
					}
				}

				glEnd();
			}

			if (1 == SHOW_TRANGLES)
			{
				//glBegin(GL_TRIANGLES);
				//PolygonMesh::Point p1, p2, p3;
				//for (auto f_it = (&Mymesh)->faces_begin(); f_it != (&Mymesh)->faces_end(); f_it++)
				//{
				//	auto fv_it = (&Mymesh)->fv_iter(*f_it);
				//	p1 = (&Mymesh)->point(*fv_it);
				//	fv_it++;
				//	p2 = (&Mymesh)->point(*fv_it);
				//	fv_it++;
				//	p3 = (&Mymesh)->point(*fv_it);
				//	glVertex3d(p1[0], p1[1], p1[2]);
				//	glVertex3d(p2[0], p2[1], p2[2]);
				//	glVertex3d(p3[0], p3[1], p3[2]);
				//}
				//glEnd();

				PolygonMesh::Point p;
				if (1 == is_Triangle)
				{
					for (auto f_it = (&Mymesh)->faces_begin(); f_it != (&Mymesh)->faces_end(); f_it++)
					{
						glBegin(GL_POLYGON);
						for (auto fv_it = Mymesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
							p = (&Mymesh)->point(*fv_it);
							glVertex3d(p[0], p[1], p[2]);
						}
						glEnd();
					}
				}
				else
				{
					for (auto f_it = (&TMymesh)->faces_begin(); f_it != (&TMymesh)->faces_end(); f_it++)
					{
						glBegin(GL_POLYGON);
						for (auto fv_it = TMymesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
							p = (&TMymesh)->point(*fv_it);
							glVertex3d(p[0], p[1], p[2]);
						}
						glEnd();
					}
				}

				
			}

		}

	}

	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
}

void viewwidget::resizeGL(int w, int h)
{
	if (h == 0) h = 1; // height == 0 not allowed

	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Calculate The Perspective of the winodw
	gluPerspective(45.0f, (GLfloat)w / (GLfloat)h, 0.1f, 100.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void viewwidget::mouseMoveEvent(QMouseEvent *e)
{
	float dx = e->x() - lastPos.x();
	float dy = e->y() - lastPos.y();

	if (e->buttons() & Qt::LeftButton)
	{
		azimuth_ += 180 * dx / 800;
		elevation_ += 180 * dy / 600;
		update();
	}
	lastPos = e->pos();
}

void viewwidget::mousePressEvent(QMouseEvent *e)
{
	lastPos = e->pos();
}

void viewwidget::wheelEvent(QWheelEvent*wle)
{
	if (wle->delta() > 0)
	{
		dist_ += 0.1;
	}
	else if (wle->delta() < 0)
	{
		dist_ -= 0.1;
	}
	update();

}

