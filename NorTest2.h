#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_NorTest2.h"
#include "viewwidget.h"

class NorTest2 : public QMainWindow
{
    Q_OBJECT

public:
    NorTest2(QWidget *parent = nullptr);
    ~NorTest2();

private:
    Ui::NorTest2Class ui;

	viewwidget *op;
	QMenu *file;
	QMenu *subdivision;
	QMenu *mode;
	QMenu *show_or_hide;

	QAction *open_for_line;
	QAction *open_for_face;  //Tri
	QAction *open_for_quad_face;  //quad
	QAction *subdivide_line;
	QAction *subdivide_mesh;
	QAction *subdivide_quad_mesh;
	QAction *subdivide_mesh_sqrt3;
	QAction *subdivide_mesh_sqrt3_PN_method;
	QAction *subdivide_mesh_sqrt2;

	QAction* sqrt3;
	QAction* sqrt2;
	QAction* Isqrt3;
	QAction* Isqrt2;

	QAction *line_mode;
	QAction *change_mode;
	QAction *surface_mode;
	QAction *nor_show;
	QAction *nor_hide;
	QAction *trangles_show;
	QAction *trangles_hide;

	QAction *cac_quad_normals;
	QAction *cac_tri_normals;

	void CreateMenus();
	void CreateActions();

protected:
	private slots:
		void file_open_txt();
		void file_open_obj();
		void file_open_quad_obj();
		void sub_for_1_time();
		void sub_for_1_time_mesh();
		void sub_local_quad();
		void sub_for_1_time_mesh_sqrt3();
		void sub_for_1_time_mesh_sqrt3_PN_method();
		void sub_for_1_time_mesh_sqrt2();
		void sub_sqrt3();
		void sub_sqrt2();
		void sub_Isqrt3();
		void sub_Isqrt2();
		void line_mode_change();
		void surface_mode_change();
		void show_normals();
		void hide_normals();
		void show_trangles();
		void hide_trangles();
		void mode_change();
		void cac_normal_and_write_tri_mesh();
		void cac_normal_and_write_quad_mesh();
		//void SurfaceSubdivision();

};
