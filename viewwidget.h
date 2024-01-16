#pragma once
#include <QtWidgets/QOpenGLWidget>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QOpenGLBuffer>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QMatrix4x4>
#include <QtGui/QWindow>
#include <Qt3DInput/QMouseEvent>  
#include<GL\GLU.h>
#include "Kernel.h"

class viewwidget :public QOpenGLWidget, protected QOpenGLFunctions
{
public:
	viewwidget(QWidget* parent);
	~viewwidget();

	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();
	void paintGLtest();
	void mouseMoveEvent(QMouseEvent *e);
	void mousePressEvent(QMouseEvent *e);
	void wheelEvent(QWheelEvent *wle);
	QPoint lastPos;

private:
	float elevation_ = 3.0;
	float azimuth_ = 0.0;
	float dist_ = 3.0;
};

