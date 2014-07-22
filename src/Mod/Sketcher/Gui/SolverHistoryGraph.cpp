/*
 * SolverHistoryGraph.cpp
 *
 *  Created on: 16 jul 2014
 *      Author: johan
 */

#include <QtGui>
#include <QtGui/qimage.h>
#include "SolverHistoryGraph.h"


namespace SketcherGui {

SolverHistoryGraph::SolverHistoryGraph(QWidget* parent):
		QWidget(parent),
		imageData(1024,150),
		graph(imageData.data(),250,imageData.height(),imageData.byte_stride(),QImage::Format_RGB888),
		current_update(imageData.width()/2),
		current_update_exp(current_update)
{
	this->setMinimumHeight(200);
	this->setMaximumWidth(imageData.width()/2);
	imageData.fill(rgb8(250,230,230));

//	connectionSolveStatusUpdate = sketchView->signalConstraintsChanged.connect(
//        boost::bind(&SketcherGui::TaskSketcherConstrains::slotConstraintsChanged, this));
//    connectionSolveStatusUpdate_exp = sketchView->signalConstraintsChanged.connect(
//        boost::bind(&SketcherGui::TaskSketcherConstrains::slotConstraintsChanged, this));
}

SolverHistoryGraph::~SolverHistoryGraph(){}

void SolverHistoryGraph::paintEvent(
		QPaintEvent* event
){
	QPainter painter(this);
	int width = std::min( this->width()-10, imageData.width()/2 );

	int index = current_update-width;
	graph=QImage(imageData.data()+sizeof(rgb8)*(index),width,imageData.height(),imageData.byte_stride(),QImage::Format_RGB888);

	for( int i=graph.height()-1; i >= 0; i-=10 )
		painter.drawLine(QPoint(7,i),QPoint(9,i));
	for( int i=10; i < graph.width(); i+=10 )
		painter.drawLine(QPoint(i,graph.height()),QPoint(i,graph.height()+3));
	painter.drawImage(QPoint(10,0),graph);
}

void SolverHistoryGraph::statusUpdate(		double solveTime, int solver, bool succeeded ){
	int graph_value = 20*log( 1 + 1e2 * solveTime );
	imageData(imageData.height() - std::min(imageData.height()-2,(1+graph_value)),current_update) = rgb8(50,50,192);;
	imageData(0,current_update_exp) = rgb8(solver == 0 ? 255:0, solver == 1 ? 255:0, solver==2? 255:0);
	++current_update;
	if(current_update >= imageData.width()){
		current_update = imageData.width()/2;
		for(int y=0;y<imageData.height();y++)
			memcpy(imageData.row_data(y),imageData.row_data(y)+( imageData.width() - current_update ), sizeof(rgb8)*( current_update) );
		for( int x = current_update;x < imageData.width(); x++ )
			imageData(0,x) = rgb8(250,230,230);
		for( int y = 1; y < imageData.height(); y++ )
			memcpy( imageData.row_data(y)+current_update, imageData.row_data(0)+current_update, sizeof(rgb8)*( imageData.width() - current_update ) );
	}
	this->repaint();
}

void SolverHistoryGraph::statusUpdate_exp(	double solveTime, int solver, bool succeeded ){
	int graph_value = 20*log( 1 + 1e2 * solveTime );
	imageData(imageData.height() - std::min(imageData.height()-2,(1+graph_value)), current_update_exp) = rgb8(192,192,50);
	imageData(1,current_update_exp) = rgb8(solver == 0 ? 255:0, solver == 1 ? 255:0, solver==2? 255:0);
	++current_update_exp;
	if(current_update_exp >= imageData.width())current_update_exp = imageData.width()/2;
}


}

#include "moc_SolverHistoryGraph.cpp"
