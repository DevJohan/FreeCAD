/***************************************************************************
 *   Copyright (c) 2011 J�rgen Riegel <juergen.riegel@web.de>              *
 *                                                                         *
 *   This file is part of the FreeCAD CAx development system.              *
 *                                                                         *
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Library General Public           *
 *   License as published by the Free Software Foundation; either          *
 *   version 2 of the License, or (at your option) any later version.      *
 *                                                                         *
 *   This library  is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Library General Public License for more details.                  *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this library; see the file COPYING.LIB. If not,    *
 *   write to the Free Software Foundation, Inc., 59 Temple Place,         *
 *   Suite 330, Boston, MA  02111-1307, USA                                *
 *                                                                         *
 ***************************************************************************/


#include "PreCompiled.h"

#ifndef _PreComp_
#endif

#include "ui_TaskSolverStatusHistory.h"
#include "TaskSolverStatusHistory.h"
#include <Gui/Application.h>
#include <Gui/Document.h>
#include <Gui/BitmapFactory.h>
#include <Gui/ViewProvider.h>
#include <Gui/WaitCursor.h>
#include <Gui/Selection.h>

#include <boost/bind.hpp>

#include "ViewProviderSketch.h"

using namespace SketcherGui;
using namespace Gui::TaskView;

TaskSolverStatusHistory::TaskSolverStatusHistory(ViewProviderSketch *sketchView)
    : TaskBox(Gui::BitmapFactory().pixmap("document-new"),
    		tr("Solver status history"),true, 0),
    		sketchView( sketchView )
{
    // we need a separate container widget to add all controls to
    proxy = new QWidget(this);
    ui = new Ui_TaskSolverStatusHistory();
    ui->setupUi(proxy);
    QMetaObject::connectSlotsByName(this);

    this->groupLayout()->addWidget(proxy);

    connectionSolveStatusUpdate =
    		sketchView->signalSolveStatusUpdate.connect(
    				boost::bind(&SketcherGui::TaskSolverStatusHistory::solveStatusUpdate, this,_1,_2,_3));
    connectionSolveStatusUpdate_exp =
    		sketchView->signalSolveStatusUpdate_exp.connect(
    				boost::bind(&SketcherGui::TaskSolverStatusHistory::solveStatusUpdate_exp, this,_1,_2,_3));
}

TaskSolverStatusHistory::~TaskSolverStatusHistory()
{
    connectionSolveStatusUpdate.disconnect();
    connectionSolveStatusUpdate_exp.disconnect();
    delete ui;
}

void TaskSolverStatusHistory::solveStatusUpdate(double solve_time,int solver,bool succeeded ){
	ui->graph_widget->statusUpdate(solve_time, solver, succeeded);
}

void TaskSolverStatusHistory::solveStatusUpdate_exp(double solve_time,int solver,bool succeeded ){
	ui->graph_widget->statusUpdate_exp(solve_time, solver, succeeded);
}

#include "moc_TaskSolverStatusHistory.cpp"
