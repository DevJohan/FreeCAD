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


#ifndef GUI_TASKVIEW_TaskSolverStatusHistory_H
#define GUI_TASKVIEW_TaskSolverStatusHistory_H

#include <Gui/TaskView/TaskView.h>
#include <Gui/Selection.h>
#include <boost/signals2.hpp>

class Ui_TaskSolverStatusHistory;

namespace App {
class Property;
}

namespace SketcherGui { 

class ViewProviderSketch;

class TaskSolverStatusHistory : public Gui::TaskView::TaskBox
{
    Q_OBJECT

public:
    TaskSolverStatusHistory(ViewProviderSketch *sketchView);
    ~TaskSolverStatusHistory();

    void solveStatusUpdate(		double solve_time, int solver, bool succeeded );
    void solveStatusUpdate_exp(	double solve_time, int solver, bool succeeded );

private Q_SLOTS:
    
protected:
    ViewProviderSketch *sketchView;
    Connection connectionSolveStatusUpdate;
    Connection connectionSolveStatusUpdate_exp;

private:
    QWidget* proxy;
    Ui_TaskSolverStatusHistory* ui;
};

} //namespace SketcherGui

#endif // GUI_TASKVIEW_TaskSolverStatusHistory_H
