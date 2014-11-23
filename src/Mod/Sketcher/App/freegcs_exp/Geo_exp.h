/***************************************************************************
 *   Copyright (c) Konstantinos Poulios      (logari81@gmail.com) 2011     *
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

#ifndef FREEGCS_GEO_EXP_H
#define FREEGCS_GEO_EXP_H

namespace GCS_EXP
{

    ///////////////////////////////////////
    // Geometries
    ///////////////////////////////////////

    class Point
    {
    public:
        Point(){x = 0; y = 0;}
        double *x;
        double *y;
    };

    class Line
    {
    public:
        Line(){}
        Point p1;
        Point p2;
    };

    class Arc
    {
    public:
        Arc(){startAngle=0;endAngle=0;radius=0;}
        double *startAngle;
        double *endAngle;
        double *radius;
        Point start;
        Point end;
        Point center;
    };

    class Circle
    {
    public:
        Circle(){radius = 0;}
        Point center;
        double *radius;
    };

    class Ellipse
    {
    public:
        Ellipse(){ radmin = 0; focus1X = 0; focus1Y = 0; }
        Point center;
        double *focus1X; //+u
        double *focus1Y;
        double *radmin;
    };

    class ArcOfEllipse
    {
    public:
        ArcOfEllipse(){ startAngle=0; endAngle=0; radmin = 0; focus1X = 0; focus1Y = 0; }
        double *startAngle;
        double *endAngle;
        double *radmin;
        Point start;
        Point end;
        Point center;
        double *focus1X; //+u
        double *focus1Y;
    };


    ///////////////////////////////////////
    // Geometries info
    ///////////////////////////////////////
    typedef int index_type;

    class PointInfo
    {
    public:
    	PointInfo(index_type x_i=0, index_type y_i=0):x_index(x_i),y_index(y_i){}
    	index_type x_index;
    	index_type y_index;
    };

    class LineInfo
    {
    public:
    	LineInfo(PointInfo pa, PointInfo pb):p1(pa),p2(pb){}
    	PointInfo p1;
    	PointInfo p2;
    };

    class ArcInfo
    {
    public:
    	ArcInfo(	index_type startAngle,
    			index_type endAngle,
    			index_type rad,
    			PointInfo start,
    			PointInfo end,
    			PointInfo center
    	):
    		startAngle(startAngle),
    		endAngle(endAngle),
    		rad(rad),
    		start(start),
    		end(end),
    		center(center)
    	{}
    	index_type startAngle;
    	index_type endAngle;
    	index_type rad;
    	PointInfo start;
    	PointInfo end;
    	PointInfo center;
    };

    class CircleInfo
    {
    public:
    	CircleInfo(PointInfo c_in, index_type ra):center(c_in),radius(ra){ }
    	PointInfo center;
    	index_type radius;
    };


    class EllipseInfo
    {
    public:
        EllipseInfo(
        		PointInfo _center,
        		index_type _focus1X,
        		index_type _focus1Y,
        		index_type _radmin
        ):
        	center(_center),
        	focus1X(_focus1X),
        	focus1Y(_focus1Y),
        	radmin(_radmin)
        {}
        PointInfo center;
        index_type focus1X; //+u
        index_type focus1Y;
        index_type radmin;
    };

    class ArcOfEllipseInfo
    {
    public:
    	ArcOfEllipseInfo(
    			index_type _startAngle,
    			index_type _endAngle,
    			index_type _radmin,
    			PointInfo _start,
    			PointInfo _end,
    			PointInfo _center,
    			index_type _focus1X,
    			index_type _focus1Y
    	):
    		startAngle(_startAngle),
    		endAngle(_endAngle),
    		radmin(radmin),
    		start(_start),
    		end(_end),
    		center(_center),
    		focus1X(_focus1X),
    		focus1Y(_focus1Y)
    	{}

        index_type startAngle;
        index_type endAngle;
        index_type radmin;
        PointInfo start;
        PointInfo end;
        PointInfo center;
        index_type focus1X; //+u
        index_type focus1Y;
    };

} //namespace GCS_EXP

#endif // FREEGCS_GEO_EXP_H
