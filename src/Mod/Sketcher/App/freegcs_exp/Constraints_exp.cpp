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

#include <cmath>
#include "Constraints_exp.h"

namespace GCS_EXP
{

///////////////////////////////////////
// Constraints
///////////////////////////////////////

Constraint::Constraint(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count
): 	variables(paramenters),
	variable_indices(indices),
	dependent_variable_count( dependent_var_count ),
	scale(1.),
	tag(0)
{
}

ConstraintType Constraint::getTypeId() const
{
    return None;
}

void Constraint::rescale( double coef )
{
    scale = coef * 1.;
}

double Constraint::maxStep( const std::vector<double>& dir, double lim )
{
    return lim;
}

// Equal
ConstraintEqual::ConstraintEqual(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintEqual::getTypeId() const
{
    return Equal;
}

Constraint* ConstraintEqual::clone() const
{
    return new ConstraintEqual(*this);
}

void ConstraintEqual::rescale(double coef)
{
    scale = coef * 1.;
}

double ConstraintEqual::error()
{
    return scale * ( value<param1>() - value<param2>() );
}

double ConstraintEqual::grad(index_type param)
{
    double deriv=0.;
    if (param == index<param1>()) deriv += 1;
    if (param == index<param2>()) deriv += -1;
    return scale * deriv;
}

void ConstraintEqual::grad( std::vector< grad_component_t >& gradVec )
{
	gradVec.clear();
    dependent_insert<param1>( gradVec,  scale );
    dependent_insert<param2>( gradVec, -scale );
}

// Difference
ConstraintDifference::ConstraintDifference(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters,indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintDifference::getTypeId() const
{
    return Difference;
}

Constraint* ConstraintDifference::clone() const
{
    return new ConstraintDifference(*this);
}

void ConstraintDifference::rescale(double coef)
{
    scale = coef * 1.;
}

double ConstraintDifference::error()
{
    return scale * ( value<param2>() - value<param1>() - value<difference>() );
}

double ConstraintDifference::grad(index_type param)
{
    double deriv=0.;
    if (param == index<param1>()) deriv += -1;
    if (param == index<param2>()) deriv += 1;
    if (param == index<difference>()) deriv += -1;
    return scale * deriv;
}

void ConstraintDifference::grad( std::vector< grad_component_t >& gradVec )
{
    dependent_insert<param1>( gradVec, -scale );
    dependent_insert<param2>( gradVec,  scale );
    dependent_insert<difference>( gradVec, -scale );
}

// P2PDistance
ConstraintP2PDistance::ConstraintP2PDistance(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintP2PDistance::getTypeId() const
{
    return P2PDistance;
}

Constraint* ConstraintP2PDistance::clone() const
{
    return new ConstraintP2PDistance(*this);
}

void ConstraintP2PDistance::rescale(double coef)
{
    scale = coef * 1.;
}

double ConstraintP2PDistance::error()
{
    double dx = (value<p1x>() - value<p2x>());
    double dy = (value<p1y>() - value<p2y>());
    double d = sqrt(dx*dx + dy*dy);
    double dist  = value<distance>();
    return scale * (d - dist);
}

double ConstraintP2PDistance::grad(index_type param)
{
    double deriv=0.;
    if (param == index<p1x>() || param == index<p1y>() ||
        param == index<p2x>() || param == index<p2y>()) {
        double dx = (value<p1x>() - value<p2x>());
        double dy = (value<p1y>() - value<p2y>());
        double d = sqrt(dx*dx + dy*dy);
        if (param == index<p1x>()) deriv += dx/d;
        if (param == index<p1y>()) deriv += dy/d;
        if (param == index<p2x>()) deriv += -dx/d;
        if (param == index<p2y>()) deriv += -dy/d;
    }
    if (param == value<distance>()) deriv += -1.;

    return scale * deriv;
}

void ConstraintP2PDistance::grad( std::vector< grad_component_t >& gradVec )
{
    if( 	is_dependent<p1x>() || is_dependent<p1y>() ||
    		is_dependent<p2x>() || is_dependent<p2y>()) {
        double dx = (value<p1x>() - value<p2x>());
        double dy = (value<p1y>() - value<p2y>());
        double d = sqrt(dx*dx + dy*dy);
        double hdx = scale * dx/d;
        double hdy = scale * dy/d;
        dependent_insert<p1x>( gradVec,  hdx );
        dependent_insert<p1y>( gradVec,  hdy );
        dependent_insert<p2x>( gradVec, -hdx );
        dependent_insert<p2y>( gradVec, -hdy );
    }
    dependent_insert<distance>( gradVec, -scale );
}


double ConstraintP2PDistance::maxStep( const std::vector<double>& dir, double lim )
{
    if( is_dependent<distance>() && dir[index<distance>()] < 0. ) {
            lim = std::min(lim, -(value<distance>()) / dir[ index<distance>() ] );
    }
    // restrict actual value_distance change
    double ddx=0.,ddy=0.;
    if( is_dependent<p1x>() ) ddx += dir[ index<p1x>() ];
    if( is_dependent<p1y>() ) ddy += dir[ index<p1y>() ];
    if( is_dependent<p2x>() ) ddx -= dir[ index<p2x>() ];
    if( is_dependent<p2y>() ) ddy -= dir[ index<p2y>() ];
    double dd = sqrt(ddx*ddx+ddy*ddy);
    double dist  = value<distance>();
    if (dd > dist) {
        double dx = (value<p1x>() - value<p2x>());
        double dy = (value<p1y>() - value<p2y>());
        double d = sqrt(dx*dx + dy*dy);
        if (dd > d)
            lim = std::min(lim, std::max(d,dist)/dd);
    }
    return lim;
}

// P2PAngle
ConstraintP2PAngle::ConstraintP2PAngle(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef,
		double da_
):Constraint(paramenters,indices,dependent_var_count), da(da_){
    rescale( scale_coef );
}

ConstraintType ConstraintP2PAngle::getTypeId() const
{
    return P2PAngle;
}

Constraint* ConstraintP2PAngle::clone() const
{
    return new ConstraintP2PAngle(*this);
}


void ConstraintP2PAngle::rescale(double coef)
{
    scale = coef * 1.;
}

double ConstraintP2PAngle::error()
{
    double dx = (value<p2x>() - value<p1x>());
    double dy = (value<p2y>() - value<p1y>());
    double a = value<angle>() + da;
    double ca = cos(a);
    double sa = sin(a);
    double x = dx*ca + dy*sa;
    double y = -dx*sa + dy*ca;
    return scale * atan2(y,x);
}

double ConstraintP2PAngle::grad(index_type param)
{
    double deriv=0.;
    if (param == index<p1x>() || param == index<p1y>() ||
        param == index<p2x>() || param == index<p2y>()) {
        double dx = (value<p2x>() - value<p1x>());
        double dy = (value<p2y>() - value<p1y>());
        double a = value<angle>() + da;
        double ca = cos(a);
        double sa = sin(a);
        double x = dx*ca + dy*sa;
        double y = -dx*sa + dy*ca;
        double r2 = dx*dx+dy*dy;
        dx = -y/r2;
        dy = x/r2;
        if (param == index<p1x>()) deriv += (-ca*dx + sa*dy);
        if (param == index<p1y>()) deriv += (-sa*dx - ca*dy);
        if (param == index<p2x>()) deriv += ( ca*dx - sa*dy);
        if (param == index<p2y>()) deriv += ( sa*dx + ca*dy);
    }
    if (param == index<angle>()) deriv += -1;

    return scale * deriv;
}

void ConstraintP2PAngle::grad( std::vector< grad_component_t >& gradVec )
{
    if(		is_dependent<p1x>() || is_dependent<p1y>() ||
    		is_dependent<p2x>() || is_dependent<p2y>()) {
        double diff_x = (value<p2x>() - value<p1x>());
        double diff_y = (value<p2y>() - value<p1y>());
        double a = value<angle>() + da;
        double ca = cos(a);
        double sa = sin(a);
        double x =  diff_x * ca + diff_y * sa;
        double y = -diff_x * sa + diff_y * ca;
        double r2 = diff_x*diff_x + diff_y*diff_y;
        double dx = -y/r2;
        double dy =  x/r2;
        dependent_insert<p1x>( gradVec, scale * (-ca*dx + sa*dy) );
        dependent_insert<p1y>( gradVec, scale * (-sa*dx - ca*dy) );
        dependent_insert<p2x>( gradVec, scale * ( ca*dx - sa*dy) );
        dependent_insert<p2y>( gradVec, scale * ( sa*dx + ca*dy) );
    }
    dependent_insert<angle>( gradVec, -scale );

}
double ConstraintP2PAngle::maxStep( const std::vector<double>& dir, double lim)
{
    // step(value_angle()) <= pi/18 = 10°
    if( is_dependent<angle>() ) {
        double step = std::abs( dir[index<angle>()] );
        if( step > M_PI/18. )
            lim = std::min( lim, (M_PI/18.) / step );
    }
    return lim;
}

// P2LDistance
ConstraintP2LDistance::ConstraintP2LDistance(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint(paramenters,indices,dependent_var_count){
    rescale( scale_coef );
}

ConstraintType ConstraintP2LDistance::getTypeId() const
{
    return P2LDistance;
}

Constraint* ConstraintP2LDistance::clone() const
{
    return new ConstraintP2LDistance(*this);
}

void ConstraintP2LDistance::rescale(double coef)
{
    scale = coef;
}

double ConstraintP2LDistance::error()
{
    double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
    double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
    double dist =  value<distance>();
    double dx = x2-x1;
    double dy = y2-y1;
    double d = sqrt(dx*dx+dy*dy);
    double area = std::abs(-x0*dy+y0*dx+x1*y2-x2*y1); // = x1y2 - x2y1 - x0y2 + x2y0 + x0y1 - x1y0 = 2*(triangle area)
    return scale * (area/d - dist);
}

double ConstraintP2LDistance::grad(index_type param)
{
    double deriv=0.;
    // darea/dx0 = (y1-y2)      darea/dy0 = (x2-x1)
    // darea/dx1 = (y2-y0)      darea/dy1 = (x0-x2)
    // darea/dx2 = (y0-y1)      darea/dy2 = (x1-x0)
    if (param == index<px>() || param == index<py>() ||
        param == index<l_p1x>() || param == index<l_p1y>() ||
        param == index<l_p2x>() || param == index<l_p2y>()) {
        double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
        double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
        double dx = x2-x1;
        double dy = y2-y1;
        double d2 = dx*dx+dy*dy;
        double d = sqrt(d2);
        double area = -x0*dy+y0*dx+x1*y2-x2*y1;
        if (param == index<px>()) deriv += (y1-y2) / d;
        if (param == index<py>()) deriv += (x2-x1) / d ;
        if (param == index<l_p1x>()) deriv += ((y2-y0)*d + (dx/d)*area) / d2;
        if (param == index<l_p1y>()) deriv += ((x0-x2)*d + (dy/d)*area) / d2;
        if (param == index<l_p2x>()) deriv += ((y0-y1)*d - (dx/d)*area) / d2;
        if (param == index<l_p2y>()) deriv += ((x1-x0)*d - (dy/d)*area) / d2;
        if (area < 0)
            deriv *= -1;
    }
    if (param == index<distance>()) deriv += -1;

    return scale * deriv;
}

void ConstraintP2LDistance::grad( std::vector< grad_component_t >& gradVec )
{
    // darea/dx0 = (y1-y2)      darea/dy0 = (x2-x1)
    // darea/dx1 = (y2-y0)      darea/dy1 = (x0-x2)
    // darea/dx2 = (y0-y1)      darea/dy2 = (x1-x0)
    if (	is_dependent<px>() 		|| is_dependent<py>() 		||
    		is_dependent<l_p1x>() 	|| is_dependent<l_p1y>() 	||
    		is_dependent<l_p2x>() 	|| is_dependent<l_p2y>()) 	{
        double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
        double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
        double dx = x2-x1;
        double dy = y2-y1;
        double d2 = dx*dx+dy*dy;
        double d = sqrt(d2);
        double area = -x0*dy + y0*dx + x1*y2 - x2*y1;
        double scale_factor = ( area < 0 ) ? -scale : scale;
        double scale_d = scale_factor / d;
        dependent_insert<px>( gradVec, scale_d * (y1-y2) );
        dependent_insert<py>( gradVec, scale_d * (x2-x1) );
        double scale_d2 = scale_factor / d2;
        dependent_insert<l_p1x>( gradVec, scale_d2 * ((y2-y0)*d + (dx/d)*area) );
        dependent_insert<l_p1y>( gradVec, scale_d2 * ((x0-x2)*d + (dy/d)*area) );
        dependent_insert<l_p2x>( gradVec, scale_d2 * ((y0-y1)*d - (dx/d)*area) );
        dependent_insert<l_p2y>( gradVec, scale_d2 * ((x1-x0)*d - (dy/d)*area) );
    }
    dependent_insert<distance>( gradVec, -scale );
}

double ConstraintP2LDistance::maxStep( const std::vector<double>& dir, double lim)
{
	std::map<index_type, double>::iterator it;
    // value_distance() >= 0
    if ( is_dependent<distance>() ) {
        if ( dir[ index<distance>() ] < 0. )
            lim = std::min(lim, -(value<distance>()) / dir[ index<distance>() ]);
    }
    // restrict actual area change
    double darea=0.;
    double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
    double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
    if( is_dependent<px>() ) darea += (y1-y2) * dir[ index<px>() ];
    if( is_dependent<py>() ) darea += (x2-x1) * dir[ index<py>() ];
    if( is_dependent<l_p1x>() ) darea += (y2-y0) * dir[ index<l_p1x>() ];
    if( is_dependent<l_p1y>() ) darea += (x0-x2) * dir[ index<l_p1y>() ];
    if( is_dependent<l_p2x>() ) darea += (y0-y1) * dir[ index<l_p2x>() ];
    if( is_dependent<l_p2y>() ) darea += (x1-x0) * dir[ index<l_p2y>() ];

    darea = std::abs(darea);
    if (darea > 0.) {
        double dx = x2-x1;
        double dy = y2-y1;
        double area = 0.3*( value<distance>() )*sqrt(dx*dx+dy*dy);
        if (darea > area) {
            area = std::max(area, 0.3*std::abs(-x0*dy+y0*dx+x1*y2-x2*y1));
            if (darea > area)
                lim = std::min(lim, area/darea);
        }
    }
    return lim;
}

// PointOnLine
ConstraintPointOnLine::ConstraintPointOnLine(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint(paramenters,indices,dependent_var_count){
    rescale( scale_coef );
}

ConstraintType ConstraintPointOnLine::getTypeId() const
{
    return PointOnLine;
}

Constraint* ConstraintPointOnLine::clone() const
{
    return new ConstraintPointOnLine(*this);
}

void ConstraintPointOnLine::rescale(double coef)
{
    scale = coef;
}

double ConstraintPointOnLine::error()
{
    double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
    double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
    double dx = x2-x1;
    double dy = y2-y1;
    double d = sqrt(dx*dx+dy*dy);
    double area = -x0*dy+y0*dx+x1*y2-x2*y1; // = x1y2 - x2y1 - x0y2 + x2y0 + x0y1 - x1y0 = 2*(triangle area)
    return scale * area/d;
}

double ConstraintPointOnLine::grad(index_type param)
{
    double deriv=0.;
    // darea/dx0 = (y1-y2)      darea/dy0 = (x2-x1)
    // darea/dx1 = (y2-y0)      darea/dy1 = (x0-x2)
    // darea/dx2 = (y0-y1)      darea/dy2 = (x1-x0)
    if (param == index<px>() || param == index<py>() ||
        param == index<l_p1x>() || param == index<l_p1y>() ||
        param == index<l_p2x>() || param == index<l_p2y>()) {
        double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
        double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
        double dx = x2-x1;
        double dy = y2-y1;
        double d2 = dx*dx+dy*dy;
        double d = sqrt(d2);
        double area = -x0*dy+y0*dx+x1*y2-x2*y1;
        if (param == index<px>()) deriv += (y1-y2) / d;
        if (param == index<py>()) deriv += (x2-x1) / d ;
        if (param == index<l_p1x>()) deriv += ((y2-y0)*d + (dx/d)*area) / d2;
        if (param == index<l_p1y>()) deriv += ((x0-x2)*d + (dy/d)*area) / d2;
        if (param == index<l_p2x>()) deriv += ((y0-y1)*d - (dx/d)*area) / d2;
        if (param == index<l_p2y>()) deriv += ((x1-x0)*d - (dy/d)*area) / d2;
    }
    return scale * deriv;
}

void ConstraintPointOnLine::grad( std::vector< grad_component_t >& gradVec )
{
    // darea/dx0 = (y1-y2)      darea/dy0 = (x2-x1)
    // darea/dx1 = (y2-y0)      darea/dy1 = (x0-x2)
    // darea/dx2 = (y0-y1)      darea/dy2 = (x1-x0)
    if(
    		is_dependent<px>() 		|| is_dependent<py>() 		||
    		is_dependent<l_p1x>() 	|| is_dependent<l_p1y>() 	||
    		is_dependent<l_p2x>() 	|| is_dependent<l_p2y>())
    {
        double x0= value<px>(), x1= value<l_p1x>(), x2= value<l_p2x>();
        double y0= value<py>(), y1= value<l_p1y>(), y2= value<l_p2y>();
        double dx = x2-x1;
        double dy = y2-y1;
        double d2 = dx*dx+dy*dy;
        double d = sqrt(d2);
        double area = -x0*dy+y0*dx+x1*y2-x2*y1;
        double scale_d = scale / d;
        dependent_insert<px>( gradVec, scale_d * (y1-y2) );
        dependent_insert<py>( gradVec, scale_d * (x2-x1) );
        double scale_d2 = scale / d2;
        dependent_insert<l_p1x>( gradVec, scale_d2 * ((y2-y0)*d + (dx/d)*area) );
        dependent_insert<l_p1y>( gradVec, scale_d2 * ((x0-x2)*d + (dy/d)*area) );
        dependent_insert<l_p2x>( gradVec, scale_d2 * ((y0-y1)*d - (dx/d)*area) );
        dependent_insert<l_p2y>( gradVec, scale_d2 * ((x1-x0)*d - (dy/d)*area) );
    }
}

// PointOnPerpBisector
ConstraintPointOnPerpBisector::ConstraintPointOnPerpBisector(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintPointOnPerpBisector::getTypeId() const
{
    return PointOnPerpBisector;
}

Constraint* ConstraintPointOnPerpBisector::clone() const
{
    return new ConstraintPointOnPerpBisector(*this);
}

void ConstraintPointOnPerpBisector::rescale(double coef)
{
    scale = coef;
}

double ConstraintPointOnPerpBisector::error()
{
    double dx1 = value<p1x>() - value<p0x>();
    double dy1 = value<p1y>() - value<p0y>();
    double dx2 = value<p2x>() - value<p0x>();
    double dy2 = value<p2y>() - value<p0y>();
    return scale * (sqrt(dx1*dx1+dy1*dy1) - sqrt(dx2*dx2+dy2*dy2));
}

double ConstraintPointOnPerpBisector::grad(index_type param)
{
    double deriv=0.;
    if (param == index<p0x>() || param == index<p0y>() ||
        param == index<p1x>() || param == index<p1y>()) {
        double dx1 = value<p1x>() - value<p0x>();
        double dy1 = value<p1y>() - value<p0y>();
        if (param == index<p0x>()) deriv -= dx1/sqrt(dx1*dx1+dy1*dy1);
        if (param == index<p0y>()) deriv -= dy1/sqrt(dx1*dx1+dy1*dy1);
        if (param == index<p1x>()) deriv += dx1/sqrt(dx1*dx1+dy1*dy1);
        if (param == index<p1y>()) deriv += dy1/sqrt(dx1*dx1+dy1*dy1);
    }
    if (param == index<p0x>() || param == index<p0y>() ||
        param == index<p2x>() || param == index<p2y>()) {
        double dx2 = value<p2x>() - value<p0x>();
        double dy2 = value<p2y>() - value<p0y>();
        if (param == index<p0x>()) deriv += dx2/sqrt(dx2*dx2+dy2*dy2);
        if (param == index<p0y>()) deriv += dy2/sqrt(dx2*dx2+dy2*dy2);
        if (param == index<p2x>()) deriv -= dx2/sqrt(dx2*dx2+dy2*dy2);
        if (param == index<p2y>()) deriv -= dy2/sqrt(dx2*dx2+dy2*dy2);
    }
    return scale * deriv;
}

void ConstraintPointOnPerpBisector::grad( std::vector< grad_component_t >& gradVec )
{
    if(
    		is_dependent<p0x>() || is_dependent<p0y>() ||
    		is_dependent<p1x>() || is_dependent<p1y>())
    {
        double dx1 = value<p1x>() - value<p0x>();
        double dy1 = value<p1y>() - value<p0y>();
        double scale_factor = scale / sqrt(dx1*dx1+dy1*dy1);
        dependent_insert<p0x>( gradVec, -dx1 * scale_factor);
        dependent_insert<p0y>( gradVec, -dy1 * scale_factor);
        dependent_insert<p1x>( gradVec,  dx1 * scale_factor);
        dependent_insert<p1y>( gradVec,  dy1 * scale_factor);
    }
    if (
    		is_dependent<p0x>() || is_dependent<p0y>() ||
    		is_dependent<p2x>() || is_dependent<p2y>())
    {
        double dx2 = value<p2x>() - value<p0x>();
        double dy2 = value<p2y>() - value<p0y>();
        double scale_factor = scale / sqrt(dx2*dx2+dy2*dy2);
        dependent_insert<p0x>( gradVec,  dx2 * scale_factor );
        dependent_insert<p0y>( gradVec,  dy2 * scale_factor );
        dependent_insert<p2x>( gradVec, -dx2 * scale_factor );
        dependent_insert<p2y>( gradVec, -dy2 * scale_factor );
    }
}
// Parallel
ConstraintParallel::ConstraintParallel(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintParallel::getTypeId() const
{
    return Parallel;
}

Constraint* ConstraintParallel::clone() const
{
    return new ConstraintParallel(*this);
}

void ConstraintParallel::rescale(double coef)
{
    double dx1 = (value<l1p1x>() - value<l1p2x>());
    double dy1 = (value<l1p1y>() - value<l1p2y>());
    double dx2 = (value<l2p1x>() - value<l2p2x>());
    double dy2 = (value<l2p1y>() - value<l2p2y>());
    scale = coef / sqrt((dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2));
}

double ConstraintParallel::error()
{
    double dx1 = (value<l1p1x>() - value<l1p2x>());
    double dy1 = (value<l1p1y>() - value<l1p2y>());
    double dx2 = (value<l2p1x>() - value<l2p2x>());
    double dy2 = (value<l2p1y>() - value<l2p2y>());
    return scale * (dx1*dy2 - dy1*dx2);
}

double ConstraintParallel::grad(index_type param)
{
    double deriv=0.;
    if (param == index<l1p1x>()) deriv +=  (value<l2p1y>() - value<l2p2y>()); // = dy2
    if (param == index<l1p2x>()) deriv += -(value<l2p1y>() - value<l2p2y>()); // = -dy2
    if (param == index<l1p1y>()) deriv += -(value<l2p1x>() - value<l2p2x>()); // = -dx2
    if (param == index<l1p2y>()) deriv +=  (value<l2p1x>() - value<l2p2x>()); // = dx2

    if (param == index<l2p1x>()) deriv += -(value<l1p1y>() - value<l1p2y>()); // = -dy1
    if (param == index<l2p2x>()) deriv +=  (value<l1p1y>() - value<l1p2y>()); // = dy1
    if (param == index<l2p1y>()) deriv +=  (value<l1p1x>() - value<l1p2x>()); // = dx1
    if (param == index<l2p2y>()) deriv += -(value<l1p1x>() - value<l1p2x>()); // = -dx1

    return scale * deriv;
}

void ConstraintParallel::grad( std::vector< grad_component_t >& gradVec )
{
		double sdy2 = scale * (value<l2p1y>() - value<l2p2y>());
		dependent_insert<l1p1x>( gradVec,  sdy2 ); // =  dy2
		dependent_insert<l1p2x>( gradVec, -sdy2 ); // = -dy2

    	double sdx2 = scale * (value<l2p1x>() - value<l2p2x>());
    	dependent_insert<l1p1y>( gradVec, -sdx2 ); // = -dx2
    	dependent_insert<l1p2y>( gradVec,  sdx2 ); // =  dx2

    	double sdy1 = scale * (value<l1p1y>() - value<l1p2y>());
    	dependent_insert<l2p1x>( gradVec, -sdy1 ); // = -dy1
    	dependent_insert<l2p2x>( gradVec,  sdy1 ); // =  dy1

    	double sdx1 = scale * (value<l1p1x>() - value<l1p2x>());
    	dependent_insert<l2p1y>( gradVec,  sdx1 ); // = dx1
    	dependent_insert<l2p2y>( gradVec, -sdx1 ); // = -dx1
}

// Perpendicular
ConstraintPerpendicular::ConstraintPerpendicular(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintPerpendicular::getTypeId() const
{
    return Perpendicular;
}

Constraint* ConstraintPerpendicular::clone() const
{
    return new ConstraintPerpendicular(*this);
}

void ConstraintPerpendicular::rescale(double coef)
{
    double dx1 = (value<l1p1x>() - value<l1p2x>());
    double dy1 = (value<l1p1y>() - value<l1p2y>());
    double dx2 = (value<l2p1x>() - value<l2p2x>());
    double dy2 = (value<l2p1y>() - value<l2p2y>());
    scale = coef / sqrt((dx1*dx1+dy1*dy1)*(dx2*dx2+dy2*dy2));
}

double ConstraintPerpendicular::error()
{
    double dx1 = (value<l1p1x>() - value<l1p2x>());
    double dy1 = (value<l1p1y>() - value<l1p2y>());
    double dx2 = (value<l2p1x>() - value<l2p2x>());
    double dy2 = (value<l2p1y>() - value<l2p2y>());
    return scale * (dx1*dx2 + dy1*dy2);
}

double ConstraintPerpendicular::grad(index_type param)
{
    double deriv=0.;
    if (param == index<l1p1x>()) deriv +=  (value<l2p1x>() - value<l2p2x>()); // = dx2
    if (param == index<l1p2x>()) deriv += -(value<l2p1x>() - value<l2p2x>()); // = -dx2
    if (param == index<l1p1y>()) deriv +=  (value<l2p1y>() - value<l2p2y>()); // = dy2
    if (param == index<l1p2y>()) deriv += -(value<l2p1y>() - value<l2p2y>()); // = -dy2

    if (param == index<l2p1x>()) deriv +=  (value<l1p1x>() - value<l1p2x>()); // = dx1
    if (param == index<l2p2x>()) deriv += -(value<l1p1x>() - value<l1p2x>()); // = -dx1
    if (param == index<l2p1y>()) deriv +=  (value<l1p1y>() - value<l1p2y>()); // = dy1
    if (param == index<l2p2y>()) deriv += -(value<l1p1y>() - value<l1p2y>()); // = -dy1

    return scale * deriv;
}

void ConstraintPerpendicular::grad( std::vector< grad_component_t >& gradVec )
{
	double sdx2 = scale * (value<l2p1x>() - value<l2p2x>());
    dependent_insert<l1p1x>( gradVec,  sdx2 ); // = dx2
    dependent_insert<l1p2x>( gradVec, -sdx2 ); // = -dx2

    double sdy2 = scale * (value<l2p1y>() - value<l2p2y>()); // = dy2
	dependent_insert<l1p1y>( gradVec,  sdy2 ); // = dy2
	dependent_insert<l1p2y>( gradVec, -sdy2 ); // = -dy2

    double sdx1 = scale * (value<l1p1x>() - value<l1p2x>()); // = dx1
    dependent_insert<l2p1x>( gradVec,  sdx1 ); // = dx1
    dependent_insert<l2p2x>( gradVec, -sdx1 ); // = -dx1

    double sdy1 = scale * (value<l1p1y>() - value<l1p2y>()); // = dy1
    dependent_insert<l2p1y>( gradVec,  sdy1 ); // = dy1
    dependent_insert<l2p2y>( gradVec, -sdy1 ); // = -dy1
}
// L2LAngle
ConstraintL2LAngle::ConstraintL2LAngle(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintL2LAngle::getTypeId() const
{
    return L2LAngle;
}

Constraint* ConstraintL2LAngle::clone() const
{
    return new ConstraintL2LAngle(*this);
}

void ConstraintL2LAngle::rescale(double coef)
{
    scale = coef * 1.;
}

double ConstraintL2LAngle::error()
{
    double dx1 = (value<l1p2x>() - value<l1p1x>());
    double dy1 = (value<l1p2y>() - value<l1p1y>());
    double dx2 = (value<l2p2x>() - value<l2p1x>());
    double dy2 = (value<l2p2y>() - value<l2p1y>());
    double a = atan2(dy1,dx1) + value<angle>();
    double ca = cos(a);
    double sa = sin(a);
    double x2 = dx2*ca + dy2*sa;
    double y2 = -dx2*sa + dy2*ca;
    return scale * atan2(y2,x2);
}

double ConstraintL2LAngle::grad(index_type param)
{
    double deriv=0.;
    if (param == index<l1p1x>() || param == index<l1p1y>() ||
        param == index<l1p2x>() || param == index<l1p2y>()) {
        double dx1 = (value<l1p2x>() - value<l1p1x>());
        double dy1 = (value<l1p2y>() - value<l1p1y>());
        double r2 = dx1*dx1+dy1*dy1;
        if (param == index<l1p1x>()) deriv += -dy1/r2;
        if (param == index<l1p1y>()) deriv += dx1/r2;
        if (param == index<l1p2x>()) deriv += dy1/r2;
        if (param == index<l1p2y>()) deriv += -dx1/r2;
    }
    if (param == index<l2p1x>() || param == index<l2p1y>() ||
        param == index<l2p2x>() || param == index<l2p2y>()) {
        double dx1 = (value<l1p2x>() - value<l1p1x>());
        double dy1 = (value<l1p2y>() - value<l1p1y>());
        double dx2 = (value<l2p2x>() - value<l2p1x>());
        double dy2 = (value<l2p2y>() - value<l2p1y>());
        double a = atan2(dy1,dx1) + value<angle>();
        double ca = cos(a);
        double sa = sin(a);
        double x2 = dx2*ca + dy2*sa;
        double y2 = -dx2*sa + dy2*ca;
        double r2 = dx2*dx2+dy2*dy2;
        dx2 = -y2/r2;
        dy2 = x2/r2;
        if (param == index<l2p1x>()) deriv += (-ca*dx2 + sa*dy2);
        if (param == index<l2p1y>()) deriv += (-sa*dx2 - ca*dy2);
        if (param == index<l2p2x>()) deriv += ( ca*dx2 - sa*dy2);
        if (param == index<l2p2y>()) deriv += ( sa*dx2 + ca*dy2);
    }
    if (param == index<angle>()) deriv += -1;

    return scale * deriv;
}

void ConstraintL2LAngle::grad( std::vector< grad_component_t >& gradVec )
{
    if(		is_dependent<l1p1x>() 	|| is_dependent<l1p1y>() ||
    		is_dependent<l1p2x>() 	|| is_dependent<l1p2y>())
    {
        double dx1 = (value<l1p2x>() - value<l1p1x>());
        double dy1 = (value<l1p2y>() - value<l1p1y>());
        double r2 = dx1*dx1+dy1*dy1;
        double scale_factor = scale / r2;
        dependent_insert<l1p1x>( gradVec, -scale_factor * dy1 );
        dependent_insert<l1p1y>( gradVec,  scale_factor * dx1 );
        dependent_insert<l1p2x>( gradVec,  scale_factor * dy1 );
        dependent_insert<l1p2y>( gradVec, -scale_factor * dx1 );
    }
    if(		is_dependent<l2p1x>() 	|| is_dependent<l2p1y>() ||
    		is_dependent<l2p2x>() 	|| is_dependent<l2p2y>())
    {
        double diff_x1 = (value<l1p2x>() - value<l1p1x>());
        double diff_y1 = (value<l1p2y>() - value<l1p1y>());
        double diff_x2 = (value<l2p2x>() - value<l2p1x>());
        double diff_y2 = (value<l2p2y>() - value<l2p1y>());
        double a = atan2(diff_y1,diff_x1) + value<angle>();
        double ca = cos(a);
        double sa = sin(a);
        double x2 =  diff_x2 * ca + diff_y2 * sa;
        double y2 = -diff_x2 * sa + diff_y2 * ca;
        double r2 = diff_x2*diff_x2 + diff_y2*diff_y2;
        double dx2 = -y2/r2;
        double dy2 =  x2/r2;
        dependent_insert<l2p1x>( gradVec, scale * (-ca*dx2 + sa*dy2) );
        dependent_insert<l2p1y>( gradVec, scale * (-sa*dx2 - ca*dy2) );
        dependent_insert<l2p2x>( gradVec, scale * ( ca*dx2 - sa*dy2) );
        dependent_insert<l2p2y>( gradVec, scale * ( sa*dx2 + ca*dy2) );
    }
    dependent_insert<angle>( gradVec, -scale );
}

double ConstraintL2LAngle::maxStep( const std::vector<double>& dir, double lim)
{
    // step(value_angle()) <= pi/18 = 10°
    if( is_dependent<angle>() ) {
        double step = std::abs( dir[ index<angle>() ]);
        if (step > M_PI/18.)
            lim = std::min(lim, (M_PI/18.) / step);
    }
    return lim;
}

// MidpointOnLine
ConstraintMidpointOnLine::ConstraintMidpointOnLine(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef
):Constraint( paramenters, indices, dependent_var_count ){
    rescale( scale_coef );
}

ConstraintType ConstraintMidpointOnLine::getTypeId() const
{
    return MidpointOnLine;
}

Constraint* ConstraintMidpointOnLine::clone() const
{
    return new ConstraintMidpointOnLine(*this);
}

void ConstraintMidpointOnLine::rescale(double coef)
{
    scale = coef * 1;
}

double ConstraintMidpointOnLine::error()
{
    double x0=( value<l1p1x>() + value<l1p2x>() )/2;
    double y0=( value<l1p1y>() + value<l1p2y>() )/2;
    double x1= value<l2p1x>(), x2= value<l2p2x>();
    double y1= value<l2p1y>(), y2= value<l2p2y>();
    double dx = x2-x1;
    double dy = y2-y1;
    double d = sqrt(dx*dx+dy*dy);
    double area = -x0*dy+y0*dx+x1*y2-x2*y1; // = x1y2 - x2y1 - x0y2 + x2y0 + x0y1 - x1y0 = 2*(triangle area)
    return scale * area/d;
}

double ConstraintMidpointOnLine::grad(index_type param)
{
    double deriv=0.;
    // darea/dx0 = (y1-y2)      darea/dy0 = (x2-x1)
    // darea/dx1 = (y2-y0)      darea/dy1 = (x0-x2)
    // darea/dx2 = (y0-y1)      darea/dy2 = (x1-x0)
    if (param == index<l1p1x>() || param == index<l1p1y>() ||
        param == index<l1p2x>() || param == index<l1p2y>()||
        param == index<l2p1x>() || param == index<l2p1y>() ||
        param == index<l2p2x>() || param == index<l2p2y>()) {
        double x0=((value<l1p1x>())+(value<l1p2x>()))/2;
        double y0=((value<l1p1y>())+(value<l1p2y>()))/2;
        double x1=value<l2p1x>(), x2=value<l2p2x>();
        double y1=value<l2p1y>(), y2=value<l2p2y>();
        double dx = x2-x1;
        double dy = y2-y1;
        double d2 = dx*dx+dy*dy;
        double d = sqrt(d2);
        double area = -x0*dy+y0*dx+x1*y2-x2*y1;
        if (param == index<l1p1x>()) deriv += (y1-y2) / (2*d);
        if (param == index<l1p1y>()) deriv += (x2-x1) / (2*d);
        if (param == index<l1p2x>()) deriv += (y1-y2) / (2*d);
        if (param == index<l1p2y>()) deriv += (x2-x1) / (2*d);
        if (param == index<l2p1x>()) deriv += ((y2-y0)*d + (dx/d)*area) / d2;
        if (param == index<l2p1y>()) deriv += ((x0-x2)*d + (dy/d)*area) / d2;
        if (param == index<l2p2x>()) deriv += ((y0-y1)*d - (dx/d)*area) / d2;
        if (param == index<l2p2y>()) deriv += ((x1-x0)*d - (dy/d)*area) / d2;
    }
    return scale * deriv;
}

void ConstraintMidpointOnLine::grad( std::vector< grad_component_t >& gradVec )
{
    double deriv=0.;
    // darea/dx0 = (y1-y2)      darea/dy0 = (x2-x1)
    // darea/dx1 = (y2-y0)      darea/dy1 = (x0-x2)
    // darea/dx2 = (y0-y1)      darea/dy2 = (x1-x0)
    if(
    		is_dependent<l1p1x>() || is_dependent<l1p1y>() ||
    		is_dependent<l1p2x>() || is_dependent<l1p2y>()||
    		is_dependent<l2p1x>() || is_dependent<l2p1y>() ||
    		is_dependent<l2p2x>() || is_dependent<l2p2y>())
    {
        double x0=((value<l1p1x>())+(value<l1p2x>()))/2;
        double y0=((value<l1p1y>())+(value<l1p2y>()))/2;
        double x1=value<l2p1x>(), x2=value<l2p2x>();
        double y1=value<l2p1y>(), y2=value<l2p2y>();
        double dx = x2-x1;
        double dy = y2-y1;
        double d2 = dx*dx + dy*dy;
        double d = sqrt(d2);
        double area = -x0*dy + y0*dx + x1*y2 - x2*y1;

        double scale_d = scale / (2*d);
        dependent_insert<l1p1x>( gradVec, scale_d * (y1-y2) );
        dependent_insert<l1p1y>( gradVec, scale_d * (x2-x1) );
        dependent_insert<l1p2x>( gradVec, scale_d * (y1-y2) );
        dependent_insert<l1p2y>( gradVec, scale_d * (x2-x1) );

        double scale_d2 = scale / d2;
        dependent_insert<l2p1x>( gradVec, scale_d2 * ((y2-y0)*d + (dx/d)*area) );
        dependent_insert<l2p1y>( gradVec, scale_d2 * ((x0-x2)*d + (dy/d)*area) );
        dependent_insert<l2p2x>( gradVec, scale_d2 * ((y0-y1)*d - (dx/d)*area) );
        dependent_insert<l2p2y>( gradVec, scale_d2 * ((x1-x0)*d - (dy/d)*area) );
    }
}

// TangentCircumf
ConstraintTangentCircumf::ConstraintTangentCircumf(
		const std::vector<double>& paramenters,
		const std::vector<index_type> indices,
		index_type dependent_var_count,
		double scale_coef,
		bool internal_
):Constraint( paramenters, indices, dependent_var_count ),
		internal(internal_){
//    internal = internal_;
    rescale(scale_coef);
}

ConstraintType ConstraintTangentCircumf::getTypeId() const
{
    return TangentCircumf;
}

Constraint* ConstraintTangentCircumf::clone() const
{
    return new ConstraintTangentCircumf(*this);
}

void ConstraintTangentCircumf::rescale(double coef)
{
    scale = coef * 1;
}

double ConstraintTangentCircumf::error()
{
    double dx = (value<c1x>() - value<c2x>());
    double dy = (value<c1y>() - value<c2y>());
    if (internal)
        return scale * (sqrt(dx*dx + dy*dy) - std::abs( value<r1>() - value<r2>() ));
    else
        return scale * (sqrt(dx*dx + dy*dy) - ( value<r1>() + value<r2>() ));
}

double ConstraintTangentCircumf::grad( index_type param )
{
    double deriv=0.;
    if (param == index<c1x>() || param == index<c1y>() ||
        param == index<c2x>() || param == index<c2y>()||
        param == index<r1>() || param == index<r2>()) {
        double dx = (value<c1x>() - value<c2x>());
        double dy = (value<c1y>() - value<c2y>());
        double d = sqrt(dx*dx + dy*dy);
        if (param == index<c1x>()) deriv += dx/d;
        if (param == index<c1y>()) deriv += dy/d;
        if (param == index<c2x>()) deriv += -dx/d;
        if (param == index<c2y>()) deriv += -dy/d;
        if (internal) {
            if (param == index<r1>()) deriv += (value<r1>() > value<r2>()) ? -1 : 1;
            if (param == index<r2>()) deriv += (value<r1>() > value<r2>()) ? 1 : -1;
        }
        else {
            if (param == index<r1>()) deriv += -1;
            if (param == index<r2>()) deriv += -1;
        }
    }
    return scale * deriv;
}

void ConstraintTangentCircumf::grad( std::vector< grad_component_t >& gradVec )
{
    double deriv=0.;
    if(
    		is_dependent<c1x>() || is_dependent<c1y>() ||
    		is_dependent<c2x>() || is_dependent<c2y>()||
    		is_dependent<r1>() 	|| is_dependent<r2>())
    {
        double dx = (value<c1x>() - value<c2x>());
        double dy = (value<c1y>() - value<c2y>());
        double d = sqrt(dx*dx + dy*dy);
        double scale_factor = scale / d;
        dependent_insert<c1x>( gradVec,  scale_factor * dx);
        dependent_insert<c1y>( gradVec,  scale_factor * dy);
        dependent_insert<c2x>( gradVec, -scale_factor * dx);
        dependent_insert<c2y>( gradVec, -scale_factor * dy);
        if (internal) {
        	dependent_insert<r1>( gradVec, (value<r1>() > value<r2>()) ? -scale :  scale );
        	dependent_insert<r2>( gradVec, (value<r1>() > value<r2>()) ?  scale : -scale );
        }
        else {
        	dependent_insert<r1>( gradVec, -scale);
        	dependent_insert<r2>( gradVec, -scale);
        }
    }
}

} //namespace GCS
