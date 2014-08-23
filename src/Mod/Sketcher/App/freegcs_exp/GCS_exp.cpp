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
#include "GCS_exp.h"

#include <iostream>
#include <algorithm>
#include <cfloat>
#include <iterator>


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#define BaseExport
#include <Base/Console.h>
#include <Base/TimeInfo.h>



namespace GCS_EXP
{

typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;


struct SubSystemInfo{
	std::set<double*> variables;
	std::set<index_type> constraint_index;
};

///////////////////////////////////////
// Solver
///////////////////////////////////////

// System
System::System():
		constraints_list(0),
//		p2subsystem(),
		dofs(-1),
		subSystems(0),
		hasUnknowns(false), hasDiagnosis(false), isInit(false)
{
}

System::System( std::vector<ConstraintInfo *> clist_):
		constraints_list(0),
//		p2subsystem(),
		dofs(-1),
		subSystems(0),
		hasUnknowns(false), hasDiagnosis(false), isInit(false)
{
	constraints_list.reserve(clist_.size());
    // create own (shallow) copy of constraints
    for (std::vector<ConstraintInfo *>::iterator constr=clist_.begin();
         constr != clist_.end(); ++constr) {
        ConstraintInfo *newconstr = (*constr)->clone();
        if (newconstr)
            addConstraint(newconstr);
    }
}

System::~System()
{
    clear();
}

void System::clear()
{
    hasUnknowns = false;
    hasDiagnosis = false;

    redundant.clear();
    conflictingTags.clear();
    redundantTags.clear();

    clearSubSystems();
    free(constraints_list);
//    p2subsystem.clear();
}

void System::clearByTag(int tagId)
{
    std::vector<ConstraintInfo *> constrvec;
    for (std::vector<ConstraintInfo *>::const_iterator
         constr=constraints_list.begin(); constr != constraints_list.end(); ++constr) {
        if ((*constr)->getTag() == tagId)
            constrvec.push_back( *constr );
    }
    for (std::vector<ConstraintInfo *>::const_iterator
         constr = constrvec.begin(); constr != constrvec.end(); ++constr) {
        removeConstraint(*constr);
    }
    free( constrvec );
}

int System::addConstraint( ConstraintInfo *constr )
{
    isInit = false;
    if ( constr->getTag() >= 0 ) // negatively tagged constraints have no impact
        hasDiagnosis = false;    // on the diagnosis

    constraints_list.push_back(constr);
    return constraints_list.size()-1;
}

void System::removeConstraint(ConstraintInfo *constr)
{
    std::vector<ConstraintInfo *>::iterator it;
    it = std::find( constraints_list.begin(), constraints_list.end(), constr );
    if( it == constraints_list.end() )
        return;

    constraints_list.erase(it);
    if (constr->getTag() >= 0)
        hasDiagnosis = false;
    clearSubSystems();
}

// basic constraints

int System::addConstraintEqual(
		double* param1,
		double* param2,
		int tagId
){
	ConstraintInfo *constr = new ConstraintInfoEqual(
			param1,
			param2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintDifference(
		double* param1,
		double* param2,
		double* difference,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoDifference(
    		param1,
    		param2,
    		difference );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintP2PDistance(
		Point& p1,
		Point& p2,
		double* distance,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoP2PDistance(
    		p1,
    		p2,
    		distance );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintP2PAngle(
		Point &p1,
		Point &p2,
		double* angle,
		double incrAngle,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoP2PAngle(
    		p1,
    		p2,
    		angle,
    		incrAngle);
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintP2PAngle(
		Point &p1,
		Point &p2,
		double* angle,
		int tagId
){
    return addConstraintP2PAngle( p1, p2, angle, 0.);
}

int System::addConstraintP2LDistance(
		Point &p,
		Line &l,
		double* distance,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoP2LDistance(
    		p,
    		l,
    		distance );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintPointOnLine(
		Point &p,
		Line &l,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoPointOnLine(
    		p,
    		l );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintPointOnLine(
		Point& p,
		Point& lp1,
		Point& lp2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoPointOnLine(
    		p,
    		lp1,
    		lp2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintPointOnPerpBisector(
		Point& p,
		Line&  l,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoPointOnPerpBisector(
    		p,
    		l );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintPointOnPerpBisector(
		Point& p,
		Point& lp1,
		Point& lp2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoPointOnPerpBisector(
    		p,
    		lp1,
    		lp2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintParallel(
		Line& l1,
		Line& l2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoParallel(
    		l1,
    		l2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintPerpendicular(
		Line& l1,
		Line& l2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoPerpendicular(
    		l1,
    		l2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintPerpendicular(
		Point& l1p1,
		Point& l1p2,
		Point& l2p1,
		Point& l2p2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoPerpendicular(
    		l1p1,
    		l1p2,
    		l2p1,
    		l2p2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintL2LAngle(
		Line&l1,
		Line&l2,
		double* angle,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoL2LAngle(
    		l1,
    		l2,
    		angle );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintL2LAngle(
		Point& l1p1,
		Point& l1p2,
		Point& l2p1,
		Point& l2p2,
		double* angle,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoL2LAngle(
    		l1p1,
    		l1p2,
    		l2p1,
    		l2p2,
    		angle);
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintMidpointOnLine(
		Line& l1,
		Line& l2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoMidpointOnLine(
    		l1,
    		l2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintMidpointOnLine(
		Point& l1p1,
		Point& l1p2,
		Point& l2p1,
		Point& l2p2,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoMidpointOnLine(
    		l1p1,
    		l1p2,
    		l2p1,
    		l2p2 );
    constr->setTag(tagId);
    return addConstraint(constr);
}

int System::addConstraintTangentCircumf(
		Point& p1,
		Point& p2,
		double* rad1,
		double* rad2,
		bool internal,
		int tagId
){
    ConstraintInfo *constr = new ConstraintInfoTangentCircumf(
    		p1,
    		p2,
    		rad1,
    		rad2,
    		internal );
    constr->setTag(tagId);
    return addConstraint(constr);
}

// derived constraints

int System::addConstraintP2PCoincident(
		Point& p1,
		Point& p2,
		int tagId
){
           addConstraintEqual( p1.x, p2.x, tagId );
    return addConstraintEqual( p1.y, p2.y, tagId );
}

int System::addConstraintHorizontal(
		Line& l,
		int tagId
){
    return addConstraintEqual( l.p1.y, l.p2.y, tagId );
}

int System::addConstraintHorizontal(
		Point& p1,
		Point& p2,
		int tagId
){
    return addConstraintEqual( p1.y, p2.y, tagId );
}

int System::addConstraintVertical(
		Line& l,
		int tagId
){
    return addConstraintEqual(l.p1.x, l.p2.x, tagId);
}

int System::addConstraintVertical(
		Point& p1,
		Point& p2,
		int tagId
){
    return addConstraintEqual(p1.x, p2.x, tagId);
}

int System::addConstraintCoordinateX(
		Point& p,
		double* x,
		int tagId
){
    return addConstraintEqual(p.x, x, tagId);
}

int System::addConstraintCoordinateY(
		Point&p,
		double* y,
		int tagId
){
    return addConstraintEqual( p.y, y, tagId );
}

int System::addConstraintArcRules(
		Arc& a,
		int tagId
){
           addConstraintP2PAngle( a.center, a.start, a.startAngle, tagId );
           addConstraintP2PAngle( a.center, a.end, a.endAngle, tagId );
           addConstraintP2PDistance( a.center, a.start, a.radius, tagId );
    return addConstraintP2PDistance( a.center, a.end, a.radius, tagId );
}

int System::addConstraintPointOnCircle(
		Point&  p,
		Circle& c,
		int tagId
){
    return addConstraintP2PDistance( p, c.center, c.radius, tagId );
}

int System::addConstraintPointOnArc(
		Point& p,
		Arc& a,
		int tagId
){
    return addConstraintP2PDistance( p, a.center, a.radius, tagId );
}

int System::addConstraintPerpendicularLine2Arc(
		Point& p1,
		Point& p2,
		Arc& a,
		int tagId
){
    addConstraintP2PCoincident(p2, a.start, tagId);
    double dx = *(p2.x) - *(p1.x);
    double dy = *(p2.y) - *(p1.y);
    if (dx * cos(*(a.startAngle)) + dy * sin(*(a.startAngle)) > 0)
        return addConstraintP2PAngle( p1, p2, a.startAngle, 0, tagId );
    else
        return addConstraintP2PAngle( p1, p2, a.startAngle, M_PI, tagId );
}

int System::addConstraintPerpendicularArc2Line(
		Arc& a,
		Point& p1,
		Point& p2,
		int tagId
){
    addConstraintP2PCoincident(p1, a.end, tagId);
    double dx = *(p2.x) - *(p1.x);
    double dy = *(p2.y) - *(p1.y);
    if (dx * cos(*(a.endAngle)) + dy * sin(*(a.endAngle)) > 0)
        return addConstraintP2PAngle(p1, p2, a.endAngle, 0, tagId);
    else
        return addConstraintP2PAngle(p1, p2, a.endAngle, M_PI, tagId);
}

int System::addConstraintPerpendicularCircle2Arc(
		Point& center,
		double *radius,
		Arc& a,
		int tagId
){
    addConstraintP2PDistance(a.start, center, radius, tagId);
    double incrAngle = *(a.startAngle) < *(a.endAngle) ? M_PI/2 : -M_PI/2;
    double tangAngle = *(a.startAngle) + incrAngle;
    double dx = *(a.start.x) - *(center.x);
    double dy = *(a.start.y) - *(center.y);
    if (dx * cos(tangAngle) + dy * sin(tangAngle) > 0)
        return addConstraintP2PAngle(center, a.start, a.startAngle, incrAngle, tagId);
    else
        return addConstraintP2PAngle(center, a.start, a.startAngle, -incrAngle, tagId);
}

int System::addConstraintPerpendicularArc2Circle(
		Arc& a,
		Point& center,
		double *radius,
		int tagId
){
    addConstraintP2PDistance(a.end, center, radius, tagId);
    double incrAngle = *(a.startAngle) < *(a.endAngle) ? -M_PI/2 : M_PI/2;
    double tangAngle = *(a.endAngle) + incrAngle;
    double dx = *(a.end.x) - *(center.x);
    double dy = *(a.end.y) - *(center.y);
    if (dx * cos(tangAngle) + dy * sin(tangAngle) > 0)
        return addConstraintP2PAngle(center, a.end, a.endAngle, incrAngle, tagId);
    else
        return addConstraintP2PAngle(center, a.end, a.endAngle, -incrAngle, tagId);
}

int System::addConstraintPerpendicularArc2Arc(
		Arc& a1,
		bool reverse1,
		Arc& a2,
		bool reverse2,
		int tagId
){
    Point& p1 = reverse1 ? a1.start : a1.end;
    Point& p2 = reverse2 ? a2.end : a2.start;
    addConstraintP2PCoincident(p1, p2, tagId);
    return addConstraintPerpendicular(a1.center, p1, a2.center, p2, tagId);
}

int System::addConstraintTangent(
		Line& l,
		Circle& c,
		int tagId
){
    return addConstraintP2LDistance(c.center, l, c.radius, tagId);
}

int System::addConstraintTangent(
		Line& l,
		Arc& a,
		int tagId
){
    return addConstraintP2LDistance(a.center, l, a.radius, tagId);
}

int System::addConstraintTangent(
		Circle& c1,
		Circle& c2,
		int tagId
){
    double dx = *(c2.center.x) - *(c1.center.x);
    double dy = *(c2.center.y) - *(c1.center.y);
    double d = sqrt(dx*dx + dy*dy);
    return addConstraintTangentCircumf(c1.center, c2.center, c1.radius, c2.radius,
                                       (d < *c1.radius || d < *c2.radius), tagId);
}

int System::addConstraintTangent(
		Arc& a1,
		Arc& a2,
		int tagId
){
    double dx = *(a2.center.x) - *(a1.center.x);
    double dy = *(a2.center.y) - *(a1.center.y);
    double d = sqrt(dx*dx + dy*dy);
    return addConstraintTangentCircumf(a1.center, a2.center, a1.radius, a2.radius,
                                       (d < *(a1.radius) || d < *(a2.radius)), tagId);
}

int System::addConstraintTangent(
		Circle& c,
		Arc& a,
		int tagId
){
    double dx = *(a.center.x) - *(c.center.x);
    double dy = *(a.center.y) - *(c.center.y);
    double d = sqrt(dx*dx + dy*dy);
    return addConstraintTangentCircumf(c.center, a.center, c.radius, a.radius,
                                       (d < *(c.radius) || d < *(a.radius)), tagId);
}

int System::addConstraintTangentLine2Arc(
		Point& p1,
		Point& p2,
		Arc& a,
		int tagId
){
    addConstraintP2PCoincident(p2, a.start, tagId);
    double incrAngle = *(a.startAngle) < *(a.endAngle) ? M_PI/2 : -M_PI/2;
    return addConstraintP2PAngle(p1, p2, a.startAngle, incrAngle, tagId);
}

int System::addConstraintTangentArc2Line(
		Arc& a,
		Point& p1,
		Point& p2,
		int tagId
){
    addConstraintP2PCoincident(p1, a.end, tagId);
    double incrAngle = *(a.startAngle) < *(a.endAngle) ? M_PI/2 : -M_PI/2;
    return addConstraintP2PAngle(p1, p2, a.endAngle, incrAngle, tagId);
}

int System::addConstraintTangentCircle2Arc(
		Circle& c,
		Arc& a,
		int tagId
){
    addConstraintPointOnCircle(a.start, c, tagId);
    double dx = *(a.start.x) - *(c.center.x);
    double dy = *(a.start.y) - *(c.center.y);
    if (dx * cos(*(a.startAngle)) + dy * sin(*(a.startAngle)) > 0)
        return addConstraintP2PAngle(c.center, a.start, a.startAngle, 0, tagId);
    else
        return addConstraintP2PAngle(c.center, a.start, a.startAngle, M_PI, tagId);
}

int System::addConstraintTangentArc2Circle(
		Arc& a,
		Circle& c,
		int tagId
){
    addConstraintPointOnCircle(a.end, c, tagId);
    double dx = *(a.end.x) - *(c.center.x);
    double dy = *(a.end.y) - *(c.center.y);
    if (dx * cos( *(a.endAngle) ) + dy * sin( *(a.endAngle) ) > 0)
        return addConstraintP2PAngle(c.center, a.end, a.endAngle, 0, tagId);
    else
        return addConstraintP2PAngle(c.center, a.end, a.endAngle, M_PI, tagId);
}

int System::addConstraintTangentArc2Arc(
		Arc& a1,
		bool reverse1,
		Arc& a2,
		bool reverse2,
		int tagId
){
    Point& p1 = reverse1 ? a1.start : a1.end;
    Point& p2 = reverse2 ? a2.end : a2.start;
    addConstraintP2PCoincident(p1, p2, tagId);

    double* angle1 = reverse1 ? a1.startAngle : a1.endAngle;
    double* angle2 = reverse2 ? a2.endAngle : a2.startAngle;
    if (cos(*(angle1)) * cos(*(angle2)) + sin(*(angle1)) * sin(*(angle2)) > 0)
        return addConstraintEqual( angle1, angle2, tagId);
    else
        return addConstraintP2PAngle(p2, a2.center, angle1, 0, tagId);
}

int System::addConstraintCircleRadius(
		Circle&c,
		double *radius,
		int tagId
){
    return addConstraintEqual( c.radius, radius, tagId );
}

int System::addConstraintArcRadius(
		Arc& a,
		double* radius,
		int tagId
){
    return addConstraintEqual(a.radius, radius, tagId);
}

int System::addConstraintEqualLength(
		Line& l1,
		Line& l2,
		double *length,
		int tagId
){
           addConstraintP2PDistance(l1.p1, l1.p2, length, tagId);
    return addConstraintP2PDistance(l2.p1, l2.p2, length, tagId);
}

int System::addConstraintEqualRadius(
		Circle& c1,
		Circle& c2,
		int tagId
){
    return addConstraintEqual(c1.radius, c2.radius, tagId);
}

int System::addConstraintEqualRadius(
		Circle& c1,
		Arc& a2,
		int tagId
){
    return addConstraintEqual(c1.radius, a2.radius, tagId);
}

int System::addConstraintEqualRadius(
		Arc& a1,
		Arc& a2,
		int tagId
){
    return addConstraintEqual(a1.radius, a2.radius, tagId);
}

int System::addConstraintP2PSymmetric(
		Point& p1,
		Point& p2,
		Line& l,
		int tagId
){
    addConstraintPerpendicular(p1, p2, l.p1, l.p2, tagId);
    return addConstraintMidpointOnLine(p1, p2, l.p1, l.p2, tagId);
}

int System::addConstraintP2PSymmetric(
		Point& p1,
		Point& p2,
		Point& p,
		int tagId
){
    addConstraintPointOnPerpBisector(p, p1, p2, tagId);
    return addConstraintPointOnLine(p, p1, p2, tagId);
}

void System::rescaleConstraint( int id, double coeff )
{
    if ( id >= constraints_list.size() || id < 0 )
        return;
    constraints_list[id]->setScaleCoeff( coeff );

//    SubSystem* constraint_subsys = getConstraintSubSystem(id);
//    if ( constraint_subsys )
//    	constraint_subsys->rescaleConstraint(id,coeff);
}

void System::declareUnknowns( std::vector<double *>& params )
{
    dependent_variables.clear();
    dependent_variables.assign( params.begin(), params.end() );
    sort( dependent_variables.begin(), dependent_variables.end() );
    hasUnknowns = true;
}

void System::initSolution()
{
    typedef std::vector<ConstraintInfo*>::iterator          cinfo_iter;
    typedef std::vector<ConstraintInfo*>::const_iterator    const_cinfo_iter;
    typedef std::vector<double*>::const_iterator      		var_iter;
    // - identifies any decoupled subsystems and partitions the original
    //   system into corresponding components
    // - Stores the current parameters in the vector "reference"
    // - Identifies the equality constraints tagged with ids >= 0
    //   and prepares a corresponding system reduction
    // - Organizes the rest of constraints into two subsystems for
    //   tag ids >=0 and < 0 respectively and applies the
    //   system reduction specified in the previous step
    Base::TimeInfo start_time;

    isInit = false;
    if (!hasUnknowns)
        return;

    // make sure dependent_variables is sorted
    if( dependent_variables.size() > 1 ){
    	for(	var_iter last_dvar = dependent_variables.begin(),
    			dvar = ++dependent_variables.begin();
    			dvar != dependent_variables.end(); ++dvar)
    	{
    		assert( *last_dvar < *dvar );
    		last_dvar = dvar;
    	}
    }

    /* partitioning into decoupled components using graph flood fill coloring.
     * The first dependent_variables.size() vertices in graph correspond to
     * the dependent variables with the same index as in dependent_variables.
     */
    Graph g;
    for (int i=0; i < int(dependent_variables.size() + constraints_list.size()); i++)
        boost::add_vertex(g);

    // first dependent_variables.size() vertices in graph g represents the dependent variables
    int cvtid = dependent_variables.size();
    for (std::vector<ConstraintInfo *>::const_iterator constr=constraints_list.begin();
         constr != constraints_list.end(); ++constr, cvtid++) {
        const std::vector<double*>& constraint_variables = (*constr)->params();
        for( std::vector<double*>::const_iterator param = constraint_variables.begin();
             param != constraint_variables.end(); ++param) {
        	std::vector<double*>::iterator it = std::lower_bound(
        			dependent_variables.begin(),
        			dependent_variables.end(),*param);
            if( it != dependent_variables.end() && *it == *param )
                boost::add_edge(cvtid, std::distance( dependent_variables.begin(), it ), g );
        }
    }

    std::vector<index_type> components(boost::num_vertices(g));
    int componentsSize = 0;
    // This performs the flood fill coloring.
    if (!components.empty())
        componentsSize = boost::connected_components(g, &components[0]);

    // Save which subsystem each constraint belong to
    constraints_to_subsystem.assign(
    		( components.begin() + dependent_variables.size() ),
    		components.end());

    // Resolve which dependent variables belongs to which subsystem
    std::vector< std::vector<double*> > plists( componentsSize );
    for (int i=0; i < int(dependent_variables.size()); ++i) {
        int cid = components[i];
        plists[cid].push_back( dependent_variables[i] );
    }

    // Resolve which constraints belongs to which subsystem
    std::vector< std::vector<ConstraintInfo*> > clists( componentsSize );
    int i = int( dependent_variables.size() );
    for ( std::vector<ConstraintInfo *>::const_iterator constr = constraints_list.begin();
         constr != constraints_list.end(); ++constr) {
    	int cid = components[i++];
    	clists[cid].push_back( *constr );
    }

    // create subsystems
    clearSubSystems();
    for (int cid=0; cid < clists.size(); cid++) {
    	if( clists[cid].size() > 0 )
    		subSystems.push_back( new SubSystem( clists[cid], plists[cid] ) );
    }

    isInit = true;

    // storing reference configuration
    setReference();

    // diagnose conflicting or redundant constraints
    if (!hasDiagnosis) {
        Base::TimeInfo diag_start_time;
    	diagnose();
        Base::TimeInfo diag_end_time;
        Base::Console().Message( "GCS_EXP::diagnose took %f\n",  Base::TimeInfo::diffTimeF(diag_start_time,diag_end_time) );
        if (!hasDiagnosis)
            return;
    }

    std::vector<ConstraintInfo *> clistR;
    if ( redundant.size() > 0 ) {
        for (std::vector<ConstraintInfo *>::const_iterator constr = constraints_list.begin();
             constr != constraints_list.end(); ++constr)
            if ( redundant.count(*constr) == 0 )
                clistR.push_back(*constr);
    }
    else
        clistR = constraints_list;

    Base::TimeInfo end_time;
    Base::Console().Message( "GCS_EXP::initSolution took %f\n",  Base::TimeInfo::diffTimeF(start_time,end_time) );
}

void System::setReference()
{
	typedef std::vector<SubSystem*>::iterator iterator;
    for(iterator it = subSystems.begin();
    		it != subSystems.end(); ++it )
    {
    	if(*it)
    		(*it)->setReference();
    }
}

void System::resetToReference()
{
	typedef std::vector<SubSystem*>::iterator iterator;
    for( iterator it = subSystems.begin();
    		it != subSystems.end(); ++it )
    {
    	if(*it)
    		(*it)->resetToReference();
    }
}

int System::solve( std::vector<double *> &params, bool isFine, Algorithm alg)
{
    declareUnknowns( params );
    initSolution();
    return solve( isFine, alg );
}

int System::solve( bool isFine, Algorithm alg )
{
    if (!isInit)
        return Failed;

    bool isReset = false;
    // return success by default in order to permit coincidence constraints to be applied
    // even if no other system has to be solved
    int res = Success;
    for( int cid=0; cid < int(subSystems.size()); cid++ ){
        if( ( subSystems[cid] ) && !isReset ) {
             resetToReference();
             isReset = true;
        }
        if(subSystems[cid] ){
            res = std::max( res, subSystems[cid]->solve( isFine , alg ) );
        }
    }
    return res;
}



void System::applySolution()
{
    for ( int cid=0; cid < int(subSystems.size()); cid++ ){
        if (subSystems[cid])
            subSystems[cid]->applySolution();
    }
}

void System::undoSolution()
{
    resetToReference();
}

int System::diagnose()
{
    // Analyses the constrainess grad of the system and provides feedback
    // The vector "conflictingTags" will hold a group of conflicting constraints

    // Hint 1: Only constraints with tag >= 0 are taken into account
    // Hint 2: Constraints tagged with 0 are treated as high priority
    //         constraints and they are excluded from the returned
    //         list of conflicting constraints. Therefore, this function
    //         will provide no feedback about possible conflicts between
    //         two high priority constraints. For this reason, tagging
    //         constraints with 0 should be used carefully.
    hasDiagnosis = false;
    if (!hasUnknowns) {
        dofs = -1;
        return dofs;
    }
    dofs = 0;

    redundant.clear();
    conflictingTags.clear();
    redundantTags.clear();
    for(std::vector<SubSystem*>::iterator it = subSystems.begin();
    		it != subSystems.end(); ++it )
    {
    	int subsystem_dof = (*it)->diagnose();
    	if( dofs >= 0 && subsystem_dof >= 0 )
    		dofs += subsystem_dof;
    	else
    		dofs = -1;
    }

    hasDiagnosis = true;
    return dofs;
}

void System::clearSubSystems()
{
    isInit = false;
    free(subSystems);
    subSystems.clear();
}



void free(std::vector<double *> &doublevec)
{
    for (std::vector<double *>::iterator it = doublevec.begin();
         it != doublevec.end(); ++it)
        if (*it) delete *it;
    doublevec.clear();
}

void free( std::vector<ConstraintInfo *> &constrvec )
{
    for (std::vector<ConstraintInfo *>::iterator constr=constrvec.begin();
         constr != constrvec.end(); ++constr) {
        if (*constr) {
            delete *constr;
            *constr = 0;
        }
    }
    constrvec.clear();
}

void free(std::vector<SubSystem *> &subsysvec)
{
    for (std::vector<SubSystem *>::iterator it=subsysvec.begin();
         it != subsysvec.end(); ++it)
        if (*it) { delete *it; *it = 0; }
}


} //namespace GCS_EXP
