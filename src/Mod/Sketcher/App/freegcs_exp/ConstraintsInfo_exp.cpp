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

#include "ConstraintsInfo_exp.h"

#include <cmath>

#include "Constraints_exp.h"
#include "SubSystem_exp.h"

namespace GCS_EXP
{

///////////////////////////////////////
// ConstraintsInfo
///////////////////////////////////////

ConstraintInfo::ConstraintInfo( )
: variables(), tag(0), scale_coef(1.0)
{
}

ConstraintType ConstraintInfo::getTypeId() const
{
    return None;
}

bool ConstraintInfo::operator< ( const ConstraintInfo& other ) const{
	ConstraintType typeId = getTypeId();
	ConstraintType typeIdOther = other.getTypeId();
	if( typeId < typeIdOther )
		return true;
	else if( typeIdOther < typeId )
		return false;

	assert( variables.size() == other.variables.size() );
	for( int i = 0; i < variables.size(); ++i ){
		if( variables[i] < other.variables[i] )
			return true;
		else if( other.variables[i] < variables[i] )
			return false;
	}
	return false;
}

// EqualInfo
ConstraintInfoEqual::ConstraintInfoEqual(
		double* const p1,
		double* const p2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::param1 ] = p1;
    variables[ ConstraintVariables::param2 ] = p2;
}

Constraint* ConstraintInfoEqual::createConstraint(SubSystem& subsys) const{
	return new ConstraintEqual(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}

ConstraintType ConstraintInfoEqual::getTypeId() const
{
    return Equal;
}

ConstraintInfo* ConstraintInfoEqual::clone() const
{
    return new ConstraintInfoEqual(*this);
}


// DifferenceInfo
ConstraintInfoDifference::ConstraintInfoDifference(
		double* const p1,
		double* const p2,
		double* const d
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::param1 ] = p1;
    variables[ ConstraintVariables::param2 ] = p2;
    variables[ ConstraintVariables::difference ] = d;
}

Constraint* ConstraintInfoDifference::createConstraint(SubSystem& subsys) const{
	return new ConstraintDifference(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}
ConstraintType ConstraintInfoDifference::getTypeId() const
{
    return Difference;
}

ConstraintInfo* ConstraintInfoDifference::clone() const {
    return new ConstraintInfoDifference( *this );
}

// P2PDistance
ConstraintInfoP2PDistance::ConstraintInfoP2PDistance(
		const Point &p1,
		const Point &p2,
		double* const d
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::p1x ] = p1.x;
    variables[ ConstraintVariables::p1y ] = p1.y;
    variables[ ConstraintVariables::p2x ] = p2.x;
    variables[ ConstraintVariables::p2y ] = p2.y;
    variables[ ConstraintVariables::distance ] = d;
}

Constraint* ConstraintInfoP2PDistance::createConstraint(SubSystem& subsys) const{
	return new ConstraintP2PDistance(
			subsys.getVariables(),
			subsys.getIndices( variables ),
			subsys.getDependentVariableCount(),
			scale_coef );
}
ConstraintType ConstraintInfoP2PDistance::getTypeId() const
{
    return P2PDistance;
}

ConstraintInfo* ConstraintInfoP2PDistance::clone() const
{
    return new ConstraintInfoP2PDistance(*this);
}


// P2PAngleInfo
ConstraintInfoP2PAngle::ConstraintInfoP2PAngle(
		const Point &p1,
		const Point &p2,
		double* const a,
		double da_
):ConstraintInfo(), da(da_){
	variables.resize( variable_count );
    variables[ ConstraintVariables::p1x ] = p1.x;
    variables[ ConstraintVariables::p1y ] = p1.y;
    variables[ ConstraintVariables::p2x ] = p2.x;
    variables[ ConstraintVariables::p2y ] = p2.y;
    variables[ ConstraintVariables::angle ] = a;
}

Constraint* ConstraintInfoP2PAngle::createConstraint(SubSystem& subsys) const{
	return new ConstraintP2PAngle(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef,
			da );
}

ConstraintType ConstraintInfoP2PAngle::getTypeId() const
{
    return P2PAngle;
}

ConstraintInfo* ConstraintInfoP2PAngle::clone() const
{
    return new ConstraintInfoP2PAngle(*this);
}


// P2LDistanceInfo
ConstraintInfoP2LDistance::ConstraintInfoP2LDistance(
		const Point &p,
		const Line &l,
		double* const d
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::px ] = p.x;
    variables[ ConstraintVariables::py ] = p.y;
    variables[ ConstraintVariables::l_p1x ] = l.p1.x;
    variables[ ConstraintVariables::l_p1y ] = l.p1.y;
    variables[ ConstraintVariables::l_p2x ] = l.p2.x;
    variables[ ConstraintVariables::l_p2y ] = l.p2.y;
    variables[ ConstraintVariables::distance ] = d;
}

Constraint* ConstraintInfoP2LDistance::createConstraint(SubSystem& subsys) const{
	return new ConstraintP2LDistance(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}

ConstraintType ConstraintInfoP2LDistance::getTypeId() const
{
    return P2LDistance;
}

ConstraintInfo* ConstraintInfoP2LDistance::clone() const
{
    return new ConstraintInfoP2LDistance(*this);
}

// PointOnLine
ConstraintInfoPointOnLine::ConstraintInfoPointOnLine(
		const Point &p,
		const Line &l
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::px ] = p.x;
    variables[ ConstraintVariables::py ] = p.y;
    variables[ ConstraintVariables::l_p1x ] = l.p1.x;
    variables[ ConstraintVariables::l_p1y ] = l.p1.y;
    variables[ ConstraintVariables::l_p2x ] = l.p2.x;
    variables[ ConstraintVariables::l_p2y ] = l.p2.y;
}

ConstraintInfoPointOnLine::ConstraintInfoPointOnLine(
		const Point &p,
		const Point &lp1,
		const Point &lp2
):ConstraintInfo(){
	variables.resize(variable_count);
    variables[ ConstraintVariables::px ] = p.x;
    variables[ ConstraintVariables::py ] = p.y;
    variables[ ConstraintVariables::l_p1x ] = lp1.x;
    variables[ ConstraintVariables::l_p1y ] = lp1.y;
    variables[ ConstraintVariables::l_p2x ] = lp2.x;
    variables[ ConstraintVariables::l_p2y ] = lp2.y;
}

Constraint* ConstraintInfoPointOnLine::createConstraint(SubSystem& subsys) const{
	return new ConstraintPointOnLine(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}
ConstraintType ConstraintInfoPointOnLine::getTypeId() const
{
    return PointOnLine;
}

ConstraintInfo* ConstraintInfoPointOnLine::clone() const
{
    return new ConstraintInfoPointOnLine(*this);
}

// PointOnPerpBisectorInfo
ConstraintInfoPointOnPerpBisector::ConstraintInfoPointOnPerpBisector(
		const Point &p,
		const Line &l
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::p0x ] = p.x;
    variables[ ConstraintVariables::p0y ] = p.y;
    variables[ ConstraintVariables::p1x ] = l.p1.x;
    variables[ ConstraintVariables::p1y ] = l.p1.y;
    variables[ ConstraintVariables::p2x ] = l.p2.x;
    variables[ ConstraintVariables::p2y ] = l.p2.y;
}

ConstraintInfoPointOnPerpBisector::ConstraintInfoPointOnPerpBisector(
		const Point &p,
		const Point &lp1,
		const Point &lp2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::p0x ] = p.x;
    variables[ ConstraintVariables::p0y ] = p.y;
    variables[ ConstraintVariables::p1x ] = lp1.x;
    variables[ ConstraintVariables::p1y ] = lp1.y;
    variables[ ConstraintVariables::p2x ] = lp2.x;
    variables[ ConstraintVariables::p2y ] = lp2.y;
}

Constraint* ConstraintInfoPointOnPerpBisector::createConstraint(SubSystem& subsys) const{
	return new ConstraintPointOnPerpBisector(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}
ConstraintType ConstraintInfoPointOnPerpBisector::getTypeId() const
{
    return PointOnPerpBisector;
}

ConstraintInfo* ConstraintInfoPointOnPerpBisector::clone() const
{
    return new ConstraintInfoPointOnPerpBisector(*this);
}

// ParallelInfo
ConstraintInfoParallel::ConstraintInfoParallel(
		const Line &l1,
		const Line &l2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1.p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1.p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1.p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1.p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2.p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2.p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2.p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2.p2.y;
}

Constraint* ConstraintInfoParallel::createConstraint(SubSystem& subsys) const{
	return new ConstraintParallel(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}
ConstraintType ConstraintInfoParallel::getTypeId() const
{
    return Parallel;
}

ConstraintInfo* ConstraintInfoParallel::clone() const
{
    return new ConstraintInfoParallel(*this);
}

// Perpendicular
ConstraintInfoPerpendicular::ConstraintInfoPerpendicular(
		const Line &l1,
		const Line &l2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1.p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1.p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1.p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1.p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2.p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2.p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2.p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2.p2.y;
}


ConstraintInfoPerpendicular::ConstraintInfoPerpendicular(
		const Point &l1p1,
		const Point &l1p2,
		const Point &l2p1,
		const Point &l2p2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2p2.y;
}

Constraint* ConstraintInfoPerpendicular::createConstraint(SubSystem& subsys) const{
	return new ConstraintPerpendicular(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}

ConstraintType ConstraintInfoPerpendicular::getTypeId() const
{
    return Perpendicular;
}

ConstraintInfo* ConstraintInfoPerpendicular::clone() const
{
    return new ConstraintInfoPerpendicular(*this);
}

// L2LAngleInfo
ConstraintInfoL2LAngle::ConstraintInfoL2LAngle(
		const Line &l1,
		const Line &l2,
		double* const a
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1.p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1.p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1.p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1.p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2.p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2.p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2.p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2.p2.y;
    variables[ ConstraintVariables::angle ] = a;
}

ConstraintInfoL2LAngle::ConstraintInfoL2LAngle(
		const Point &l1p1,
		const Point &l1p2,
		const Point &l2p1,
		const Point &l2p2,
		double* const a
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2p2.y;
    variables[ ConstraintVariables::angle ] = a;
}

Constraint* ConstraintInfoL2LAngle::createConstraint(SubSystem& subsys) const{
	return new ConstraintL2LAngle(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}

ConstraintType ConstraintInfoL2LAngle::getTypeId() const
{
    return L2LAngle;
}

ConstraintInfo* ConstraintInfoL2LAngle::clone() const
{
    return new ConstraintInfoL2LAngle(*this);
}


// MidpointOnLineInfo
ConstraintInfoMidpointOnLine::ConstraintInfoMidpointOnLine(
		const Line &l1,
		const Line &l2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1.p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1.p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1.p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1.p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2.p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2.p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2.p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2.p2.y;
}

ConstraintInfoMidpointOnLine::ConstraintInfoMidpointOnLine(
		const Point &l1p1,
		const Point &l1p2,
		const Point &l2p1,
		const Point &l2p2
):ConstraintInfo(){
	variables.resize( variable_count );
    variables[ ConstraintVariables::l1p1x ] = l1p1.x;
    variables[ ConstraintVariables::l1p1y ] = l1p1.y;
    variables[ ConstraintVariables::l1p2x ] = l1p2.x;
    variables[ ConstraintVariables::l1p2y ] = l1p2.y;
    variables[ ConstraintVariables::l2p1x ] = l2p1.x;
    variables[ ConstraintVariables::l2p1y ] = l2p1.y;
    variables[ ConstraintVariables::l2p2x ] = l2p2.x;
    variables[ ConstraintVariables::l2p2y ] = l2p2.y;
}

Constraint* ConstraintInfoMidpointOnLine::createConstraint(SubSystem& subsys) const{
	return new ConstraintMidpointOnLine(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef );
}

ConstraintType ConstraintInfoMidpointOnLine::getTypeId() const
{
    return MidpointOnLine;
}

ConstraintInfo* ConstraintInfoMidpointOnLine::clone() const
{
    return new ConstraintInfoMidpointOnLine(*this);
}

// TangentCircumfInfo
ConstraintInfoTangentCircumf::ConstraintInfoTangentCircumf(
		const Point &p1,
		const Point &p2,
		double* const rad1,
		double* const rad2,
		bool internal_
):ConstraintInfo(), internal(internal_){
	variables.resize( variable_count );
    variables[ ConstraintVariables::c1x ] = p1.x;
    variables[ ConstraintVariables::c1y ] = p1.y;
    variables[ ConstraintVariables::c2x ] = p2.x;
    variables[ ConstraintVariables::c2y ] = p2.y;
    variables[ ConstraintVariables::r1 ] = rad1;
    variables[ ConstraintVariables::r2 ] = rad2;
}

Constraint* ConstraintInfoTangentCircumf::createConstraint(SubSystem& subsys) const{
	return new ConstraintTangentCircumf(
			subsys.getVariables(),
			subsys.getIndices(variables),
			subsys.getDependentVariableCount(),
			scale_coef,
			internal );
}

ConstraintType ConstraintInfoTangentCircumf::getTypeId() const
{
    return TangentCircumf;
}

ConstraintInfo* ConstraintInfoTangentCircumf::clone() const
{
    return new ConstraintInfoTangentCircumf(*this);
}


} //namespace GCS
