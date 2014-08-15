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

#ifndef FREEGCS_CONSTRAINTS_EXP_H
#define FREEGCS_CONSTRAINTS_EXP_H

#include "Geo_exp.h"
#include "Util_exp.h"

namespace GCS_EXP
{

    ///////////////////////////////////////
    // Constraints
    ///////////////////////////////////////

    enum ConstraintType {
        None = 0,
        Equal = 1,
        Difference = 2,
        P2PDistance = 3,
        P2PAngle = 4,
        P2LDistance = 5,
        PointOnLine = 6,
        PointOnPerpBisector = 7,
        Parallel = 8,
        Perpendicular = 9,
        L2LAngle = 10,
        MidpointOnLine = 11,
        TangentCircumf = 12
    };

    template <ConstraintType CT>
    class ConstraintVariables;

    template <>
    struct ConstraintVariables<Equal> {
    	enum Variables{ param1, param2, variable_count }; };

    template <>
    struct ConstraintVariables<Difference> {
    	enum Variables{ param1, param2, difference, variable_count }; };

    template <>
    struct ConstraintVariables<P2PDistance> {
    	enum Variables{ p1x, p1y, p2x, p2y, distance, variable_count }; };

    template <>
    struct ConstraintVariables<P2PAngle> {
    	enum Variables{ p1x, p1y, p2x, p2y, angle, variable_count }; };

    template <>
    struct ConstraintVariables<P2LDistance> {
    	enum Variables{ px, py, l_p1x, l_p1y, l_p2x, l_p2y, distance, variable_count }; };

    template <>
    struct ConstraintVariables<PointOnLine> {
    	enum Variables{ px, py, l_p1x, l_p1y, l_p2x, l_p2y, variable_count }; };

    template <>
    struct ConstraintVariables<PointOnPerpBisector> {
    	enum Variables{ p0x, p0y, p1x, p1y, p2x, p2y, variable_count }; };

    template <>
    struct ConstraintVariables<Parallel> {
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, variable_count }; };

    template <>
    struct ConstraintVariables<Perpendicular> {
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, variable_count }; };

    template <>
    struct ConstraintVariables<L2LAngle> {
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, angle, variable_count }; };

    template <>
    struct ConstraintVariables<MidpointOnLine> {
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, variable_count }; };

    template <>
    struct ConstraintVariables<TangentCircumf> {
    	enum Variables{ c1x, c1y, c2x, c2y, r1, r2, variable_count }; };




    class Constraint
    {
    public:
    	typedef std::pair<variable_index_type,double> grad_component_t;
    protected:
    	const std::vector<double>& variables;
    	std::vector<index_type> variable_indices;
        variable_index_type dependent_variable_count;
        double scale;
        int tag;

        template <int i> index_type index() const {
        	return variable_indices[i]; }
        template <int i> double value() const {
        	return variables[index<i>()]; }
        template <int i> bool is_dependent() const {
        	return index<i>() < dependent_variable_count; }
        template <int i> void dependent_insert( std::vector< grad_component_t >& gradVec, double value ) const {
        	if( is_dependent<i>() ) gradVec.push_back( grad_component_t( index<i>(), value ) ); }
    public:
        Constraint(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count );
        virtual ~Constraint(){}

        inline const std::vector<index_type>& indices() const { return variable_indices; }

        void setTag(int tagId) { tag = tagId; }
        int getTag() const { return tag; }


        virtual ConstraintType getTypeId() const = 0;
        virtual Constraint* clone() const  = 0;
        virtual void rescale(double coef=1.);
        virtual double error() = 0;
        virtual double grad(index_type) = 0;
        virtual void grad( std::vector< grad_component_t >& gradVec ) = 0;
        virtual double maxStep(const std::vector<double>& dir, double lim=1.);
    };

    // Equal
    class ConstraintEqual : public Constraint, protected ConstraintVariables<Equal>
    {
    public:
        ConstraintEqual(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef  );
        virtual ConstraintType getTypeId() const;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // Difference
    class ConstraintDifference : public Constraint, protected ConstraintVariables<Difference>
    {
    public:
        ConstraintDifference(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // P2PDistance
    class ConstraintP2PDistance : public Constraint, protected ConstraintVariables<P2PDistance>
    {
    public:
        ConstraintP2PDistance(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
        virtual double maxStep( const std::vector<double>& dir, double lim=1.);
    };

    // P2PAngle
    class ConstraintP2PAngle : public Constraint, protected ConstraintVariables<P2PAngle>
    {
    private:
        double da;
    public:
        ConstraintP2PAngle(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef,
        		double da_/*=0.*/);
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
        virtual double maxStep(const std::vector<double>& dir, double lim=1.);
    };

    // P2LDistance
    class ConstraintP2LDistance : public Constraint, protected ConstraintVariables<P2LDistance>
    {
    public:
        ConstraintP2LDistance(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
        virtual double maxStep(const std::vector<double>& dir, double lim=1.);
    };

    // PointOnLine
    class ConstraintPointOnLine : public Constraint, ConstraintVariables<PointOnLine>
    {
    public:
        ConstraintPointOnLine(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // PointOnPerpBisector
    class ConstraintPointOnPerpBisector : public Constraint, protected ConstraintVariables<PointOnPerpBisector>
    {
    public:
        ConstraintPointOnPerpBisector(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // Parallel
    class ConstraintParallel : public Constraint, protected ConstraintVariables<Parallel>
    {
    public:
        ConstraintParallel(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // Perpendicular
    class ConstraintPerpendicular : public Constraint, protected ConstraintVariables<Perpendicular>
    {
    public:
        ConstraintPerpendicular(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // L2LAngle
    class ConstraintL2LAngle : public Constraint, protected ConstraintVariables<L2LAngle>
    {
    public:
        ConstraintL2LAngle(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
        virtual double maxStep( const std::vector<double>& dir, double lim=1.);
    };

    // MidpointOnLine
    class ConstraintMidpointOnLine : public Constraint, protected ConstraintVariables<MidpointOnLine>
    {
    public:
        ConstraintMidpointOnLine(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef );
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };

    // TangentCircumf
    class ConstraintTangentCircumf : public Constraint, protected ConstraintVariables<TangentCircumf>
    {
    private:
        bool internal;
    public:
        ConstraintTangentCircumf(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count,
        		double scale_coef,
        		bool internal_/*=false*/);
        virtual ConstraintType getTypeId() const ;
        virtual Constraint* clone() const;
        virtual void rescale(double coef=1.);
        virtual double error();
        virtual double grad(index_type);
        virtual void grad( std::vector< grad_component_t >& gradVec );
    };



} //namespace GCS_EXP

#endif // FREEGCS_CONSTRAINTS_EXP_H
