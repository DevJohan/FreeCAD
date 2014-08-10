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





    class Constraint
    {
    protected:
    	const std::vector<double>& variables;
    	std::vector<index_type> paramenters_indices;
        variable_index_type dependent_variable_count;
        double scale;
        int tag;

        template <int i> index_type index() const { return paramenters_indices[i]; }
        template <int i> double value() const { return variables[index<i>()]; }
        template <int i> bool is_dependent() const { return index<i>() < dependent_variable_count; }
    public:
        Constraint(
        		const std::vector<double>& parameters,
        		const std::vector<index_type> indices,
        		index_type dependent_var_count );
        virtual ~Constraint(){}

        inline const std::vector<index_type>& params() const { return paramenters_indices; }

        void setTag(int tagId) { tag = tagId; }
        int getTag() const { return tag; }


        virtual ConstraintType getTypeId() const = 0;
        virtual Constraint* clone() const  = 0;
        virtual void rescale(double coef=1.);
        virtual double error() = 0;
        virtual double grad(index_type) = 0;
        // virtual void grad(MAP_pD_D &deriv);  --> TODO: vectorized grad version
        virtual double maxStep(const std::vector<double>& dir, double lim=1.);
    };

    // Equal
    class ConstraintEqual : public Constraint
    {
    private:
    	enum Variables{ param1, param2, variable_count };
//        inline index_type param1_i() const { return paramenters_indices[0]; }
//        inline index_type param2_i() const { return paramenters_indices[1]; }
//
//        inline double value_param1() const { return variables[param1_i()]; }
//        inline double value_param2() const { return variables[param2_i()]; }
//
//        inline bool param1_dependent() const { return is_dependent<param1>(); }
//        inline bool param2_dependent() const { return is_dependent<param2>(); }
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
    };

    // Difference
    class ConstraintDifference : public Constraint
    {
    private:
    	enum Variables{ param1, param2, difference, variable_count };
//        inline index_type param1_i() const { return paramenters_indices[0]; }
//        inline index_type param2_i() const { return paramenters_indices[1]; }
//        inline index_type difference_i() const { return paramenters_indices[2]; }
//
//        inline double value_param1() const { return variables[param1_i()]; }
//        inline double value_param2() const { return variables[param2_i()]; }
//        inline double value_difference() const { return variables[difference_i()]; }
//
//        inline bool param1_dependent() const { return isDependentVariable( param1_i() ); }
//        inline bool param2_dependent() const { return isDependentVariable( param2_i() ); }
//        inline bool difference_dependent() const { return isDependentVariable( difference_i() ); }
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
    };

    // P2PDistance
    class ConstraintP2PDistance : public Constraint
    {
    private:
    	enum Variables{ p1x, p1y, p2x, p2y, distance, variable_count };
//        inline index_type p1x_i() const   { return paramenters_indices[0]; }
//        inline index_type p1y_i() const   { return paramenters_indices[1]; }
//        inline index_type p2x_i() const   { return paramenters_indices[2]; }
//        inline index_type p2y_i() const   { return paramenters_indices[3]; }
//        inline index_type distance_i() const { return paramenters_indices[4]; }
//
//        inline double value_p1x() const   { return variables[p1x_i()]; }
//        inline double value_p1y() const   { return variables[p1y_i()]; }
//        inline double value_p2x() const   { return variables[p2x_i()]; }
//        inline double value_p2y() const   { return variables[p2y_i()]; }
//        inline double value_distance() const { return variables[distance_i()]; }
//
//        inline bool p1x_dependent() const { return isDependentVariable( p1x_i() ); }
//        inline bool p1y_dependent() const { return isDependentVariable( p1y_i() ); }
//        inline bool p2x_dependent() const { return isDependentVariable( p2x_i() ); }
//        inline bool p2y_dependent() const { return isDependentVariable( p2y_i() ); }
//        inline bool distance_dependent() const { return isDependentVariable( distance_i() ); }
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
        virtual double maxStep( const std::vector<double>& dir, double lim=1.);
    };

    // P2PAngle
    class ConstraintP2PAngle : public Constraint
    {
    private:
    	enum Variables{ p1x, p1y, p2x, p2y, angle, variable_count };
//        inline index_type p1x_i() const   { return paramenters_indices[0]; }
//        inline index_type p1y_i() const   { return paramenters_indices[1]; }
//        inline index_type p2x_i() const   { return paramenters_indices[2]; }
//        inline index_type p2y_i() const   { return paramenters_indices[3]; }
//        inline index_type angle_i() const { return paramenters_indices[4]; }
//
//        inline double value_p1x() const   { return variables[p1x_i()]; }
//        inline double value_p1y() const   { return variables[p1y_i()]; }
//        inline double value_p2x() const   { return variables[p2x_i()]; }
//        inline double value_p2y() const   { return variables[p2y_i()]; }
//        inline double value_angle() const { return variables[angle_i()]; }
//
//        inline bool p1x_dependent() const { return isDependentVariable( p1x_i() ); }
//        inline bool p1y_dependent() const { return isDependentVariable( p1y_i() ); }
//        inline bool p2x_dependent() const { return isDependentVariable( p2x_i() ); }
//        inline bool p2y_dependent() const { return isDependentVariable( p2y_i() ); }
//        inline bool angle_dependent() const { return isDependentVariable( angle_i() ); }

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
        virtual double maxStep(const std::vector<double>& dir, double lim=1.);
    };

    // P2LDistance
    class ConstraintP2LDistance : public Constraint
    {
    private:
    	enum Variables{ px, py, l_p1x, l_p1y, l_p2x, l_p2y, distance, variable_count };
//        inline index_type p0x_i() const { return paramenters_indices[0]; }
//        inline index_type p0y_i() const { return paramenters_indices[1]; }
//        inline index_type p1x_i() const { return paramenters_indices[2]; }
//        inline index_type p1y_i() const { return paramenters_indices[3]; }
//        inline index_type p2x_i() const { return paramenters_indices[4]; }
//        inline index_type p2y_i() const { return paramenters_indices[5]; }
//        inline index_type distance_i() const { return paramenters_indices[6]; }
//
//        inline double value_p0x() const { return variables[p0x_i()]; }
//        inline double value_p0y() const { return variables[p0y_i()]; }
//        inline double value_p1x() const { return variables[p1x_i()]; }
//        inline double value_p1y() const { return variables[p1y_i()]; }
//        inline double value_p2x() const { return variables[p2x_i()]; }
//        inline double value_p2y() const { return variables[p2y_i()]; }
//        inline double value_distance() const { return variables[distance_i()]; }
//
//        inline bool p0x_dependent() const { return isDependentVariable( p0x_i() ); }
//        inline bool p0y_dependent() const { return isDependentVariable( p0y_i() ); }
//        inline bool p1x_dependent() const { return isDependentVariable( p1x_i() ); }
//        inline bool p1y_dependent() const { return isDependentVariable( p1y_i() ); }
//        inline bool p2x_dependent() const { return isDependentVariable( p2x_i() ); }
//        inline bool p2y_dependent() const { return isDependentVariable( p2y_i() ); }
//        inline bool distance_dependent() const { return isDependentVariable( distance_i() ); }
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
        virtual double maxStep(const std::vector<double>& dir, double lim=1.);
    };

    // PointOnLine
    class ConstraintPointOnLine : public Constraint
    {
    private:
    	enum Variables{ px, py, l_p1x, l_p1y, l_p2x, l_p2y, variable_count };
//        inline index_type p0x_i() const { return paramenters_indices[0]; }
//        inline index_type p0y_i() const { return paramenters_indices[1]; }
//        inline index_type p1x_i() const { return paramenters_indices[2]; }
//        inline index_type p1y_i() const { return paramenters_indices[3]; }
//        inline index_type p2x_i() const { return paramenters_indices[4]; }
//        inline index_type p2y_i() const { return paramenters_indices[5]; }
//
//        inline double value_p0x() const { return variables[p0x_i()]; }
//        inline double value_p0y() const { return variables[p0y_i()]; }
//        inline double value_p1x() const { return variables[p1x_i()]; }
//        inline double value_p1y() const { return variables[p1y_i()]; }
//        inline double value_p2x() const { return variables[p2x_i()]; }
//        inline double value_p2y() const { return variables[p2y_i()]; }
//
//        inline bool p0x_dependent() const { return isDependentVariable( p0x_i() ); }
//        inline bool p0y_dependent() const { return isDependentVariable( p0y_i() ); }
//        inline bool p1x_dependent() const { return isDependentVariable( p1x_i() ); }
//        inline bool p1y_dependent() const { return isDependentVariable( p1y_i() ); }
//        inline bool p2x_dependent() const { return isDependentVariable( p2x_i() ); }
//        inline bool p2y_dependent() const { return isDependentVariable( p2y_i() ); }
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
    };

    // PointOnPerpBisector
    class ConstraintPointOnPerpBisector : public Constraint
    {
    private:
    	enum Variables{ p0x, p0y, p1x, p1y, p2x, p2y, variable_count };
//        inline index_type p0x_i() const { return paramenters_indices[0]; }
//        inline index_type p0y_i() const { return paramenters_indices[1]; }
//        inline index_type p1x_i() const { return paramenters_indices[2]; }
//        inline index_type p1y_i() const { return paramenters_indices[3]; }
//        inline index_type p2x_i() const { return paramenters_indices[4]; }
//        inline index_type p2y_i() const { return paramenters_indices[5]; }
//
//        inline double value_p0x() const { return variables[p0x_i()]; }
//        inline double value_p0y() const { return variables[p0y_i()]; }
//        inline double value_p1x() const { return variables[p1x_i()]; }
//        inline double value_p1y() const { return variables[p1y_i()]; }
//        inline double value_p2x() const { return variables[p2x_i()]; }
//        inline double value_p2y() const { return variables[p2y_i()]; }
//
//        inline bool p0x_dependent() const { return isDependentVariable( p0x_i() ); }
//        inline bool p0y_dependent() const { return isDependentVariable( p0y_i() ); }
//        inline bool p1x_dependent() const { return isDependentVariable( p1x_i() ); }
//        inline bool p1y_dependent() const { return isDependentVariable( p1y_i() ); }
//        inline bool p2x_dependent() const { return isDependentVariable( p2x_i() ); }
//        inline bool p2y_dependent() const { return isDependentVariable( p2y_i() ); }
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
    };

    // Parallel
    class ConstraintParallel : public Constraint
    {
    private:
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, variable_count };
//        inline index_type l1p1x_i() const { return paramenters_indices[0]; }
//        inline index_type l1p1y_i() const { return paramenters_indices[1]; }
//        inline index_type l1p2x_i() const { return paramenters_indices[2]; }
//        inline index_type l1p2y_i() const { return paramenters_indices[3]; }
//        inline index_type l2p1x_i() const { return paramenters_indices[4]; }
//        inline index_type l2p1y_i() const { return paramenters_indices[5]; }
//        inline index_type l2p2x_i() const { return paramenters_indices[6]; }
//        inline index_type l2p2y_i() const { return paramenters_indices[7]; }
//
//        inline double value_l1p1x() const { return variables[l1p1x_i()]; }
//        inline double value_l1p1y() const { return variables[l1p1y_i()]; }
//        inline double value_l1p2x() const { return variables[l1p2x_i()]; }
//        inline double value_l1p2y() const { return variables[l1p2y_i()]; }
//        inline double value_l2p1x() const { return variables[l2p1x_i()]; }
//        inline double value_l2p1y() const { return variables[l2p1y_i()]; }
//        inline double value_l2p2x() const { return variables[l2p2x_i()]; }
//        inline double value_l2p2y() const { return variables[l2p2y_i()]; }
//
//        inline bool l1p1x_dependent() const { return isDependentVariable( l1p1x_i() ); }
//        inline bool l1p1y_dependent() const { return isDependentVariable( l1p1y_i() ); }
//        inline bool l1p2x_dependent() const { return isDependentVariable( l1p2x_i() ); }
//        inline bool l1p2y_dependent() const { return isDependentVariable( l1p2y_i() ); }
//        inline bool l2p1x_dependent() const { return isDependentVariable( l2p1x_i() ); }
//        inline bool l2p1y_dependent() const { return isDependentVariable( l2p1y_i() ); }
//        inline bool l2p2x_dependent() const { return isDependentVariable( l2p2x_i() ); }
//        inline bool l2p2y_dependent() const { return isDependentVariable( l2p2y_i() ); }
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
    };

    // Perpendicular
    class ConstraintPerpendicular : public Constraint
    {
    private:
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, variable_count };
//        inline index_type l1p1x_i() const { return paramenters_indices[0]; }
//        inline index_type l1p1y_i() const { return paramenters_indices[1]; }
//        inline index_type l1p2x_i() const { return paramenters_indices[2]; }
//        inline index_type l1p2y_i() const { return paramenters_indices[3]; }
//        inline index_type l2p1x_i() const { return paramenters_indices[4]; }
//        inline index_type l2p1y_i() const { return paramenters_indices[5]; }
//        inline index_type l2p2x_i() const { return paramenters_indices[6]; }
//        inline index_type l2p2y_i() const { return paramenters_indices[7]; }
//
//        inline double value_l1p1x() const { return variables[l1p1x_i()]; }
//        inline double value_l1p1y() const { return variables[l1p1y_i()]; }
//        inline double value_l1p2x() const { return variables[l1p2x_i()]; }
//        inline double value_l1p2y() const { return variables[l1p2y_i()]; }
//        inline double value_l2p1x() const { return variables[l2p1x_i()]; }
//        inline double value_l2p1y() const { return variables[l2p1y_i()]; }
//        inline double value_l2p2x() const { return variables[l2p2x_i()]; }
//        inline double value_l2p2y() const { return variables[l2p2y_i()]; }
//
//        inline bool l1p1x_dependent() const { return isDependentVariable( l1p1x_i() ); }
//        inline bool l1p1y_dependent() const { return isDependentVariable( l1p1y_i() ); }
//        inline bool l1p2x_dependent() const { return isDependentVariable( l1p2x_i() ); }
//        inline bool l1p2y_dependent() const { return isDependentVariable( l1p2y_i() ); }
//        inline bool l2p1x_dependent() const { return isDependentVariable( l2p1x_i() ); }
//        inline bool l2p1y_dependent() const { return isDependentVariable( l2p1y_i() ); }
//        inline bool l2p2x_dependent() const { return isDependentVariable( l2p2x_i() ); }
//        inline bool l2p2y_dependent() const { return isDependentVariable( l2p2y_i() ); }
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
    };

    // L2LAngle
    class ConstraintL2LAngle : public Constraint
    {
    private:
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, angle, variable_count };
//        inline index_type l1p1x_i() const { return paramenters_indices[0]; }
//        inline index_type l1p1y_i() const { return paramenters_indices[1]; }
//        inline index_type l1p2x_i() const { return paramenters_indices[2]; }
//        inline index_type l1p2y_i() const { return paramenters_indices[3]; }
//        inline index_type l2p1x_i() const { return paramenters_indices[4]; }
//        inline index_type l2p1y_i() const { return paramenters_indices[5]; }
//        inline index_type l2p2x_i() const { return paramenters_indices[6]; }
//        inline index_type l2p2y_i() const { return paramenters_indices[7]; }
//        inline index_type angle_i() const { return paramenters_indices[8]; }
//
//        inline double value_l1p1x() const { return variables[l1p1x_i()]; }
//        inline double value_l1p1y() const { return variables[l1p1y_i()]; }
//        inline double value_l1p2x() const { return variables[l1p2x_i()]; }
//        inline double value_l1p2y() const { return variables[l1p2y_i()]; }
//        inline double value_l2p1x() const { return variables[l2p1x_i()]; }
//        inline double value_l2p1y() const { return variables[l2p1y_i()]; }
//        inline double value_l2p2x() const { return variables[l2p2x_i()]; }
//        inline double value_l2p2y() const { return variables[l2p2y_i()]; }
//        inline double value_angle() const { return variables[angle_i()]; }
//
//        inline bool l1p1x_dependent() const { return isDependentVariable( l1p1x_i() ); }
//        inline bool l1p1y_dependent() const { return isDependentVariable( l1p1y_i() ); }
//        inline bool l1p2x_dependent() const { return isDependentVariable( l1p2x_i() ); }
//        inline bool l1p2y_dependent() const { return isDependentVariable( l1p2y_i() ); }
//        inline bool l2p1x_dependent() const { return isDependentVariable( l2p1x_i() ); }
//        inline bool l2p1y_dependent() const { return isDependentVariable( l2p1y_i() ); }
//        inline bool l2p2x_dependent() const { return isDependentVariable( l2p2x_i() ); }
//        inline bool l2p2y_dependent() const { return isDependentVariable( l2p2y_i() ); }
//        inline bool angle_dependent() const { return isDependentVariable( angle_i() ); }
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
        virtual double maxStep( const std::vector<double>& dir, double lim=1.);
    };

    // MidpointOnLine
    class ConstraintMidpointOnLine : public Constraint
    {
    private:
    	enum Variables{ l1p1x, l1p1y, l1p2x, l1p2y, l2p1x, l2p1y, l2p2x, l2p2y, variable_count };
//        inline index_type l1p1x_i() const { return paramenters_indices[0]; }
//        inline index_type l1p1y_i() const { return paramenters_indices[1]; }
//        inline index_type l1p2x_i() const { return paramenters_indices[2]; }
//        inline index_type l1p2y_i() const { return paramenters_indices[3]; }
//        inline index_type l2p1x_i() const { return paramenters_indices[4]; }
//        inline index_type l2p1y_i() const { return paramenters_indices[5]; }
//        inline index_type l2p2x_i() const { return paramenters_indices[6]; }
//        inline index_type l2p2y_i() const { return paramenters_indices[7]; }
//
//        inline double value_l1p1x() const { return variables[l1p1x_i()]; }
//        inline double value_l1p1y() const { return variables[l1p1y_i()]; }
//        inline double value_l1p2x() const { return variables[l1p2x_i()]; }
//        inline double value_l1p2y() const { return variables[l1p2y_i()]; }
//        inline double value_l2p1x() const { return variables[l2p1x_i()]; }
//        inline double value_l2p1y() const { return variables[l2p1y_i()]; }
//        inline double value_l2p2x() const { return variables[l2p2x_i()]; }
//        inline double value_l2p2y() const { return variables[l2p2y_i()]; }
//
//        inline bool l1p1x_dependent() const ;// { return variables[l1p1x_i()]; }
//        inline bool l1p1y_dependent() const ;// { return variables[l1p1y_i()]; }
//        inline bool l1p2x_dependent() const ;// { return variables[l1p2x_i()]; }
//        inline bool l1p2y_dependent() const ;// { return variables[l1p2y_i()]; }
//        inline bool l2p1x_dependent() const ;// { return variables[l2p1x_i()]; }
//        inline bool l2p1y_dependent() const ;// { return variables[l2p1y_i()]; }
//        inline bool l2p2x_dependent() const ;// { return variables[l2p2x_i()]; }
//        inline bool l2p2y_dependent() const ;// { return variables[l2p2y_i()]; }
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
    };

    // TangentCircumf
    class ConstraintTangentCircumf : public Constraint
    {
    private:
    	enum Variables{ c1x, c1y, c2x, c2y, r1, r2, variable_count };
//        inline index_type c1x_i() const { return paramenters_indices[0]; }
//        inline index_type c1y_i() const { return paramenters_indices[1]; }
//        inline index_type c2x_i() const { return paramenters_indices[2]; }
//        inline index_type c2y_i() const { return paramenters_indices[3]; }
//        inline index_type r1_i() const  { return paramenters_indices[4]; }
//        inline index_type r2_i() const  { return paramenters_indices[5]; }
//
//        inline double value_c1x() const { return variables[c1x_i()]; }
//        inline double value_c1y() const { return variables[c1y_i()]; }
//        inline double value_c2x() const { return variables[c2x_i()]; }
//        inline double value_c2y() const { return variables[c2y_i()]; }
//        inline double value_r1() const  { return variables[r1_i()]; }
//        inline double value_r2() const  { return variables[r2_i()]; }
//
//        inline bool c1x_dependent() const ;// { return variables[c1x_i()]; }
//        inline bool c1y_dependent() const ;// { return variables[c1y_i()]; }
//        inline bool c2x_dependent() const ;// { return variables[c2x_i()]; }
//        inline bool c2y_dependent() const ;// { return variables[c2y_i()]; }
//        inline bool r1_dependent() const ;//  { return variables[r1_i()]; }
//        inline bool r2_dependent() const ;//  { return variables[r2_i()]; }
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
    };



} //namespace GCS_EXP

#endif // FREEGCS_CONSTRAINTS_EXP_H
