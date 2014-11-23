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

#ifndef FREEGCS_CONSTRAINTSINFO_EXP_H
#define FREEGCS_CONSTRAINTSINFO_EXP_H

#include "Constraints_exp.h"

namespace GCS_EXP
{
	class SubSystem;

    ///////////////////////////////////////
    // ConstraintsInfo
    ///////////////////////////////////////

    class ConstraintInfo
    {
    protected:
    	std::vector<double*> variables;
        int tag;
        double scale_coef;
    public:
        ConstraintInfo( );
        virtual ~ConstraintInfo(){}

        inline const std::vector<double*>& params() const { return variables; }

        void setTag(int tagId) { tag = tagId; }
        int getTag() const { return tag; }
        bool isPriorityConstraint() const { return tag >= 0; }

        void setScaleCoeff( double coeff ){ scale_coef = coeff; }
        double getScaleCoeff() const { return scale_coef; }

        bool operator< ( const ConstraintInfo& ) const;

        virtual Constraint* createConstraint( SubSystem& sub_system ) const = 0;
        virtual ConstraintType getTypeId() const;
        virtual ConstraintInfo* clone() const  = 0;
    };

    // EqualInfo
    class ConstraintInfoEqual : public ConstraintInfo, protected ConstraintVariables<Equal>
    {
    public:
        ConstraintInfoEqual( double* const p1, double* const p2 );
        virtual Constraint* createConstraint( SubSystem& sub_system ) const;
        virtual ConstraintType getTypeId() const;
        virtual ConstraintInfo* clone() const;
    };

    // DifferenceInfo
    class ConstraintInfoDifference : public ConstraintInfo, protected ConstraintVariables<Difference>
    {
    public:
        ConstraintInfoDifference( double* const p1, double* const p2, double* const d);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // P2PDistanceInfo
    class ConstraintInfoP2PDistance : public ConstraintInfo, protected ConstraintVariables<P2PDistance>
    {
    public:
        ConstraintInfoP2PDistance(const Point &p1, const Point &p2, double* const d );
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // P2PAngleInfo
    class ConstraintInfoP2PAngle : public ConstraintInfo, protected ConstraintVariables<P2PAngle>
    {
    private:
        double da;
    public:
        ConstraintInfoP2PAngle(
        		const Point &p1,
        		const Point &p2,
        		double* const a,
        		double da_=0.);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // P2LDistanceInfo
    class ConstraintInfoP2LDistance : public ConstraintInfo, protected ConstraintVariables<P2LDistance>
    {
    public:
        ConstraintInfoP2LDistance(
        		const Point &p,
        		const Line &l,
        		double* const d);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // PointOnLineInfo
    class ConstraintInfoPointOnLine : public ConstraintInfo, protected ConstraintVariables<PointOnLine>
    {
    public:
        ConstraintInfoPointOnLine(
        		const Point &p,
        		const Line &l);
        ConstraintInfoPointOnLine(
        		const Point &p,
        		const Point &lp1,
        		const Point &lp2);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // PointOnPerpBisector
    class ConstraintInfoPointOnPerpBisector : public ConstraintInfo, protected ConstraintVariables<PointOnPerpBisector>
    {
    public:
        ConstraintInfoPointOnPerpBisector(
        		const Point &p,
        		const Line &l);
        ConstraintInfoPointOnPerpBisector(
        		const Point &p,
        		const Point &lp1,
        		const Point &lp2);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // ParallelInfo
    class ConstraintInfoParallel : public ConstraintInfo, protected ConstraintVariables<Parallel>
    {
    public:
        ConstraintInfoParallel(
        		const Line &l1,
        		const Line &l2);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // PerpendicularInfo
    class ConstraintInfoPerpendicular : public ConstraintInfo, protected ConstraintVariables<Perpendicular>
    {
    public:
        ConstraintInfoPerpendicular(
        		const Line &l1,
        		const Line &l2);
        ConstraintInfoPerpendicular(
        		const Point &l1p1,
        		const Point &l1p2,
        		const Point &l2p1,
        		const Point &l2p2);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // L2LAngleInfo
    class ConstraintInfoL2LAngle : public ConstraintInfo, protected ConstraintVariables<L2LAngle>
    {
    public:
        ConstraintInfoL2LAngle(
        		const Line &l1,
        		const Line &l2,
        		double* const a);
        ConstraintInfoL2LAngle(
        		const Point &l1p1,
        		const Point &l1p2,
        		const Point &l2p1,
        		const Point &l2p2,
        		double* const a);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // MidpointOnLineInfo
    class ConstraintInfoMidpointOnLine : public ConstraintInfo, protected ConstraintVariables<MidpointOnLine>
    {
    public:
        ConstraintInfoMidpointOnLine(
        		const Line &l1,
        		const Line &l2);
        ConstraintInfoMidpointOnLine(
        		const Point &l1p1,
        		const Point &l1p2,
        		const Point &l2p1,
        		const Point &l2p2);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    // TangentCircumfInfo
    class ConstraintInfoTangentCircumf : public ConstraintInfo, protected ConstraintVariables<TangentCircumf>
    {
    private:
    	bool internal;
    public:
        ConstraintInfoTangentCircumf(
        		const Point &p1,
        		const Point &p2,
        		double* const rd1,
        		double* const rd2,
        		bool internal_=false);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };


    // PointOnEllipse
    class ConstraintInfoPointOnEllipse: public ConstraintInfo, protected ConstraintVariables<PointOnEllipse>
    {
    public:
        ConstraintInfoPointOnEllipse(Point &p, Ellipse &e);
        ConstraintInfoPointOnEllipse(Point &p, ArcOfEllipse &a);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    class ConstraintInfoEllipseTangentLine: public ConstraintInfo, protected ConstraintVariables<TangentEllipseLine>
    {
    public:
        ConstraintInfoEllipseTangentLine(Line &l, Ellipse &e);
        ConstraintInfoEllipseTangentLine(Line &l, ArcOfEllipse &a);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    class ConstraintInfoInternalAlignmentPoint2Ellipse : public ConstraintInfo, protected ConstraintVariables<InternalAlignmentPoint2Ellipse>
    {
    private:
        InternalAlignmentType AlignmentType;
    public:
        ConstraintInfoInternalAlignmentPoint2Ellipse(Ellipse &e, Point &p1, InternalAlignmentType alignmentType);
        ConstraintInfoInternalAlignmentPoint2Ellipse(ArcOfEllipse &e, Point &p1, InternalAlignmentType alignmentType);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    class ConstraintInfoEqualMajorAxesEllipse: public ConstraintInfo, protected ConstraintVariables<EqualMajorAxesEllipse>
    {
    public:
        ConstraintInfoEqualMajorAxesEllipse(Ellipse &e1, Ellipse &e2);
        ConstraintInfoEqualMajorAxesEllipse(ArcOfEllipse &a1, Ellipse &e2);
        ConstraintInfoEqualMajorAxesEllipse(ArcOfEllipse &a1, ArcOfEllipse &a2);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

    class ConstraintInfoEllipticalArcRangeToEndPoints: public ConstraintInfo, protected ConstraintVariables<EllipticalArcRangeToEndPoints>
    {
    public:
        ConstraintInfoEllipticalArcRangeToEndPoints(Point &p, ArcOfEllipse &a, double* const angle_t);
        virtual Constraint* createConstraint(SubSystem& sub_system) const;
        virtual ConstraintType getTypeId() const ;
        virtual ConstraintInfo* clone() const;
    };

} //namespace GCS_EXP

#endif // FREEGCS_CONSTRAINTSINFO_EXP_H
