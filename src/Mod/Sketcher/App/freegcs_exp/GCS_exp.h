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

#ifndef FREEGCS_GCS_EXP_H
#define FREEGCS_GCS_EXP_H

#include "SubSystem_exp.h"
#include "ConstraintsInfo_exp.h"
namespace GCS_EXP
{

    ///////////////////////////////////////
    // Solver
    ///////////////////////////////////////


    class System
    {
    // This is the main class. It holds all constraints and information
    // about partitioning into subsystems and solution strategies
    public:
    	typedef std::vector<ConstraintInfo *> constraits_container;
    private:
    	constraits_container constraints_list;
        constraits_container reduced_constraints_list;
        std::vector<index_type> constraints_to_subsystem;
        /* dependent_variables is a vector rather than a set because
         * in initSolution() the distance is taken between iterators
         * which is slower in set. O(1) in vector O(n)? in set.
         ***/
        std::vector< double*> dependent_variables;
//        std::map< double*, index_type > p2subsystem;

        std::vector<SubSystem *> subSystems;
        void clearSubSystems();

        SubSystem* getConstraintSubSystem( int id ){
        	if( isInit && id >= 0 && id < constraints_to_subsystem.size() )
        		return subSystems[ constraints_to_subsystem[id]];
        	return static_cast<SubSystem*>(0);
        }

        void setReference();     // copies the current parameter values to reference
        void resetToReference(); // reverts all parameter values to the stored reference

        int dofs;
        std::set<ConstraintInfo *> redundant;
        VEC_I conflictingTags, redundantTags;

        bool hasUnknowns;  // if plist is filled with the unknown parameters
        bool hasDiagnosis; // if dofs, conflictingTags, redundantTags are up to date
        bool isInit;       // if plists, clists, reductionmaps are up to date


    public:
        System();
        System(constraits_container clist_);
        ~System();

        void clear();
        void clearByTag(int tagId);

        int addReducedConstraint( ConstraintInfo *constr );

        int addConstraint(ConstraintInfo *constr);
        void removeConstraint(ConstraintInfo *constr);



        // basic constraints
        int addConstraintEqual( double* param1, double* param2, int tagId=0 );
        int addConstraintDifference( double* param1, double* param2,
        		double* difference, int tagId=0 );
        int addConstraintP2PDistance( Point& p1, Point& p2, double* distance, int tagId=0 );
        int addConstraintP2PAngle( Point& p1, Point& p2, double* angle,
                                  double incrAngle, int tagId=0 );
        int addConstraintP2PAngle( Point& p1, Point& p2, double* angle, int tagId=0 );
        int addConstraintP2LDistance( Point& p, Line& l, double* distance, int tagId=0 );
        int addConstraintPointOnLine( Point& p, Line& l, int tagId=0 );
        int addConstraintPointOnLine( Point& p, Point& lp1, Point& lp2, int tagId=0 );
        int addConstraintPointOnPerpBisector( Point& p, Line& l, int tagId=0 );
        int addConstraintPointOnPerpBisector( Point& p, Point& lp1, Point& lp2, int tagId=0 );
        int addConstraintParallel( Line& l1, Line& l2, int tagId=0 );
        int addConstraintPerpendicular( Line& l1, Line& l2, int tagId=0 );
        int addConstraintPerpendicular( Point& l1p1, Point& l1p2,
                                       Point& l2p1, Point& l2p2, int tagId=0 );
        int addConstraintL2LAngle( Line& l1, Line& l2, double* angle, int tagId=0 );
        int addConstraintL2LAngle( Point& l1p1, Point& l1p2, Point& l2p1, Point& l2p2,
        		double* angle, int tagId=0 );
        int addConstraintMidpointOnLine( Line& l1, Line& l2, int tagId=0);
        int addConstraintMidpointOnLine( Point& l1p1, Point& l1p2, Point& l2p1, Point& l2p2,
                                        int tagId=0 );
        int addConstraintTangentCircumf( Point& p1, Point& p2, double* rd1, double* rd2,
                                        bool internal=false, int tagId=0 );

        // derived constraints
        int addConstraintP2PCoincident( Point& p1, Point& p2, int tagId=0 );
        int addConstraintHorizontal( Line& l, int tagId=0 );
        int addConstraintHorizontal( Point& p1, Point& p2, int tagId=0 );
        int addConstraintVertical( Line& l, int tagId=0 );
        int addConstraintVertical( Point& p1, Point& p2, int tagId=0 );
        int addConstraintCoordinateX( Point& p, double* x, int tagId=0 );
        int addConstraintCoordinateY( Point& p, double* y, int tagId=0 );
        int addConstraintArcRules( Arc& a, int tagId=0 );
        int addConstraintPointOnCircle( Point& p, Circle& c, int tagId=0);
        int addConstraintPointOnArc( Point& p, Arc& a, int tagId=0 );
        int addConstraintPerpendicularLine2Arc( Point& p1, Point& p2, Arc& a,
                                               int tagId=0 );
        int addConstraintPerpendicularArc2Line( Arc& a, Point& p1, Point& p2,
                                               int tagId=0 );
        int addConstraintPerpendicularCircle2Arc( Point& center, double *radius, Arc& a,
                                                 int tagId=0 );
        int addConstraintPerpendicularArc2Circle( Arc& a, Point& center, double *radius,
                                                 int tagId=0 );
        int addConstraintPerpendicularArc2Arc( Arc& a1, bool reverse1,
                                              Arc& a2, bool reverse2, int tagId=0 );
        int addConstraintTangent( Line& l, Circle& c, int tagId=0 );
        int addConstraintTangent( Line& l, Arc& a, int tagId=0 );
        int addConstraintTangent( Circle& c1, Circle& c2, int tagId=0 );
        int addConstraintTangent( Arc& a1, Arc& a2, int tagId=0 );
        int addConstraintTangent( Circle& c, Arc& a, int tagId=0 );
        int addConstraintTangentLine2Arc( Point& p1, Point& p2, Arc& a, int tagId=0 );
        int addConstraintTangentArc2Line( Arc& a, Point& p1, Point& p2, int tagId=0 );
        int addConstraintTangentCircle2Arc( Circle& c, Arc& a, int tagId=0 );
        int addConstraintTangentArc2Circle( Arc& a, Circle& c, int tagId=0 );
        int addConstraintTangentArc2Arc( Arc& a1, bool reverse1, Arc& a2, bool reverse2,
                                        int tagId=0 );
        int addConstraintCircleRadius( Circle&c, double *radius, int tagId=0 );
        int addConstraintArcRadius( Arc&a, double* radius, int tagId=0 );
        int addConstraintEqualLength( Line& l1, Line& l2, double* length, int tagId=0 );
        int addConstraintEqualRadius( Circle& c1, Circle& c2, int tagId=0 );
        int addConstraintEqualRadius( Circle& c1, Arc& a2, int tagId=0 );
        int addConstraintEqualRadius( Arc& a1, Arc& a2, int tagId=0 );
        int addConstraintP2PSymmetric( Point& p1, Point& p2, Line& l, int tagId=0 );
        int addConstraintP2PSymmetric( Point& p1, Point& p2, Point& p, int tagId=0 );

        void rescaleConstraint( int id, double coeff );

        void declareUnknowns( std::vector<double *> &params );
        void initSolution();

        int solve(bool isFine=true, Algorithm alg=DogLeg);
        int solve(std::vector<double *> &params, bool isFine=true, Algorithm alg=DogLeg);

        void applySolution();
        void undoSolution();

        int diagnose();
        int dofsNumber() { return hasDiagnosis ? dofs : -1; }
        void getConflicting(VEC_I &conflictingOut) const
          { conflictingOut = hasDiagnosis ? conflictingTags : VEC_I(0); }
        void getRedundant(VEC_I &redundantOut) const
          { redundantOut = hasDiagnosis ? redundantTags : VEC_I(0); }
    };


    ///////////////////////////////////////
    // Helper elements
    ///////////////////////////////////////

    void free(std::vector<double *> &doublevec);
    void free(std::vector<ConstraintInfo *> &constrvec);
    void free(std::vector<SubSystem *> &subsysvec);

} //namespace GCS_EXP

#endif // FREEGCS_GCS_EXP_H
