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

#ifndef FREEGCS_SUBSYSTEM_EXP_H
#define FREEGCS_SUBSYSTEM_EXP_H

#undef min
#undef max

#include <Eigen/Core>
#include "Constraints_exp.h"
#include "ConstraintsInfo_exp.h"

namespace GCS_EXP
{

enum SolveStatus {
    Success = 0,   // Found a solution zeroing the error function
    Converged = 1, // Found a solution minimizing the error function
    Failed = 2     // Failed to find any solution
};

enum Algorithm {
    BFGS = 0,
    LevenbergMarquardt = 1,
    DogLeg = 2
};


    class SubSystem
    {
    private:
        int dependent_variable_count;
        int priority_constraints_count;
        // current variables vector
    	std::vector<double> variables;
    	// map of the external variables to local variables
        std::map< double*, variable_index_type > parameter_indices;
        std::vector<Constraint *> constraints_list;
//        // variable to constraints adjacency list
//        std::vector< std::vector<index_type> > p2c;
        // reference values of the variables vector
        std::vector<double> reference_values;

        // Diagnose system
        int dofs;
        bool hasDiagnosis;
//        bool hasUnknowns;
        std::set<Constraint*> redundant;
        std::vector<index_type> conflictingTags, redundantTags;

        void initialize();

        variable_index_type getIndex( double* param ){
        	typedef std::map< double*, variable_index_type > map_type;
        	map_type::iterator param_it =
        			parameter_indices.find( param );
        	if( param_it != parameter_indices.end())
        		return param_it->second;

        	parameter_indices.insert( map_type::value_type( param, variables.size() ));
        	variables.push_back( *param );
        	return parameter_indices[ param ];
        }

        bool isDependentVariable( variable_index_type i ){
        	return ( i < dependent_variable_count );
        }
        double& getDependentVariable( variable_index_type i ){
        	assert( isDependentVariable(i) );
        	return variables[ i ];
        }

        bool isPriorityConstraint( index_type index){
        	return index < priority_constraints_count;
        }
        bool isAuxiliaryConstraint( index_type index ){
        	return index >= priority_constraints_count;
        }

        void addDependentVariables( const std::vector< double* >& dependent_vars );
        void addConstraints( const std::vector<ConstraintInfo *> &clist_ );

        template <typename IteratorType>
        void calcGrad( IteratorType it, const IteratorType end, Eigen::VectorXd &grad );
        template <typename IteratorType>
        void calcJacobi( IteratorType it, const IteratorType it_end, Eigen::MatrixXd &jacobi );

        int solve_BFGS( bool isFine);
        int solve_LM();
        int solve_DL();

    public:
        SubSystem( const std::vector<ConstraintInfo *>& clist_, const std::vector<double*>& dependent_variables );
        SubSystem( const std::vector<Constraint *>& clist_ , const std::vector<variable_index_type>& dep_var );
        ~SubSystem();

        void setReference();
        void resetToReference();

        void rescaleConstraint(int id, double coeff)/*{
        if ( id >= constraints_list.size() || id < 0 )
            return;
        SubSystem* constraint_subsys = getSubsystem(id);
        if ( constraint_subsys )
        	constraint_subsys->rescaleConstraint(id,coeff);
        	}*/;

        const std::vector<double>& getVariables() const { return variables; }
        int getDependentVariableCount() const { return dependent_variable_count; }
        std::vector<variable_index_type> getIndices(
        		const std::vector<double*>& original_variables);

//        void getParamMap(MAP_pD_pD &pmapOut);
//        void getParamList(std::vector<double *> &plistOut);

        void getParams(Eigen::VectorXd &xOut);
        void setParams(Eigen::VectorXd &xIn);

        double error();
        void calcResidual(Eigen::VectorXd &r);
        void calcResidual(Eigen::VectorXd &r, double &err);
        void calcJacobi(Eigen::MatrixXd &jacobi);
        void calcGrad( Eigen::VectorXd& grad );
        double maxStep( Eigen::VectorXd &xdir);

        // For priority constraints of the system ( tag >= 0 )
        double errorPriority();
        void calcResidualPriority(Eigen::VectorXd &r);
        void calcJacobiPriority(Eigen::MatrixXd &jacobi);
        void calcGradPriority( Eigen::VectorXd& grad );
        double maxStepPriority( Eigen::VectorXd &xdir);

        // For auxiliary constraints of the system ( tag < 0 )
        double errorAuxiliary();
        void calcResidualAuxiliary(Eigen::VectorXd &r);
        void calcJacobiAuxiliary(Eigen::MatrixXd &jacobi);
        void calcGradAuxiliary( Eigen::VectorXd& grad );
        double maxStepAuxiliary( Eigen::VectorXd &xdir);

        void applySolution();
        void updateSystemParameters();

        void analyse(Eigen::MatrixXd &J, Eigen::MatrixXd &ker, Eigen::MatrixXd &img);
        void report();
        int diagnose();

        int solve( bool isFine = true, Algorithm alg = DogLeg );

        void printResidual();
    };

    ///////////////////////////////////////
    // BFGS Solver parameters
    ///////////////////////////////////////
	#define XconvergenceRough 1e-8
	#define XconvergenceFine  1e-10
	#define smallF            1e-20
	#define MaxIterations     100 //Note that the total number of iterations allowed is MaxIterations *xLength
//    static const double XconvergenceRough = 1e-8;
//    static const double XconvergenceFine  = 1e-10;
//    static const double smallF            = 1e-20;
//    static const size_t MaxIterations     = 100; //Note that the total number of iterations allowed is MaxIterations *xLength

    double lineSearch(SubSystem *subsys, Eigen::VectorXd &xdir);

} //namespace GCS_EXP

#endif // FREEGCS_SUBSYSTEM_EXP_H
