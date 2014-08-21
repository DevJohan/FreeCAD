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

#include "SubSystem_exp.h"

#include <iostream>
#include <iterator>
#include <cfloat>

#include "qp_eq_exp.h"
#include <Eigen/QR>


namespace GCS_EXP
{

// SubSystem
SubSystem::SubSystem(
		const std::vector<ConstraintInfo *>& clist_,
		const std::vector<double*>& dependent_vars
):
		dependent_variable_count(0),
		priority_constraints_count(0),
		variables(),
		parameter_indices(),
//		dependent_variables(),
		constraints_list(),
		dofs(-1),
		hasDiagnosis(false)
//		hasUnknowns(false)
{
	/***
	 * Dependent variables are given indices first so that all dependent variables
	 * are in the beginning of the variables vector. Which variables are dependent
	 * is determined by checking index smaller than dependent_variable_count.
	 */
	addDependentVariables( dependent_vars );

	addConstraints( clist_ );
	initialize();
}

SubSystem::~SubSystem()
{
}

void SubSystem::addDependentVariables( const std::vector< double* >& dependent_vars ){
	typedef std::vector<double*>::const_iterator dpv_iterator;
    for( dpv_iterator dp_iter = dependent_vars.begin();
    		dp_iter != dependent_vars.end(); ++dp_iter ) {
    	getIndex( *dp_iter );
    }
    dependent_variable_count = variables.size();
}

void SubSystem::addConstraints( const std::vector<ConstraintInfo *>& c_info_list ){
	typedef std::vector<Constraint*>::iterator front_iter;
	typedef std::vector<Constraint*>::reverse_iterator back_iter;

	constraints_list.resize( c_info_list.size() );
	front_iter prio_constraints = constraints_list.begin();
	back_iter aux_constraints = constraints_list.rbegin();
    for (std::vector<ConstraintInfo *>::const_iterator c_info_it = c_info_list.begin();
    		c_info_it != c_info_list.end(); ++c_info_it ) {
    	if( (*c_info_it)->isPriorityConstraint() ){
    		*prio_constraints = (*c_info_it)->createConstraint(*this);
    		prio_constraints++;
    	}else{
    		*aux_constraints = (*c_info_it)->createConstraint(*this);
    		aux_constraints++;
    	}
    }
    priority_constraints_count =
    		std::distance( constraints_list.begin(), prio_constraints );

    assert( std::distance( constraints_list.rbegin(), aux_constraints ) +
    		priority_constraints_count == constraints_list.size() );
}

void SubSystem::initialize()
{
//	p2c.resize( dependent_variable_count );
//    for (std::vector<Constraint *>::iterator constr = constraints_list.begin();
//    		constr != constraints_list.end(); ++constr) {
//    	std::vector<variable_index_type> constr_params = (*constr)->indices();
//    	for(std::vector<variable_index_type>::iterator param_it = constr_params.begin();
//    			param_it != constr_params.end(); ++param_it ){
//    		const variable_index_type index = *param_it;
//    		if( isDependentVariable( index ) )
//    			p2c[ index ].push_back( std::distance(constraints_list.begin(), constr) );
//    	}
//    }
}


std::vector<variable_index_type> SubSystem::getIndices(
		const std::vector<double*>& original_variables
){
	std::vector<variable_index_type> variable_indices;
	for(std::vector<double*>::const_iterator var_it = original_variables.begin();
			var_it != original_variables.end(); ++var_it ){
		variable_indices.push_back( getIndex(*var_it) );
	}
	assert( variable_indices.size() == original_variables.size() );
	return variable_indices;
}

void SubSystem::setReference(){
	reference_values.assign( variables.begin(), variables.end() );
}

void SubSystem::resetToReference(){
	assert( reference_values.size() == variables.size() );
	variables.assign(reference_values.begin(), reference_values.end());
}

void SubSystem::getParams( Eigen::VectorXd &xOut )
{
	const int psize = dependent_variable_count;
    if ( xOut.size() != psize )
        xOut.setZero( psize );

    for (int i=0; i < psize; i++)
        xOut[i] = getDependentVariable(i);
}

void SubSystem::setParams( Eigen::VectorXd &xIn )
{
	const int psize = dependent_variable_count;
    assert(xIn.size() == psize);
    for (int i=0; i < psize; i++)
        getDependentVariable(i) = xIn[i];
}

double SubSystem::error()
{
    double err = 0.;
    for (std::vector<Constraint *>::const_iterator constr=constraints_list.begin();
         constr != constraints_list.end(); ++constr) {
        double tmp = (*constr)->error();
        err += tmp*tmp;
    }
    err *= 0.5;
    return err;
}

double SubSystem::errorPriority()
{
    double err = 0.;
    const int priority_constraints = priority_constraints_count;
    for( int i = 0; i < priority_constraints; ++i){
        double tmp = constraints_list[i]->error();
        err += tmp*tmp;
    }
    err *= 0.5;
    return err;
}

// auxiliary
double SubSystem::errorAuxiliary()
{
    double err = 0.;
    const int end_index = constraints_list.size();
    for( int i = priority_constraints_count; i < end_index; ++i ){
        double tmp = constraints_list[i]->error();
        err += tmp*tmp;
    }
    err *= 0.5;
    return err;
}

void SubSystem::calcResidual(Eigen::VectorXd &r)
{
    assert(r.size() == constraints_list.size());

    int i=0;
    for (std::vector<Constraint *>::const_iterator constr=constraints_list.begin();
         constr != constraints_list.end(); ++constr, i++) {
        r[i] = (*constr)->error();
    }
}

void SubSystem::calcResidual(Eigen::VectorXd &r, double &err)
{
    assert(r.size() == constraints_list.size());

    int i=0;
    err = 0.;
    for (std::vector<Constraint *>::const_iterator constr=constraints_list.begin();
         constr != constraints_list.end(); ++constr, i++) {
        r[i] = (*constr)->error();
        err += r[i]*r[i];
    }
    err *= 0.5;
}

void SubSystem::calcResidualPriority( Eigen::VectorXd &r )
{
	const int prio_count = priority_constraints_count ;
    assert( r.size() == prio_count );

    for( int i=0; i < prio_count; ++i ){
        r[i] = constraints_list[i]->error();
    }
}

void SubSystem::calcResidualAuxiliary( Eigen::VectorXd &r )
{
	const int end_constraint = constraints_list.size();
	const int aux_count = end_constraint - priority_constraints_count;
    assert( r.size() == aux_count );

    for( int i = priority_constraints_count; i < end_constraint; ++i ){
        r[i] = constraints_list[i]->error();
    }
}

template <typename IteratorType>
void SubSystem::calcJacobi( IteratorType it, const IteratorType it_end, Eigen::MatrixXd &jacobi )
{
    typedef Constraint::grad_component_t grad_component_t;
    typedef std::vector< Constraint* >::const_iterator constrait_iterator;
    typedef std::vector< grad_component_t >::iterator grad_iterator;

    std::vector< grad_component_t > gradComp;
    int cid = 0;
    for ( ; it != it_end; ++it, ++cid ){
    	gradComp.clear();
    	(*it)->grad( gradComp );
    	for( grad_iterator comp_it = gradComp.begin();
    			comp_it != gradComp.end(); ++comp_it )
    	{
    		jacobi( cid, comp_it->first ) += comp_it->second;
    	}
    }
}

void SubSystem::calcJacobi( Eigen::MatrixXd &jacobi )
{
    jacobi.setZero( constraints_list.size(), dependent_variable_count );
    calcJacobi( constraints_list.begin(), constraints_list.end(), jacobi );
}

void SubSystem::calcJacobiPriority( Eigen::MatrixXd &jacobi )
{
    jacobi.setZero( priority_constraints_count, dependent_variable_count );
	std::vector<Constraint*>::iterator it_end = constraints_list.begin();
	std::advance( it_end, priority_constraints_count );
	calcJacobi( constraints_list.begin(), it_end, jacobi );
}

void SubSystem::calcJacobiAuxiliary( Eigen::MatrixXd &jacobi )
{
    jacobi.setZero( constraints_list.size()-priority_constraints_count, dependent_variable_count );
	std::vector<Constraint*>::iterator it = constraints_list.begin();
	std::advance( it, priority_constraints_count );
	calcJacobi( it, constraints_list.end(), jacobi );
}

template <typename IteratorType>
void SubSystem::calcGrad( IteratorType it, const IteratorType it_end, Eigen::VectorXd &grad )
{
    typedef Constraint::grad_component_t grad_component_t;
    typedef std::vector< Constraint* >::const_iterator constrait_iterator;
    typedef std::vector< grad_component_t >::iterator grad_iterator;

    assert( grad.size() == dependent_variable_count );

    grad.setZero();
    std::vector< grad_component_t > gradComp;
    for ( ; it != it_end; ++it ){
    	const double error = (*it)->error();
    	gradComp.clear();
    	(*it)->grad( gradComp );
    	for( grad_iterator comp_it = gradComp.begin();
    			comp_it != gradComp.end(); ++comp_it )
    		grad[ comp_it->first ] += error * comp_it->second;
    }
}

void SubSystem::calcGrad( Eigen::VectorXd &grad )
{
    calcGrad( constraints_list.begin(), constraints_list.end(), grad );
}

void SubSystem::calcGradPriority( Eigen::VectorXd &grad )
{
	std::vector<Constraint*>::iterator it_end = constraints_list.begin();
	std::advance( it_end, priority_constraints_count );
    calcGrad( constraints_list.begin(), it_end, grad );
}

void SubSystem::calcGradAuxiliary( Eigen::VectorXd &grad )
{
	std::vector<Constraint*>::iterator it = constraints_list.begin();
	std::advance( it, priority_constraints_count );
    calcGrad( it, constraints_list.end(), grad );
}

double SubSystem::maxStep( Eigen::VectorXd &xdir )
{
	const int dv_count = dependent_variable_count;
	assert( xdir.size() == dv_count );

	std::vector<double> dir( dv_count );
	for (int j=0; j < dv_count; j++) {
		dir[j] = xdir[j];
	}

	double alpha=1e10;
	for( std::vector<Constraint *>::iterator constr = constraints_list.begin();
			constr != constraints_list.end(); ++constr)
		alpha = (*constr)->maxStep( dir, alpha );

	return alpha;
}

double SubSystem::maxStepPriority( Eigen::VectorXd &xdir )
{
	typedef std::vector<Constraint *>::iterator cop_iterator;
	const int dv_count = dependent_variable_count;
	assert( xdir.size() == dv_count );

	std::vector<double> dir( dv_count );
	for (int j=0; j < dv_count; j++) {
		dir[j] = xdir[j];
	}

	double alpha=1e10;
	cop_iterator constr_end = constraints_list.begin()+ priority_constraints_count;
	for( cop_iterator constr = constraints_list.begin();
			constr != constr_end; ++constr)
	{
		alpha = (*constr)->maxStep( dir, alpha );
	}

	return alpha;
}

void SubSystem::applySolution()
{
	typedef std::map<double*, variable_index_type>::const_iterator cpi_iter;
	for (cpi_iter it = parameter_indices.begin();
			it != parameter_indices.end(); ++it ){
		if( isDependentVariable( it->second ) )
			*(it->first) = variables[ it->second ];
	}
}

void SubSystem::updateSystemParameters()
{
	typedef std::map<double*, variable_index_type>::const_iterator cpi_iter;
	for (cpi_iter it = parameter_indices.begin();
			it != parameter_indices.end(); ++it ){
		if( !isDependentVariable( it->second ) )
			variables[ it->second ] = *(it->first);
	}
}

void SubSystem::analyse( Eigen::MatrixXd &J, Eigen::MatrixXd &ker, Eigen::MatrixXd &img )
{
}

void SubSystem::report()
{
}

int SubSystem::diagnose()
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
//    if(!hasUnknowns) {
//        dofs = -1;
//        return dofs;
//    }

    redundant.clear();
    conflictingTags.clear();
    redundantTags.clear();
    Eigen::MatrixXd jacobian;
    calcJacobiPriority( jacobian );
//    int count=0;
//    for ( std::vector<Constraint*>::iterator constr=constraints_list.begin();
//         constr != constraints_list.end(); ++constr) {
//        if( (*constr)->getTag() >= 0) {
//            count++;
//            for (int j=0; j < int(variables.size()); j++)
//                jacobian(count-1,j) = (*constr)->grad(j);
//        }
//    }

    if( jacobian.rows() > 0 ) {
        Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qrJT( jacobian.transpose() );
        Eigen::MatrixXd Q = qrJT.matrixQ ();
        int paramsNum = qrJT.rows();
        int constrNum = qrJT.cols();
        int rank = qrJT.rank();

        Eigen::MatrixXd R;
        if( constrNum >= paramsNum )
            R = qrJT.matrixQR().triangularView<Eigen::Upper>();
        else
            R = qrJT.matrixQR().topRows(constrNum)
                               .triangularView<Eigen::Upper>();

        if( constrNum > rank ) { // conflicting or redundant constraints
            for( int i=1; i < rank; i++) {
                // eliminate non zeros above pivot
                assert(R(i,i) != 0);
                for( int row=0; row < i; row++) {
                    if (R(row,i) != 0) {
                        double coef = R(row,i) / R(i,i);
                        R.block( row, i+1, 1, constrNum-i-1 ) -= coef * R.block(i,i+1,1,constrNum-i-1);
                        R(row,i) = 0;
                    }
                }
            }
            std::vector< std::vector< Constraint* > > conflictGroups( constrNum-rank );
            for( int j=rank; j < constrNum; j++ ){
                for( int row=0; row < rank; row++ ){
                    if( fabs( R(row,j) ) > 1e-10 ){
                        int origCol = qrJT.colsPermutation().indices()[row];
                        conflictGroups[j-rank].push_back(constraints_list[origCol]);
                    }
                }
                int origCol = qrJT.colsPermutation().indices()[j];
                conflictGroups[j-rank].push_back(constraints_list[origCol]);
            }

            // try to remove the conflicting constraints and solve the
            // system in order to check if the removed constraints were
            // just redundant but not really conflicting
            std::set< Constraint* > skipped;
            SET_I satisfiedGroups;
            while (1) {
                std::map< Constraint*, SET_I > conflictingMap;
                for (int i=0; i < conflictGroups.size(); i++) {
                    if (satisfiedGroups.count(i) == 0) {
                        for (int j=0; j < conflictGroups[i].size(); j++) {
                            Constraint *constr = conflictGroups[i][j];
                            if (constr->getTag() != 0) // exclude constraints tagged with zero
                                conflictingMap[constr].insert(i);
                        }
                    }
                }
                if (conflictingMap.empty())
                    break;

                int maxPopularity = 0;
                Constraint *mostPopular = NULL;
                for (std::map< Constraint *, SET_I >::const_iterator it=conflictingMap.begin();
                     it != conflictingMap.end(); it++) {
                    if (it->second.size() > maxPopularity ||
                        (it->second.size() == maxPopularity && mostPopular &&
                         it->first->getTag() > mostPopular->getTag())) {
                        mostPopular = it->first;
                        maxPopularity = it->second.size();
                    }
                }
                if (maxPopularity > 0) {
                    skipped.insert(mostPopular);
                    for (SET_I::const_iterator it=conflictingMap[mostPopular].begin();
                         it != conflictingMap[mostPopular].end(); it++)
                        satisfiedGroups.insert(*it);
                }
            }

//            std::vector< Constraint* > clistTmp;
//            clistTmp.reserve( constraints_list.size() );
//            for ( std::vector< Constraint* >::iterator constr=constraints_list.begin();
//            		constr != constraints_list.end(); ++constr)
//            	if( skipped.count( *constr ) == 0 )
//            		clistTmp.push_back(*constr);
//
//            std::vector<variable_index_type> dependent_variables( dependent_variable_count );
//            for( int i=0; i < dependent_variable_count; i++)
//            	dependent_variables[i] = i;
//            SubSystem* subSysTmp = new SubSystem( clistTmp , dependent_variables );
//            int res = subSysTmp->solve();
//            if ( res == Success ) {
//            	subSysTmp->applySolution();
//            	for( std::set<Constraint *>::const_iterator constr=skipped.begin();
//            			constr != skipped.end(); ++constr ) {
//            		double err = (*constr)->error();
//            		if( err * err < XconvergenceFine )
//            			redundant.insert( *constr );
//            	}
//            	resetToReference();
//
//            	std::vector< std::vector<Constraint *> > conflictGroupsOrig=conflictGroups;
//            	conflictGroups.clear();
//            	for( int i = conflictGroupsOrig.size()-1; i >= 0; i-- ){
//            		bool isRedundant = false;
//            		for( int j = 0; j < conflictGroupsOrig[i].size(); j++ ){
//            			if( redundant.count( conflictGroupsOrig[i][j] ) > 0 ){
//            				isRedundant = true;
//            				break;
//            			}
//            		}
//                    if( !isRedundant )
//                        conflictGroups.push_back( conflictGroupsOrig[i] );
//                    else
//                        constrNum--;
//                }
//            }
//            delete subSysTmp;

            // simplified output of conflicting tags
            std::set<index_type> conflictingTagsSet;
            for (int i=0; i < conflictGroups.size(); i++) {
                for (int j=0; j < conflictGroups[i].size(); j++) {
                    conflictingTagsSet.insert(conflictGroups[i][j]->getTag());
                }
            }
            conflictingTagsSet.erase(0); // exclude constraints tagged with zero
            conflictingTags.resize(conflictingTagsSet.size());
            std::copy(conflictingTagsSet.begin(), conflictingTagsSet.end(),
                      conflictingTags.begin());

            // output of redundant tags
            std::set<index_type> redundantTagsSet;
            for( std::set< Constraint* >::iterator constr = redundant.begin();
                 constr != redundant.end(); ++constr)
                redundantTagsSet.insert((*constr)->getTag());

            // remove tags represented at least in one non-redundant constraint
            for( std::vector< Constraint* >::iterator constr = constraints_list.begin();
                 constr != constraints_list.end(); ++constr)
                if (redundant.count(*constr) == 0)
                    redundantTagsSet.erase((*constr)->getTag());

            redundantTags.resize( redundantTagsSet.size() );
            std::copy( redundantTagsSet.begin(), redundantTagsSet.end(),
                      redundantTags.begin());

            if ( paramsNum == rank && constrNum > rank ) { // over-constrained
                hasDiagnosis = true;
                dofs = paramsNum - constrNum;
                return dofs;
            }
        }

        hasDiagnosis = true;
        dofs = paramsNum - rank;
        return dofs;
    }
    hasDiagnosis = true;
    dofs = variables.size();
    return dofs;
}


int SubSystem::solve_BFGS( bool isFine)
{
    int xsize = dependent_variable_count;
    if (xsize == 0)
        return Success;

    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(xsize, xsize);
    Eigen::VectorXd x(xsize);
    Eigen::VectorXd xdir(xsize);
    Eigen::VectorXd grad(xsize);
    Eigen::VectorXd h(xsize);
    Eigen::VectorXd y(xsize);
    Eigen::VectorXd Dy(xsize);

    // Initial unknowns vector and initial gradient vector
    getParams(x);
    calcGrad(grad);

    // Initial search direction opposed to gradient (steepest-descent)
    xdir = -grad;
    lineSearch( this, xdir );
    double err = error();

    h = x;
    getParams(x);
    h = x - h; // = x - xold

    double convergence = isFine ? XconvergenceFine : XconvergenceRough;
    int maxIterNumber = MaxIterations * xsize;
    double divergingLim = 1e6*err + 1e12;

    for (int iter=1; iter < maxIterNumber; iter++) {

        if (h.norm() <= convergence || err <= smallF)
            break;
        if (err > divergingLim || err != err) // check for diverging and NaN
            break;

        y = grad;
        calcGrad( grad );
        y = grad - y; // = grad - gradold

        double hty = h.dot(y);
        //make sure that hty is never 0
        if (hty == 0)
            hty = .0000000001;

        Dy = D * y;

        double ytDy = y.dot(Dy);

        //Now calculate the BFGS update on D
        D += (1.+ytDy/hty)/hty * h * h.transpose();
        D -= 1./hty * (h * Dy.transpose() + Dy * h.transpose());

        xdir = -D * grad;
        lineSearch( this, xdir);
        err = error();

        h = x;
        getParams(x);
        h = x - h; // = x - xold
    }

    if (err <= smallF)
        return Success;
    if (h.norm() <= convergence)
        return Converged;
    return Failed;
}

int SubSystem::solve_LM()
{
    const int xsize = dependent_variable_count;
    const int csize = constraints_list.size();

    if (xsize == 0)
        return Success;

    Eigen::VectorXd e(csize), e_new(csize); // vector of all function errors (every constraint is one function)
    Eigen::MatrixXd J(csize, xsize);        // Jacobi of the subsystem
    Eigen::MatrixXd A(xsize, xsize);
    Eigen::VectorXd x(xsize), h(xsize), x_new(xsize), g(xsize), diag_A(xsize);

    getParams(x);
    calcResidual(e);
    e*=-1;

    int maxIterNumber = MaxIterations * xsize;
    double divergingLim = 1e6*e.squaredNorm() + 1e12;

    double eps=1e-10, eps1=1e-80;
    double tau=1e-3;
    double nu=2, mu=0;
    int iter=0, stop=0;
    for (iter=0; iter < maxIterNumber && !stop; ++iter) {

        // check error
        double err = e.squaredNorm();
        if( err <= eps ) { // error is small, Success
            stop = 1;
            break;
        }else if( err > divergingLim || err != err ){ // check for diverging and NaN
            stop = 6;
            break;
        }

        // J^T J, J^T e
        calcJacobi(J);;

        A = J.transpose()*J;
        g = J.transpose()*e;

        // Compute ||J^T e||_inf
        double g_inf = g.lpNorm<Eigen::Infinity>();
        diag_A = A.diagonal(); // save diagonal entries so that augmentation can be later canceled

        // check for convergence
        if( g_inf <= eps1 ){
            stop = 2;
            break;
        }

        // compute initial damping factor
        if (iter == 0)
            mu = tau * diag_A.lpNorm<Eigen::Infinity>();

        // determine increment using adaptive damping
        int k=0;
        while (k < 50) {
            // augment normal equations A = A+uI
            for (int i=0; i < xsize; ++i)
                A(i,i) += mu;

            //solve augmented functions A*h=-g
            h = A.fullPivLu().solve(g);
            double rel_error = (A*h - g).norm() / g.norm();

            // check if solving works
            if (rel_error < 1e-5) {

                // restrict h according to maxStep
                double scale = maxStep( h );
                if (scale < 1.)
                    h *= scale;

                // compute par's new estimate and ||d_par||^2
                x_new = x + h;
                double h_norm = h.squaredNorm();

                if (h_norm <= eps1*eps1*x.norm()) { // relative change in p is small, stop
                    stop = 3;
                    break;
                }
                else if (h_norm >= (x.norm()+eps1)/(DBL_EPSILON*DBL_EPSILON)) { // almost singular
                    stop = 4;
                    break;
                }

                setParams(x_new);
                calcResidual(e_new);
                e_new *= -1;

                double dF = e.squaredNorm() - e_new.squaredNorm();
                double dL = h.dot(mu*h+g);

                if (dF>0. && dL>0.) { // reduction in error, increment is accepted
                    double tmp=2*dF/dL-1.;
                    mu *= std::max(1./3., 1.-tmp*tmp*tmp);
                    nu=2;

                    // update par's estimate
                    x = x_new;
                    e = e_new;
                    break;
                }
            }

            // if this point is reached, either the linear system could not be solved or
            // the error did not reduce; in any case, the increment must be rejected

            mu*=nu;
            nu*=2.0;
            for (int i=0; i < xsize; ++i) // restore diagonal J^T J entries
                A(i,i) = diag_A(i);

            k++;
        }
        if (k > 50) {
            stop = 7;
            break;
        }
    }

    if( iter >= maxIterNumber )
        stop = 5;

    return (stop == 1) ? Success : Failed;
}


int SubSystem::solve_DL()
{
    double tolg=1e-80, tolx=1e-80, tolf=1e-10;

    const int xsize = dependent_variable_count;
    const int csize = constraints_list.size();

    if (xsize == 0)
        return Success;

    Eigen::VectorXd x(xsize), x_new(xsize);
    Eigen::VectorXd fx(csize), fx_new(csize);
    Eigen::MatrixXd Jx(csize, xsize), Jx_new(csize, xsize);
    Eigen::VectorXd g(xsize), h_sd(xsize), h_gn(xsize), h_dl(xsize);

    double err;
    getParams(x);
    calcResidual(fx, err);
    calcJacobi(Jx);

    g = Jx.transpose()*(-fx);

    // get the infinity norm fx_inf and g_inf
    double g_inf = g.lpNorm<Eigen::Infinity>();
    double fx_inf = fx.lpNorm<Eigen::Infinity>();

    int maxIterNumber = MaxIterations * xsize;
    double divergingLim = 1e6*err + 1e12;

    double delta=0.1;
    double alpha=0.;
    double nu=2.;
    int iter=0, stop=0, reduce=0;
    while (!stop) {

        // check if finished
        if (fx_inf <= tolf) // Success
            stop = 1;
        else if (g_inf <= tolg)
            stop = 2;
        else if (delta <= tolx*(tolx + x.norm()))
            stop = 2;
        else if (iter >= maxIterNumber)
            stop = 4;
        else if (err > divergingLim || err != err) { // check for diverging and NaN
            stop = 6;
        }
        else {
            // get the steepest descent direction
            alpha = g.squaredNorm()/(Jx*g).squaredNorm();
            h_sd  = alpha*g;

            // get the gauss-newton step
            h_gn = Jx.fullPivLu().solve(-fx);
            double rel_error = (Jx*h_gn + fx).norm() / fx.norm();
            if (rel_error > 1e15)
                break;

            // compute the dogleg step
            if (h_gn.norm() < delta) {
                h_dl = h_gn;
                if  (h_dl.norm() <= tolx*(tolx + x.norm())) {
                    stop = 5;
                    break;
                }
            }
            else if (alpha*g.norm() >= delta) {
                h_dl = (delta/(alpha*g.norm()))*h_sd;
            }
            else {
                //compute beta
                double beta = 0;
                Eigen::VectorXd b = h_gn - h_sd;
                double bb = (b.transpose()*b).norm();
                double gb = (h_sd.transpose()*b).norm();
                double c = (delta + h_sd.norm())*(delta - h_sd.norm());

                if (gb > 0)
                    beta = c / (gb + sqrt(gb * gb + c * bb));
                else
                    beta = (sqrt(gb * gb + c * bb) - gb)/bb;

                // and update h_dl and dL with beta
                h_dl = h_sd + beta*b;
            }
        }

        // see if we are already finished
        if (stop)
            break;

// it didn't work in some tests
//        // restrict h_dl according to maxStep
//        double scale = subsys->maxStep(h_dl);
//        if (scale < 1.)
//            h_dl *= scale;

        // get the new values
        double err_new;
        x_new = x + h_dl;
        setParams(x_new);
        calcResidual(fx_new, err_new);
        calcJacobi(Jx_new);

        // calculate the linear model and the update ratio
        double dL = err - 0.5*(fx + Jx*h_dl).squaredNorm();
        double dF = err - err_new;
        double rho = dL/dF;

        if (dF > 0 && dL > 0) {
            x  = x_new;
            Jx = Jx_new;
            fx = fx_new;
            err = err_new;

            g = Jx.transpose()*(-fx);

            // get infinity norms
            g_inf = g.lpNorm<Eigen::Infinity>();
            fx_inf = fx.lpNorm<Eigen::Infinity>();
        }
        else
            rho = -1;

        // update delta
        if (fabs(rho-1.) < 0.2 && h_dl.norm() > delta/3. && reduce <= 0) {
            delta = 3*delta;
            nu = 2;
            reduce = 0;
        }
        else if (rho < 0.25) {
            delta = delta/nu;
            nu = 2*nu;
            reduce = 2;
        }
        else
            reduce--;

        // count this iteration and start again
        iter++;
    }

    return (stop == 1) ? Success : Failed;
}



int SubSystem::solve( bool isFine, Algorithm alg )
{
	updateSystemParameters();

	if( 	priority_constraints_count == 0 ||
			priority_constraints_count == constraints_list.size() )
	{
		if (alg == BFGS)
			return solve_BFGS(isFine);
		else if (alg == LevenbergMarquardt)
			return solve_LM();
		else // if (alg == DogLeg)
			return solve_DL();
	}

	// The following solver variant solves a system with a priority and an
	// auxiliary part.
	// Constraints with index less-than priority_constraints_count make up
	// the priority part. The rest make up the auxiliary part

    int csizeA = priority_constraints_count;
    int xsize = dependent_variable_count;

    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(xsize, xsize);
    Eigen::MatrixXd JA(csizeA, xsize);
    Eigen::MatrixXd Y,Z;

    Eigen::VectorXd resA(csizeA);
    Eigen::VectorXd lambda(csizeA), lambda0(csizeA), lambdadir(csizeA);
    Eigen::VectorXd x(xsize), x0(xsize), xdir(xsize), xdir1(xsize);
    Eigen::VectorXd grad(xsize);
    Eigen::VectorXd h(xsize);
    Eigen::VectorXd y(xsize);
    Eigen::VectorXd Bh(xsize);

    // We assume that there are no common constraints in Priority and Auxiliary

    getParams( x );

    calcGradAuxiliary( grad );
    calcJacobiPriority( JA );
    calcResidualPriority( resA );

    double convergence = isFine ? XconvergenceFine : XconvergenceRough;
    int maxIterNumber = MaxIterations * xsize;
    double divergingLim = 1e6*errorPriority() + 1e12;

    double mu = 0;
    lambda.setZero();
    for (int iter=1; iter < maxIterNumber; iter++) {
        int status = qp_eq( B, grad, JA, resA, xdir, Y, Z);
        if (status)
            break;

        x0 = x;
        lambda0 = lambda;
        lambda = Y.transpose() * (B * xdir + grad);
        lambdadir = lambda - lambda0;

        // line search
        {
            double eta=0.25;
            double tau=0.5;
            double rho=0.5;
            double alpha=1;
            alpha = std::min( alpha, maxStepPriority( xdir ) );

            // Eq. 18.32
            // double mu = lambda.lpNorm<Eigen::Infinity>() + 0.01;
            // Eq. 18.33
            // double mu =  grad.dot(xdir) / ( (1.-rho) * resA.lpNorm<1>());
            // Eq. 18.36
            mu =  std::max(mu,
                           (grad.dot(xdir) +  std::max(0., 0.5*xdir.dot(B*xdir))) /
                           ( (1. - rho) * resA.lpNorm<1>() ) );

            // Eq. 18.27
            double f0 = errorAuxiliary() + mu * resA.lpNorm<1>();

            // Eq. 18.29
            double deriv = grad.dot(xdir) - mu * resA.lpNorm<1>();

            x = x0 + alpha * xdir;
            setParams(x);
            calcResidualPriority(resA);
            double f = errorAuxiliary() + mu * resA.lpNorm<1>();

            // line search, Eq. 18.28
            bool first = true;
            while (f > f0 + eta * alpha * deriv) {
                if (first) { // try a second order step
//                    xdir1 = JA.jacobiSvd(Eigen::ComputeThinU |
//                                         Eigen::ComputeThinV).solve(-resA);
                    xdir1 = -Y*resA;
                    x += xdir1; // = x0 + alpha * xdir + xdir1
                    setParams(x);
                    calcResidualPriority(resA);
                    f = errorAuxiliary() + mu * resA.lpNorm<1>();
                    if (f < f0 + eta * alpha * deriv)
                        break;
                }
                alpha = tau * alpha;
                if (alpha < 1e-8) // let the linesearch fail
                    alpha = 0.;
                x = x0 + alpha * xdir;
                setParams(x);
                calcResidualPriority(resA);
                f = errorAuxiliary() + mu * resA.lpNorm<1>();
                if (alpha < 1e-8) // let the linesearch fail
                    break;
            }
            lambda = lambda0 + alpha * lambdadir;

        }
        h = x - x0;

        y = grad - JA.transpose() * lambda;
        {
            calcGradAuxiliary(grad);
            calcJacobiPriority(JA);
            calcResidualPriority(resA);
        }
        y = grad - JA.transpose() * lambda - y; // Eq. 18.13

        if (iter > 1) {
            double yTh = y.dot(h);
            if (yTh != 0) {
                Bh = B * h;
                //Now calculate the BFGS update on B
                B += 1./yTh * y * y.transpose();
                B -= 1./h.dot(Bh) * (Bh * Bh.transpose());
            }
        }

        double err = errorPriority();
        if (h.norm() <= convergence && err <= smallF)
            break;
        if (err > divergingLim || err != err) // check for diverging and NaN
            break;
    }

    int ret;
    if ( errorPriority() <= smallF )
        ret = Success;
    else if ( h.norm() <= convergence )
        ret = Converged;
    else
        ret = Failed;

//    revertParams();
    return ret;

}

double lineSearch(SubSystem *subsys, Eigen::VectorXd &xdir)
{
    double f1,f2,f3,alpha1,alpha2,alpha3,alphaStar;

    double alphaMax = subsys->maxStep(xdir);

    Eigen::VectorXd x0, x;

    //Save initial values
    subsys->getParams(x0);

    //Start at the initial position alpha1 = 0
    alpha1 = 0.;
    f1 = subsys->error();

    //Take a step of alpha2 = 1
    alpha2 = 1.;
    x = x0 + alpha2 * xdir;
    subsys->setParams(x);
    f2 = subsys->error();

    //Take a step of alpha3 = 2*alpha2
    alpha3 = alpha2*2;
    x = x0 + alpha3 * xdir;
    subsys->setParams(x);
    f3 = subsys->error();

    //Now reduce or lengthen alpha2 and alpha3 until the minimum is
    //Bracketed by the triplet f1>f2<f3
    while (f2 > f1 || f2 > f3) {
        if (f2 > f1) {
            //If f2 is greater than f1 then we shorten alpha2 and alpha3 closer to f1
            //Effectively both are shortened by a factor of two.
            alpha3 = alpha2;
            f3 = f2;
            alpha2 = alpha2 / 2;
            x = x0 + alpha2 * xdir;
            subsys->setParams(x);
            f2 = subsys->error();
        }
        else if (f2 > f3) {
            if (alpha3 >= alphaMax)
                break;
            //If f2 is greater than f3 then we increase alpha2 and alpha3 away from f1
            //Effectively both are lengthened by a factor of two.
            alpha2 = alpha3;
            f2 = f3;
            alpha3 = alpha3 * 2;
            x = x0 + alpha3 * xdir;
            subsys->setParams(x);
            f3 = subsys->error();
        }
    }
    //Get the alpha for the minimum f of the quadratic approximation
    alphaStar = alpha2 + ((alpha2-alpha1)*(f1-f3))/(3*(f1-2*f2+f3));

    //Guarantee that the new alphaStar is within the bracket
    if (alphaStar >= alpha3 || alphaStar <= alpha1)
        alphaStar = alpha2;

    if (alphaStar > alphaMax)
        alphaStar = alphaMax;

    if (alphaStar != alphaStar)
        alphaStar = 0.;

    //Take a final step to alphaStar
    x = x0 + alphaStar * xdir;
    subsys->setParams(x);

    return alphaStar;
}


void SubSystem::printResidual()
{
    Eigen::VectorXd r( constraints_list.size() );
    int i=0;
    double err = 0.;
    for (std::vector<Constraint *>::const_iterator constr=constraints_list.begin();
         constr != constraints_list.end(); ++constr, i++) {
        r[i] = (*constr)->error();
        err += r[i]*r[i];
    }
    err *= 0.5;
    std::cout << "Residual r = " << r.transpose() << std::endl;
    std::cout << "Residual err = " << err << std::endl;
}


} //namespace GCS_EXP
