/***************************************************************************
 *   Copyright (c) 2014 Johan Kristensen                                   *
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

#include "ProductExpression.h"
#include <sstream>


namespace SketcherExpressions{

template < bool is_quantity >
ProductExpression< is_quantity >::~ProductExpression(){

}


template < bool is_quantity >
void ProductExpression< is_quantity >::print( std::ostream& os, expr_print_modifier epm ) const {
	if( numerator_arguments.size() == 0 ){
		os << "1";
	} else{
		argument_list_type::iterator it = numerator_arguments.begin();
		os << "( ";
		if( (*it)->needParenthesis() ){
			os << "( ";
			(*it)->print(os);
			os << " )";
		}else{
			(*it)->print(os);
		}
		++it;
		const argument_list_type::iterator numer_it_end = numerator_arguments.end();
		for( ;it != numer_it_end; ++it ){
			os << "*";
			if((*it)->needParenthesis()){
				os << "( ";
				(*it)->print(os);
				os << " )";
			}else{
				(*it)->print(os);
			}
		}
		os <<")";
	} // if( numerator_arguments.size == 0 )
	if( denominator_arguments.size() > 0 ){
		argument_list_type::iterator it = denominator_arguments.begin();
		os << "/(";
		if( (*it)->needParenthesis() ){
			os << "( ";
			(*it)->print(os);
			os << " )";
		}else{
			(*it)->print(os);
		}
		++it;
		const argument_list_type::iterator denom_it_end = denominator_arguments.end();
		for( ;it != denom_it_end; ++it ){
			os << "*";
			if((*it)->needParenthesis()){
				os << "( ";
				(*it)->print(os);
				os << " )";
			}else{
				(*it)->print(os);
			}
		}
		os <<")";
	}
	return os;
}

} // namespace SketcherExpressions

