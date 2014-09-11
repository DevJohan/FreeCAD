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

#ifndef SKETCHER_PRODUCTEXPRESSION_H_
#define SKETCHER_PRODUCTEXPRESSION_H_

#include "FunctionExpressions.h"


namespace SketcherExpressions{

template < ExpressionBaseTypes data_value_type >
class ProductExpression: public ExpressionBase<ExpressionType::Product, data_value_type>{
public:
	typedef ExpressionBase<ExpressionType::Product, data_value_type> base_type;
	typedef typename base_type::value_type  value_type;
	typedef std::vector<ExpressionReference> argument_list_type;

	ProductExpression(
			const std::vector<ExpressionReference>& mult_factors,
			const std::vector<ExpressionReference>& div_factor);

	virtual ~ProductExpression();
	virtual void print(std::ostream& os, expr_print_modifier epm) const;
protected:
	argument_list_type numerator_arguments;
	argument_list_type denominator_arguments;
};


} // namespace SketcherExpressions


#endif /* SKETCHER_PRODUCTEXPRESSION_H_ */
