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

#ifndef SKETCHER_EXPRESSIONCONTEXT_H_
#define SKETCHER_EXPRESSIONCONTEXT_H_

#include <vector>
#include <set>
#include <map>
#include <QString>
#include "ExpressionBase.h"

namespace SketcherExpressions{

struct ExpressionDescriptor{
	std::string _identifier_name;
	QString _display_name;
	ExpressionReference _expression;
	std::set<expression_index_type> _dependent_variables;
};

class ExpressionContext{
public:
	ExpressionContext();
	~ExpressionContext();

	ExpressionReference getExpression( std::string );

	expression_index_type getNamedIndex( const std::string& name ) const;
private:
	// Contains all the defined expressions in this context
	std::vector<ExpressionDescriptor> expressions;
	// Contains the list of indices of expressions variable which no longer
	// corresponds to an expression and can be overwritten.
	std::set<expression_index_type> empty_expression_slots;
	// Maps user usable names to corresponding expression.
	std::map< std::string, expression_index_type > named_expression_map;
};

} // namespace SketcherExpressions

#endif /* SKETCHER_EXPRESSIONCONTEXT_H_ */
