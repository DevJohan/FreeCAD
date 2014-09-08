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

#include "PreCompiled.h"

#include "ExpressionBase.h"
#include <sstream>

#include "ExpressionContext.h"

namespace SketcherExpressions{


Expression::~Expression(){ }

std::string Expression::print() const {
	std::ostringstream oss;
	print(oss);
	return oss.str();
}

std::ostream& operator<<( std::ostream& os, const Expression& expr ){
	expr.print(os);
	return os;
}

template <bool is_quantity>
NamedExpression<is_quantity>::NamedExpression( ExpressionContext& context, const std::string& name )
:_context(context),
 _index(_context.getNamedIndex( name ))
{ }

} // namespace SketcherExpressions
