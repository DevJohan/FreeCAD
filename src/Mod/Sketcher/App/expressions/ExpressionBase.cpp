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



// ExpressionBase<_type, plain>

template <ExpressionType::Types _type>
void ExpressionBase<_type, plain>::print(std::ostream& os, expr_print_modifier epm) const{

}

template <ExpressionType::Types _type>
ExpressionType::Types ExpressionBase<_type, plain>::getType() const{
	return type;
}

template <ExpressionType::Types _type>
ExpressionReference ExpressionBase<_type, plain>::clone() const{

}

template <ExpressionType::Types _type>
ExpressionReference ExpressionBase<_type, plain>::dimensionlessCopy() const{

}

template <ExpressionType::Types _type>
bool ExpressionBase<_type, plain>::isContextDependent() const{

}

template <ExpressionType::Types _type>
bool ExpressionBase<_type, plain>::isConstant() const{

}

template <ExpressionType::Types _type>
std::vector<expression_index_type> ExpressionBase<_type, plain>::getDependencies() const{

}

template <ExpressionType::Types _type>
std::string ExpressionBase<_type, plain>::stringRepresentation() const{

}


// ExpressionBase<_type, quantity>
template <ExpressionType::Types _type>
void ExpressionBase<_type, quantity>::print(std::ostream& os, expr_print_modifier epm) const{

}

template <ExpressionType::Types _type>
ExpressionType::Types ExpressionBase<_type, quantity>::getType() const{
	return type;
}

template <ExpressionType::Types _type>
ExpressionReference ExpressionBase<_type, quantity>::clone() const{

}

template <ExpressionType::Types _type>
ExpressionReference ExpressionBase<_type, quantity>::dimensionlessCopy() const{

}

template <ExpressionType::Types _type>
bool ExpressionBase<_type, quantity>::isContextDependent() const{

}

template <ExpressionType::Types _type>
bool ExpressionBase<_type, quantity>::isConstant() const{

}

template <ExpressionType::Types _type>
std::vector<expression_index_type> ExpressionBase<_type, quantity>::getDependencies() const{

}

template <ExpressionType::Types _type>
std::string ExpressionBase<_type, quantity>::stringRepresentation() const{

}


template <ExpressionBaseTypes data_value_type>
NamedExpression<data_value_type>::NamedExpression( ExpressionContext& context, const std::string& name )
:_context(context),
 _index(_context.getNamedIndex( name ))
{ }


// ContextlessNamedExpression

template <ExpressionBaseTypes data_value_type>
ContextlessNamedExpression<data_value_type>::ContextlessNamedExpression(const std::string& name)
:_name(name)
{

}

template <ExpressionBaseTypes data_value_type>
std::vector<expression_index_type> ContextlessNamedExpression<data_value_type>::getDependencies() const {
	return std::vector<expression_index_type>();
}

// Template instantiation
template class ExpressionBase<ExpressionType::Constant, plain>;
template class ExpressionBase<ExpressionType::Constant, quantity>;

template class ExpressionBase<ExpressionType::Named, plain>;
template class ExpressionBase<ExpressionType::Named, quantity>;

template class ExpressionBase<ExpressionType::ContextlessNamed, plain>;
template class ExpressionBase<ExpressionType::ContextlessNamed, quantity>;

template class ExpressionBase<ExpressionType::Sum, plain>;
template class ExpressionBase<ExpressionType::Sum, quantity>;
template class ExpressionBase<ExpressionType::Product, plain>;
template class ExpressionBase<ExpressionType::Product, quantity>;

template class ExpressionBase<ExpressionType::Sin, plain>;
template class ExpressionBase<ExpressionType::Sin, quantity>;
template class ExpressionBase<ExpressionType::Cos, plain>;
template class ExpressionBase<ExpressionType::Cos, quantity>;
template class ExpressionBase<ExpressionType::Tan, plain>;
template class ExpressionBase<ExpressionType::Tan, quantity>;
template class ExpressionBase<ExpressionType::Arcsin, plain>;
template class ExpressionBase<ExpressionType::Arcsin, quantity>;
template class ExpressionBase<ExpressionType::Arccos, plain>;
template class ExpressionBase<ExpressionType::Arccos, quantity>;
template class ExpressionBase<ExpressionType::Arctan, plain>;
template class ExpressionBase<ExpressionType::Arctan, quantity>;
template class ExpressionBase<ExpressionType::Arctan2, plain>;
template class ExpressionBase<ExpressionType::Arctan2, quantity>;
template class ExpressionBase<ExpressionType::Abs, plain>;
template class ExpressionBase<ExpressionType::Abs, quantity>;
template class ExpressionBase<ExpressionType::Pow, plain>;
template class ExpressionBase<ExpressionType::Pow, quantity>;
template class ExpressionBase<ExpressionType::Log, plain>;
template class ExpressionBase<ExpressionType::Log, quantity>;
template class ExpressionBase<ExpressionType::Exp, plain>;
template class ExpressionBase<ExpressionType::Exp, quantity>;
template class ExpressionBase<ExpressionType::Ln, plain>;
template class ExpressionBase<ExpressionType::Ln, quantity>;
template class ExpressionBase<ExpressionType::Sinh, plain>;
template class ExpressionBase<ExpressionType::Sinh, quantity>;
template class ExpressionBase<ExpressionType::Cosh, plain>;
template class ExpressionBase<ExpressionType::Cosh, quantity>;
template class ExpressionBase<ExpressionType::Tanh, plain>;
template class ExpressionBase<ExpressionType::Tanh, quantity>;
template class ExpressionBase<ExpressionType::Sqrt, plain>;
template class ExpressionBase<ExpressionType::Sqrt, quantity>;

template class ContextlessNamedExpression<plain>;
template class ContextlessNamedExpression<quantity>;

} // namespace SketcherExpressions
