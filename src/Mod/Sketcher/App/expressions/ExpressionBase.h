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

#ifndef SKETCHER_EXPRESSIONBASE_H_
#define SKETCHER_EXPRESSIONBASE_H_

#include <string>
#include <vector>

#define BaseExport
#include <Base/Quantity.h>

namespace SketcherExpressions{
using Base::Quantity;

class ExpressionContext;

struct ExpressionType{
	enum Types{
		Constant,
		Named,
		ContextlessNamed,
		Sum,
		Product,
		Sin,
		Cos,
		Tan,
		Arcsin,
		Arccos,
		Arctan,
		Arctan2,
		Abs,
		Pow,
		Log,
		Exp,
		Ln,
		Sinh,
		Cosh,
		Tanh,
		Sqrt
	};
};

struct expr_print_modifier{
	enum Modifiers { add_outer_parenthesis };
	int flags;
	template < Modifiers mod >
	void set( bool value ){
		if(value){
			flags |= ( 1 << mod );
		}else{
			flags &= ~(1 << mod);
		}
	}
	template < Modifiers mod>
	bool get(){
		return flags & ( 1 << mod );
	}
};

class Expression;
typedef Expression* ExpressionReference;

typedef size_t expression_index_type;

class Expression{
public:
	virtual ~Expression() = 0;

	/* Return the type of this expression.
	 */
	virtual ExpressionType getType() const = 0;
	/* Creates a copy of the current expression.
	 */
	virtual ExpressionReference clone() const = 0;
	/* Creates a dimensionless copy of the current expression. Expressions
	 * which carry a Quantity should give the value in the Unit used in the
	 * internal representation.
	 */
	virtual ExpressionReference dimensionlessCopy() const = 0;
	/* Returns true if this expression is dependent on a ExpressionContext.
	 */
	virtual bool isContextDependent() const = 0;
	/* A expression is constant if it's value doesn't depend on
	 * a named non-constant value.
	 */
	virtual bool isConstant() const = 0;
	/* Returns the named variables that affect this expression.
	 */
	virtual std::vector<expression_index_type> getDependencies() const = 0;
	/* Returns a textual representation which can be used to recreate
	 * the current expression exactly.
	 */
	virtual std::string stringRepresentation() const = 0;
	/* Returns a visual description of this expression as a string.
	 * This doesn't have to be able to recreate the expression.
	 */
	std::string print() const;
	/* Prints visual description of this expression to ostream.
	 * This doesn't have to be able to recreate the expression.
	 */
	virtual void print( std::ostream& os, expr_print_modifier epm=expr_print_modifier() ) const = 0;
	/* Prints visual description of this expression to ostream.
	 * This doesn't have to be able to recreate the expression.
	 */
	virtual bool needParenthesis() const { return true; }
};

std::ostream& operator<<( std::ostream& os, const Expression& expr );

template <ExpressionType::Types type, bool is_quantity>
class ExpressionBase;

template <ExpressionType::Types _type>
class ExpressionBase<_type, true>: public Expression{
public:
	static const ExpressionType::Types type = _type;
	typedef Quantity value_type;
	void print(std::ostream& os, expr_print_modifier epm) const;
	virtual ExpressionType getType() const;
	virtual ExpressionReference clone() const;
	virtual ExpressionReference dimensionlessCopy() const;
	virtual bool isContextDependent() const;
	virtual bool isConstant() const;
	virtual std::vector<expression_index_type> getDependencies() const;
	virtual std::string stringRepresentation() const;
};

template <ExpressionType::Types _type>
class ExpressionBase<_type, false>: public Expression{
public:
	static const ExpressionType::Types type = _type;
	typedef double value_type;
	void print(std::ostream& os, expr_print_modifier epm) const;
	virtual ExpressionType getType() const;
	virtual ExpressionReference clone() const;
	virtual ExpressionReference dimensionlessCopy() const;
	virtual bool isContextDependent() const;
	virtual bool isConstant() const;
	virtual std::vector<expression_index_type> getDependencies() const;
	virtual std::string stringRepresentation() const;
};


template <bool is_quantity>
class ConstantExpression: public ExpressionBase< ExpressionType::Constant,is_quantity>{
public:
	typedef typename ExpressionBase< ExpressionType::Constant,is_quantity>::value_type  value_type;

	ConstantExpression( value_type value ): _value( value ){ }
	ConstantExpression( ExpressionReference value ): _value( 0 ){ }
	virtual ~ConstantExpression(){ }
	virtual std::vector<expression_index_type> getDependencies() const{ return std::vector<expression_index_type>();};
protected:
	value_type _value;
};

template <bool is_quantity>
class NamedExpression: public ExpressionBase<ExpressionType::Named, is_quantity>{
public:
	typedef typename ExpressionBase<ExpressionType::Named, is_quantity>::value_type  value_type;

	NamedExpression( ExpressionContext& context, const std::string& name );
	virtual ~NamedExpression(){ }
	virtual std::vector<ExpressionReference> getDependencies() const;
protected:
	ExpressionContext& _context;
	expression_index_type _index;
};

template <bool is_quantity>
class ContextlessNamedExpression
		: public ExpressionBase<ExpressionType::ContextlessNamed, is_quantity>{
public:
	typedef typename ExpressionBase<ExpressionType::ContextlessNamed, is_quantity>::value_type  value_type;

	virtual ~ContextlessNamedExpression(){ }
	virtual std::vector<ExpressionReference> getDependencies() const;
protected:
	std::string _name;
};

} // namespace SketcherExpressions

#endif /* SKETCHER_EXPRESSIONBASE_H_ */
