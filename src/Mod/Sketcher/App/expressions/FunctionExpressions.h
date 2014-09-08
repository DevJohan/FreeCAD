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

#ifndef SKETCHER_FUNCTIONEXPRESSIONS_H_
#define SKETCHER_FUNCTIONEXPRESSIONS_H_

#include <Mod/Sketcher/App/expressions/ExpressionBase.h>

namespace SketcherExpressions{

static const int dynamic = -1;

template < int arg_count, ExpressionType::Types type, bool is_quantity>
class FunctionExpression: public ExpressionBase<type, is_quantity>{
public:
	typedef typename ExpressionBase<type, is_quantity>::value_type  value_type;

	static const int argument_count = arg_count;

	virtual ~FunctionExpression();

protected:
	Expression* arguments[argument_count];
};

template < ExpressionType::Types type, bool is_quantity>
class FunctionExpression<dynamic,type,is_quantity>: public ExpressionBase<type, is_quantity>{
public:
	typedef typename ExpressionBase< type, is_quantity >::value_type  value_type;

	virtual ~FunctionExpression();

protected:
	std::vector<Expression*> arguments;
};

template < bool is_quantity>
class SinExpression: public FunctionExpression<1,ExpressionType::Sin,is_quantity>{
public:
	typedef typename ExpressionBase<ExpressionType::Sin,is_quantity>::value_type  value_type;

	SinExpression(ExpressionReference argument);
	virtual ~SinExpression();
};

template < bool is_quantity>
class CosExpression: public FunctionExpression<1,ExpressionType::Cos,is_quantity>{
public:
	typedef typename ExpressionBase<ExpressionType::Cos,is_quantity>::value_type  value_type;

	CosExpression(ExpressionReference argument);
	virtual ~CosExpression();
};

template < bool is_quantity>
class TanExpression: public FunctionExpression<1,ExpressionType::Tan,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Tan,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	TanExpression(ExpressionReference argument);
	virtual ~TanExpression();
};

template < bool is_quantity>
class ArcsinExpression: public FunctionExpression<1,ExpressionType::Arcsin,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Arcsin,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	ArcsinExpression(ExpressionReference argument);
	virtual ~ArcsinExpression();
};

template < bool is_quantity>
class ArccosExpression: public FunctionExpression<1,ExpressionType::Arccos,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Arccos,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	ArccosExpression(ExpressionReference argument);
	virtual ~ArccosExpression();
};

template < bool is_quantity>
class ArctanExpression: public FunctionExpression<1,ExpressionType::Arctan,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Arctan,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	ArctanExpression(ExpressionReference argument);
	virtual ~ArctanExpression();
};

template < bool is_quantity>
class Arctan2Expression: public FunctionExpression<2, ExpressionType::Arctan2, is_quantity>{
public:
	typedef FunctionExpression<2,ExpressionType::Arctan2,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	Arctan2Expression(ExpressionReference y, ExpressionReference x);
	virtual ~Arctan2Expression();
};

//ABS
template < bool is_quantity>
class AbsExpression: public FunctionExpression<1,ExpressionType::Abs,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Abs,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	AbsExpression(ExpressionReference argument);
	virtual ~AbsExpression();
};

//POW
template < bool is_quantity>
class PowExpression: public FunctionExpression<2,ExpressionType::Pow,is_quantity>{
public:
	typedef FunctionExpression<2,ExpressionType::Pow,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	PowExpression(ExpressionReference base, ExpressionReference argument);
	virtual ~PowExpression();
};

//LOG
template < bool is_quantity>
class LogExpression: public FunctionExpression<2,ExpressionType::Log,is_quantity>{
public:
	typedef FunctionExpression<2,ExpressionType::Log,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	LogExpression(ExpressionReference base, ExpressionReference argument);
	virtual ~LogExpression();
};

//EXP
template < bool is_quantity>
class ExpExpression: public FunctionExpression<1,ExpressionType::Exp,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Exp,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	ExpExpression(ExpressionReference argument);
	virtual ~ExpExpression();
};

//LN
template < bool is_quantity>
class LnExpression: public FunctionExpression<1,ExpressionType::Ln,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Ln,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	LnExpression(ExpressionReference argument);
	virtual ~LnExpression();
};


//SINH
template < bool is_quantity>
class SinhExpression: public FunctionExpression<1,ExpressionType::Sinh,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Sinh,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	SinhExpression(ExpressionReference argument);
	virtual ~SinhExpression();
};

//COSH
template < bool is_quantity>
class CoshExpression: public FunctionExpression<1,ExpressionType::Cosh,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Cosh,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	CoshExpression(ExpressionReference argument);
	virtual ~CoshExpression();
};

//TANH
template < bool is_quantity>
class TanhExpression: public FunctionExpression<1,ExpressionType::Tanh,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Tanh,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	TanhExpression(ExpressionReference argument);
	virtual ~TanhExpression();
};

//SQRT
template < bool is_quantity>
class SqrtExpression: public FunctionExpression<1,ExpressionType::Sqrt,is_quantity>{
public:
	typedef FunctionExpression<1,ExpressionType::Sqrt,is_quantity>  base_type;
	typedef typename base_type::value_type  value_type;

	SqrtExpression(ExpressionReference argument);
	virtual ~SqrtExpression();
};

} // namespace SketcherExpressions

#endif /* FUNCTIONEXPRESSIONS_H_ */
