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

template < int arg_count, ExpressionType::Types type, ExpressionBaseTypes data_value_type>
class FunctionExpression: public ExpressionBase<type, data_value_type>{
public:
	typedef typename ExpressionBase<type, data_value_type>::value_type  value_type;

	static const int argument_count = arg_count;

	FunctionExpression();
	virtual ~FunctionExpression();

protected:
	Expression* arguments[argument_count];
};

template < ExpressionType::Types type, ExpressionBaseTypes data_value_type>
class FunctionExpression<dynamic,type,data_value_type>: public ExpressionBase<type, data_value_type>{
public:
	typedef typename ExpressionBase< type, data_value_type >::value_type  value_type;

	FunctionExpression();
	virtual ~FunctionExpression();

protected:
	std::vector<Expression*> arguments;
};

template < ExpressionBaseTypes data_value_type>
class SinExpression: public FunctionExpression<1,ExpressionType::Sin,data_value_type>{
public:
	typedef typename ExpressionBase<ExpressionType::Sin,data_value_type>::value_type  value_type;

	SinExpression(ExpressionReference argument);
	virtual ~SinExpression();
};

template < ExpressionBaseTypes data_value_type>
class CosExpression: public FunctionExpression<1,ExpressionType::Cos,data_value_type>{
public:
	typedef typename ExpressionBase<ExpressionType::Cos,data_value_type>::value_type  value_type;

	CosExpression(ExpressionReference argument);
	virtual ~CosExpression();
};

template < ExpressionBaseTypes data_value_type>
class TanExpression: public FunctionExpression<1,ExpressionType::Tan,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Tan,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	TanExpression(ExpressionReference argument);
	virtual ~TanExpression();
};

template < ExpressionBaseTypes data_value_type>
class ArcsinExpression: public FunctionExpression<1,ExpressionType::Arcsin,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Arcsin,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	ArcsinExpression(ExpressionReference argument);
	virtual ~ArcsinExpression();
};

template < ExpressionBaseTypes data_value_type>
class ArccosExpression: public FunctionExpression<1,ExpressionType::Arccos,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Arccos,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	ArccosExpression(ExpressionReference argument);
	virtual ~ArccosExpression();
};

template < ExpressionBaseTypes data_value_type>
class ArctanExpression: public FunctionExpression<1,ExpressionType::Arctan,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Arctan,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	ArctanExpression(ExpressionReference argument);
	virtual ~ArctanExpression();
};

template < ExpressionBaseTypes data_value_type>
class Arctan2Expression: public FunctionExpression<2, ExpressionType::Arctan2, data_value_type>{
public:
	typedef FunctionExpression<2,ExpressionType::Arctan2,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	Arctan2Expression(ExpressionReference y, ExpressionReference x);
	virtual ~Arctan2Expression();
};

//ABS
template < ExpressionBaseTypes data_value_type>
class AbsExpression: public FunctionExpression<1,ExpressionType::Abs,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Abs,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	AbsExpression(ExpressionReference argument);
	virtual ~AbsExpression();
};

//POW
template < ExpressionBaseTypes data_value_type>
class PowExpression: public FunctionExpression<2,ExpressionType::Pow,data_value_type>{
public:
	typedef FunctionExpression<2,ExpressionType::Pow,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	PowExpression(ExpressionReference base, ExpressionReference argument);
	virtual ~PowExpression();
};

//LOG
template < ExpressionBaseTypes data_value_type>
class LogExpression: public FunctionExpression<2,ExpressionType::Log,data_value_type>{
public:
	typedef FunctionExpression<2,ExpressionType::Log,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	LogExpression(ExpressionReference base, ExpressionReference argument);
	virtual ~LogExpression();
};

//EXP
template < ExpressionBaseTypes data_value_type>
class ExpExpression: public FunctionExpression<1,ExpressionType::Exp,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Exp,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	ExpExpression(ExpressionReference argument);
	virtual ~ExpExpression();
};

//LN
template < ExpressionBaseTypes data_value_type>
class LnExpression: public FunctionExpression<1,ExpressionType::Ln,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Ln,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	LnExpression(ExpressionReference argument);
	virtual ~LnExpression();
};


//SINH
template < ExpressionBaseTypes data_value_type>
class SinhExpression: public FunctionExpression<1,ExpressionType::Sinh,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Sinh,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	SinhExpression(ExpressionReference argument);
	virtual ~SinhExpression();
};

//COSH
template < ExpressionBaseTypes data_value_type>
class CoshExpression: public FunctionExpression<1,ExpressionType::Cosh,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Cosh,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	CoshExpression(ExpressionReference argument);
	virtual ~CoshExpression();
};

//TANH
template < ExpressionBaseTypes data_value_type>
class TanhExpression: public FunctionExpression<1,ExpressionType::Tanh,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Tanh,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	TanhExpression(ExpressionReference argument);
	virtual ~TanhExpression();
};

//SQRT
template < ExpressionBaseTypes data_value_type>
class SqrtExpression: public FunctionExpression<1,ExpressionType::Sqrt,data_value_type>{
public:
	typedef FunctionExpression<1,ExpressionType::Sqrt,data_value_type>  base_type;
	typedef typename base_type::value_type  value_type;

	SqrtExpression(ExpressionReference argument);
	virtual ~SqrtExpression();
};

} // namespace SketcherExpressions

#endif /* FUNCTIONEXPRESSIONS_H_ */
