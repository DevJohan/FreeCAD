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

#include "ExpressionParserContext.h"

#include <Mod/Sketcher/App/expressions/ExpressionLexer.h>
#include <Mod/Sketcher/App/expressions/ExpressionParser.hpp>

#include <Mod/Sketcher/App/expressions/ExpressionBase.h>
#include <Mod/Sketcher/App/expressions/FunctionExpressions.h>


namespace SketcherExpressions{


ExpressionParserContext::LexerStateHolder::LexerStateHolder( const std::string& str ):
																				_state(0),
																				_buffer(0)
{
	if( yylex_init( &_state ) != 0 ){
		_state = 0;
		return;
	}
	_buffer = yy_scan_bytes( str.c_str(),str.size(), _state );

}
ExpressionParserContext::LexerStateHolder::~LexerStateHolder(){
	if( _state ){
		yy_delete_buffer(_buffer, _state);
		yylex_destroy( &_state );
		_state = 0;
	}
}


ExpressionParserContext::ExpressionParserContext ()
{
}

ExpressionParserContext::~ExpressionParserContext ()
{
}

int ExpressionParserContext::parse( const std::string &expression_string )
{
	LexerStateHolder lex_state( expression_string );
	ExpressionParser::location_type location;
	ExpressionParser parser( lex_state, *this, location );
	int res = parser.parse();
	return res;
}

bool isQuantityExpression(ExpressionReference number){ return true; }

ExpressionReference ExpressionParserContext::create_const(
		double number
){
	return new ConstantExpression<false>( number );
}

ExpressionReference ExpressionParserContext::create_const(
		Quantity value
){
	return new ConstantExpression<true>( value );
}

ExpressionReference ExpressionParserContext::create_const(
		ExpressionReference number
){
	return isQuantityExpression(number) ?
			static_cast<ExpressionReference>(new ConstantExpression<true>( number )):
			static_cast<ExpressionReference>(new ConstantExpression<false>( number ));
}

ExpressionReference ExpressionParserContext::create_acos_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new ArccosExpression<true>( argument )):
					static_cast<ExpressionReference>(new ArccosExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_asin_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new ArcsinExpression<true>( argument )):
					static_cast<ExpressionReference>(new ArcsinExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_atan_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new ArctanExpression<true>( argument )):
					static_cast<ExpressionReference>(new ArctanExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_atan2_expr(
		ExpressionReference y,
		ExpressionReference x
){
	return isQuantityExpression(y) ?
			static_cast<ExpressionReference>(new Arctan2Expression<true>( y, x )):
					static_cast<ExpressionReference>(new Arctan2Expression<false>( y, x ));
}

ExpressionReference ExpressionParserContext::create_abs_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new AbsExpression<true>( argument )):
					static_cast<ExpressionReference>(new AbsExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_pow_expr(
		ExpressionReference base,
		ExpressionReference exponent
){
	return isQuantityExpression(base) ?
			static_cast<ExpressionReference>(new PowExpression<true>( base, exponent )):
					static_cast<ExpressionReference>(new PowExpression<false>( base, exponent ));
}

ExpressionReference ExpressionParserContext::create_exp_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new ExpExpression<true>( argument )):
					static_cast<ExpressionReference>(new ExpExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_ln_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new LnExpression<true>( argument )):
					static_cast<ExpressionReference>(new LnExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_log_expr(
		ExpressionReference base,
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new LogExpression<true>( base, argument )):
			static_cast<ExpressionReference>(new LogExpression<false>( base, argument ));
}

ExpressionReference ExpressionParserContext::create_sin_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new SinExpression<true>( argument )):
					static_cast<ExpressionReference>(new SinExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_sinh_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new SinhExpression<true>( argument )):
					static_cast<ExpressionReference>(new SinhExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_cos_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new CosExpression<true>( argument )):
					static_cast<ExpressionReference>(new CosExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_cosh_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new CoshExpression<true>( argument )):
					static_cast<ExpressionReference>(new CoshExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_tan_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new TanExpression<true>( argument )):
					static_cast<ExpressionReference>(new TanExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_tanh_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new TanhExpression<true>( argument )):
					static_cast<ExpressionReference>(new TanhExpression<false>( argument ));
}

ExpressionReference ExpressionParserContext::create_sqrt_expr(
		ExpressionReference argument
){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new SqrtExpression<true>( argument )):
					static_cast<ExpressionReference>(new SqrtExpression<false>( argument ));
}



} // namespace SketcherExpressions
