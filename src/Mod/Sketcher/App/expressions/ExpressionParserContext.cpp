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
#include <Mod/Sketcher/App/expressions/SumExpression.h>
#include <Mod/Sketcher/App/expressions/ProductExpression.h>
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

template< template <ExpressionBaseTypes> class T >
static ExpressionReference create_expr(ExpressionReference argument){
	return isQuantityExpression(argument) ?
			static_cast<ExpressionReference>(new T<quantity>( argument )):
			static_cast<ExpressionReference>(new T<plain>( argument ));
}

template< template <ExpressionBaseTypes> class T >
static ExpressionReference create_expr(ExpressionReference argument1, ExpressionReference argument2){
	return isQuantityExpression(argument1) ?
			static_cast<ExpressionReference>(new T<quantity>( argument1, argument2 )):
			static_cast<ExpressionReference>(new T<plain>( argument1, argument2 ));
}

template< template <ExpressionBaseTypes> class T >
static ExpressionReference create_expr(
		const std::vector<ExpressionReference>& argument1,
		const std::vector<ExpressionReference>& argument2 = std::vector<ExpressionReference>()){
	return isQuantityExpression(argument1.front()) ?
			static_cast<ExpressionReference>(new T<quantity>( argument1, argument2 )):
			static_cast<ExpressionReference>(new T<plain>( argument1, argument2 ));
}

ExpressionReference ExpressionParserContext::create_parameter( std::string& name ){
	return new ContextlessNamedExpression<quantity>( name );
}

ExpressionReference ExpressionParserContext::create_const(
		double number
){
	return new ConstantExpression<plain>( number );
}

ExpressionReference ExpressionParserContext::create_const(
		Quantity value
){
	return new ConstantExpression<quantity>( value );
}

ExpressionReference ExpressionParserContext::create_const(
		ExpressionReference number
){
	return create_expr<ConstantExpression>( number );
}

ExpressionReference ExpressionParserContext::create_addition( ExpressionReference expr1, ExpressionReference expr2 ){
	std::vector<ExpressionReference> arguments( 2 );
	arguments.push_back(expr1);
	arguments.push_back(expr2);
	return create_expr<SumExpression>( arguments );
}
ExpressionReference ExpressionParserContext::create_subtraction( ExpressionReference expr1, ExpressionReference expr2 ){
	std::vector<ExpressionReference> arguments1( 1 );
	arguments1.push_back(expr1);
	std::vector<ExpressionReference> arguments2( 1 );
	arguments2.push_back(expr2);
	return create_expr<SumExpression>( arguments1, arguments2 );
}
ExpressionReference ExpressionParserContext::create_multiplication( ExpressionReference expr1, ExpressionReference expr2 ){
	std::vector<ExpressionReference> arguments( 2 );
	arguments.push_back(expr1);
	arguments.push_back(expr2);
	return create_expr<ProductExpression>( arguments );
}
ExpressionReference ExpressionParserContext::create_division( ExpressionReference expr1, ExpressionReference expr2 ){
	std::vector<ExpressionReference> arguments1( 1 );
	arguments1.push_back(expr1);
	std::vector<ExpressionReference> arguments2( 1 );
	arguments2.push_back(expr2);
	return create_expr<ProductExpression>( arguments1, arguments2 );
}
ExpressionReference ExpressionParserContext::create_negation( ExpressionReference expr){
	std::vector<ExpressionReference> arguments( 1 );
	arguments.push_back(expr);
	return create_expr<SumExpression>( std::vector<ExpressionReference>(), arguments );
}

ExpressionReference ExpressionParserContext::create_acos_expr(
		ExpressionReference argument
){
	return create_expr<ArccosExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_asin_expr(
		ExpressionReference argument
){
	return create_expr<ArcsinExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_atan_expr(
		ExpressionReference argument
){
	return create_expr<ArctanExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_atan2_expr(
		ExpressionReference y,
		ExpressionReference x
){
	return create_expr<Arctan2Expression>( y, x );
}

ExpressionReference ExpressionParserContext::create_abs_expr(
		ExpressionReference argument
){
	return create_expr<AbsExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_pow_expr(
		ExpressionReference base,
		ExpressionReference exponent
){
	return create_expr<PowExpression>( base, exponent );
}

ExpressionReference ExpressionParserContext::create_exp_expr(
		ExpressionReference argument
){
	return create_expr<ExpExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_ln_expr(
		ExpressionReference argument
){
	return create_expr<LnExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_log_expr(
		ExpressionReference base,
		ExpressionReference argument
){
	return create_expr<LogExpression>( new ConstantExpression<plain>(10.0), argument );
}

ExpressionReference ExpressionParserContext::create_sin_expr(
		ExpressionReference argument
){
	return create_expr<SinExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_sinh_expr(
		ExpressionReference argument
){
	return create_expr<SinhExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_cos_expr(
		ExpressionReference argument
){
	return create_expr<CosExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_cosh_expr(
		ExpressionReference argument
){
	return create_expr<CoshExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_tan_expr(
		ExpressionReference argument
){
	return create_expr<TanExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_tanh_expr(
		ExpressionReference argument
){
	return create_expr<TanhExpression>( argument );
}

ExpressionReference ExpressionParserContext::create_sqrt_expr(
		ExpressionReference argument
){
	return create_expr<SqrtExpression>( argument );
}



} // namespace SketcherExpressions
