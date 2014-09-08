/*
 * ExpressionParserContext.h
 *
 *  Created on: 22 aug 2014
 *      Author: johan
 */

#ifndef SKETCHER_QUANTITYPARSERCONTEXT_H_
#define SKETCHER_QUANTITYPARSERCONTEXT_H_

#include <string>
#include <map>

namespace SketcherExpressions{
	class ExpressionParserContext;
}
#include <Mod/Sketcher/App/expressions/ExpressionParser.hpp>
#include <Mod/Sketcher/App/expressions/ExpressionLexer.h>
#include <Base/Exception.h>

#include <Mod/Sketcher/App/expressions/ExpressionBase.h>


// Data from http://www.itl.nist.gov/div898/strd/univ/data/PiDigits.dat
#define MATH_CONSTANT_PI ( 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317 )
// Data from http://apod.nasa.gov/htmltest/gifcity/e.2mil
#define MATH_CONSTANT_E  ( 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729 )

namespace SketcherExpressions{


class SketcherExport ExpressionParserContext{
protected:
	class LexerStateHolder{
	public:
		LexerStateHolder( const std::string& str );
		~LexerStateHolder();
		operator yyscan_t& (){ return _state; }
	private:
		LexerStateHolder(const LexerStateHolder&);
		LexerStateHolder& operator=(const LexerStateHolder&);
		yyscan_t _state;
		YY_BUFFER_STATE _buffer;
	};
public:
	ExpressionParserContext ();
	virtual ~ExpressionParserContext ();

	std::map<std::string, int> variables;

	Quantity result;

	ExpressionReference create_parameter( std::string& name );
	ExpressionReference create_const( double number );
	ExpressionReference create_const( Quantity value );
	ExpressionReference create_const(      ExpressionReference number   );
	ExpressionReference create_addition( ExpressionReference expr1, ExpressionReference expr2 );
	ExpressionReference create_subtraction( ExpressionReference expr1, ExpressionReference expr2 );
	ExpressionReference create_multiplication( ExpressionReference expr1, ExpressionReference expr2 );
	ExpressionReference create_division( ExpressionReference expr1, ExpressionReference expr2 );
	ExpressionReference create_negation( ExpressionReference expr);
	ExpressionReference create_acos_expr(  ExpressionReference argument );
	ExpressionReference create_asin_expr(  ExpressionReference argument );
	ExpressionReference create_atan_expr(  ExpressionReference argument );
	ExpressionReference create_atan2_expr( ExpressionReference y, ExpressionReference x );
	ExpressionReference create_abs_expr(   ExpressionReference argument );
	ExpressionReference create_pow_expr(   ExpressionReference base, ExpressionReference exponent );
	ExpressionReference create_exp_expr(   ExpressionReference argument );
	ExpressionReference create_ln_expr(    ExpressionReference argument );
	ExpressionReference create_log_expr( double base, ExpressionReference argument ){ return new ConstantExpression<false>(7);}
	ExpressionReference create_log_expr( ExpressionReference base, ExpressionReference argument );
	ExpressionReference create_sin_expr(   ExpressionReference argument );
	ExpressionReference create_sinh_expr(  ExpressionReference argument );
	ExpressionReference create_cos_expr(   ExpressionReference argument );
	ExpressionReference create_cosh_expr(   ExpressionReference argument );
	ExpressionReference create_tan_expr(   ExpressionReference argument );
	ExpressionReference create_tanh_expr(  ExpressionReference argument );
	ExpressionReference create_sqrt_expr(  ExpressionReference argument );

	int parse (const std::string& f);
};

} // namespace SketcherExpressions

#endif /* SKETCHER_QUANTITYPARSERCONTEXT_H_ */
