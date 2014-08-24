/*
 * ExpressionParserContext.h
 *
 *  Created on: 22 aug 2014
 *      Author: johan
 */

#ifndef SKETCHER_QUANTITYPARSERCONTEXT_H_
#define SKETCHER_QUANTITYPARSERCONTEXT_H_

#include "ExpressionParserContext.h"
#include <string>
#include <map>
#include <Mod/Sketcher/App/expressions/ExpressionLexer.h>
#include <Mod/Sketcher/App/expressions/ExpressionParser.hpp>
#include <Base/Exception.h>

#define YY_DECL ExpressionParser::symbol_type yylex( yyscan_t& yyscanner,  ExpressionParserContext& context, ExpressionParser::location_type& location )

// Data from http://www.itl.nist.gov/div898/strd/univ/data/PiDigits.dat
#define MATH_CONSTANT_PI ( 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317 )
// Data from http://apod.nasa.gov/htmltest/gifcity/e.2mil
#define MATH_CONSTANT_E  ( 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729 )

namespace SketcherExpressions{


class SketcherExport ExpressionParserContext{
protected:
	class LexerStateHolder{
	public:
		LexerStateHolder( std::string& str );
		~LexerStateHolder();
		operator yyscan_t (){ return _state; }
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

	int parse (const std::string& f);
};

} // namespace SketcherExpressions

#endif /* SKETCHER_QUANTITYPARSERCONTEXT_H_ */
