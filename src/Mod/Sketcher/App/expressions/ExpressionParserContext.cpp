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


#include "ExpressionParserContext.h"
#include <Mod/Sketcher/App/expressions/ExpressionLexer.h>
#include <Mod/Sketcher/App/expressions/ExpressionParser.hpp>

namespace SketcherExpressions{


ExpressionParserContext::LexerStateHolder::LexerStateHolder( std::string& str ):
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

int ExpressionParserContext::parse(const std::string &expression_string)
{
	LexerStateHolder lex_state;
	QuantityParser parser( lex_state, *this );
	int res = parser.parse();
	return res;
}

} // namespace SketcherExpressions
