/* Parser for the FreeCAD  Units language           */
/* (c) 2013 Juergen Riegel  LGPL                    */

/* Bison options.  */
%defines
%language "C++"
%define api.namespace { SketcherExpressions }
%define api.value.type variant
%define api.token.constructor
%define api.location.type { linear_location }
%define parser_class_name { ExpressionParser }

%code requires {
	#include "PreCompiled.h"
	#ifndef _PreComp_
	# include <sstream>
	#endif

    #define BaseExport
    #include <Base/Unit.h>
    #include <Base/Quantity.h>
    #include <Mod/Sketcher/App/expressions/ScaledUnit.h>
    #include <Mod/Sketcher/App/expressions/ExpressionBase.h>
    using Base::Unit;
    using Base::Quantity;
    using SketcherExpressions::ScaledUnit;

	namespace SketcherExpressions {
		class ExpressionParserContext;
		typedef struct{ int begin; int end; } linear_location;
		inline void reset_location(linear_location& location){ location.begin = location.end; }
		inline void set_location(linear_location& location, int size){ location.end = location.begin + size; }
	}
	
	#define YY_DECL \
		SketcherExpressions::ExpressionParser::symbol_type yylex( \
			yyscan_t& yyscanner,  \
			SketcherExpressions::ExpressionParserContext& context, \
			SketcherExpressions::ExpressionParser::location_type& location )

	#include <Mod/Sketcher/App/expressions/ExpressionLexer.h>
			

}

%locations
%param { yyscan_t& yyscanner }
%param { ExpressionParserContext& context }
%param { location_type& location }
%initial-action{
}

/* Represents the many different ways we can access our data */
%{
        #ifndef  DOUBLE_MAX
        # define DOUBLE_MAX 1.7976931348623157E+308    /* max decimal value of a "double"*/
        #endif
        #ifndef  DOUBLE_MIN
        # define DOUBLE_MIN 2.2250738585072014E-308    /* min decimal value of a "double"*/
        #endif
%}

%code {
	#include <Mod/Sketcher/App/expressions/ExpressionParserContext.h>
	// Declare lexer prototype for parser to use
	YY_DECL;

	#include <limits>
	template <typename T>
    bool is_integer(T value){
    	const T result = abs(std::fmod(value,1.0)); 
    	return ( result > .5 ? 1.0 - result : result ) < std::numeric_limits<T>::epsilon();
    }
}

/* Bison declarations.  */
%token <ScaledUnit> UNIT "unit";
%token <double> NUM "floating point number";
%token <std::string> IDENTIFIER "identifier";
%token <char> OPERATOR_CHAR "operator char";
%token MINUSSIGN;
%token ACOS ASIN ATAN COS COSH EXP ABS LOG LOG10 SIN SINH TAN TANH SQRT;
%token END 0 "the end"
%left MINUSSIGN '+'
%left '*' '/'
%left NEG     /* negation--unary minus */
%right '^'    /* exponentiation */

%{
        #include <cmath>
               
%}

%type <ExpressionReference> expr
%type <double> num
%type <ScaledUnit> unit
%type <Quantity> quantity

%start expr

%%

     expr:    num                			{ $$ = context.create_const( $1 );     	}
             | quantity    			    	{ $$ = context.create_const( $1 );     	}
             | IDENTIFIER		    		{ $$ = context.create_parameter( $1 ); 	}
             | expr '+' expr        		{ $$ = context.create_addition( $1, $3);}
             | expr MINUSSIGN expr          { $$ = context.create_subtraction( $1, $3); }
             | expr '*' expr       			{ $$ = context.create_multiplication( $1, $3);}
             | expr '/' expr       			{ $$ = context.create_division( $1, $3);} 
             | MINUSSIGN expr  %prec NEG    { $$ = context.create_negation( $2);        	       				}
             | expr '^' expr       			{ $$ = context.create_pow_expr($1, $3);	}
             | '(' expr ')'        			{ $$ = $2;         						}
             | ACOS  '(' expr ')'  			{ $$ = context.create_acos_expr(  $3  ); }
             | ASIN  '(' expr ')'  			{ $$ = context.create_asin_expr(  $3  ); }
             | ATAN  '(' expr ')'  			{ $$ = context.create_atan_expr(  $3  ); }
             | ATAN  '(' expr '\\' expr ')'	{ $$ = context.create_atan2_expr( $5 , $3  ); }
             | ABS  '(' expr ')'   			{ $$ = context.create_abs_expr(   $3  ); }
             | EXP  '(' expr ')'   			{ $$ = context.create_exp_expr(   $3  ); }
             | LOG  '(' expr ')'			{ $$ = context.create_ln_expr(    $3  ); }
             | LOG10  '(' expr ')'			{ $$ = context.create_log_expr( 10.0, $3  ); }
             | SIN  '(' expr ')'   			{ $$ = context.create_sin_expr(   $3  ); }
             | SINH '(' expr ')'   			{ $$ = context.create_sinh_expr(  $3  ); }
             | COS  '(' expr ')'   			{ $$ = context.create_cos_expr(   $3  ); }
             | COSH  '(' expr ')'   		{ $$ = context.create_cosh_expr(  $3  ); }
             | TAN  '(' expr ')'   			{ $$ = context.create_tan_expr(   $3  ); }
             | TANH  '(' expr ')'   		{ $$ = context.create_tanh_expr(  $3  ); }
             | SQRT  '(' expr ')'   		{ $$ = context.create_sqrt_expr(  $3  ); }
;            
     
     num:      NUM                			{ $$ = $1;                             }
             | num '+' num        			{ $$ = $1 + $3;  }
             | num MINUSSIGN num      		{ $$ = $1 - $3;  }
             | num '*' num        			{ $$ = $1 * $3;  }
             | num '/' num        			{ $$ = $1 / $3;  }
             | MINUSSIGN num  %prec NEG     { $$ = -$2;        	       }
             | num '^' num        			{ $$ = pow( $1 , $3 );}
             | '(' num ')'        			{ $$ = $2;         	}
             | ACOS   '(' num ')'  			{ $$ = acos( $3 );   	}
             | ASIN   '(' num ')'  			{ $$ = asin( $3 );   	}
             | ATAN   '(' num ')'  			{ $$ = atan( $3 );   	}
             | ATAN   '(' num '\\' num ')'	{ $$ = atan2( $3, $5 );}
             | ABS    '(' num ')'   		{ $$ = fabs( $3 );   	}
             | EXP    '(' num ')'   		{ $$ = exp( $3 );    	}
             | LOG    '(' num ')'			{ $$ = log( $3 );     }
             | LOG10  '(' num ')'			{ $$ = log10( $3 );   }
             | SIN    '(' num ')'  			{ $$ = sin( $3 );     }
             | SINH   '(' num ')'  			{ $$ = sinh( $3 );    }
             | TAN    '(' num ')'  			{ $$ = tan( $3 );     }
             | TANH   '(' num ')'  			{ $$ = tanh( $3 );    }
             | SQRT   '(' num ')'  			{ $$ = sqrt( $3 );    }
             | COS    '(' num ')'   		{ $$ = cos( $3 );    }
;            

    unit:       UNIT                        { $$ = $1;         	                }
            |   unit '*' unit        	    { $$ = $1 * $3;    	                }
            |   unit '/' unit        		{ $$ = $1 / $3;    	                }
            |   unit '^' num       	        { assert(is_integer($3)); $$ = $1.pow( static_cast<int>( round($3) ) ); }
            |   '(' unit ')'                { $$ = $2;                          }
;
    quantity:   num unit                    { $$ = Quantity( $1*$2.getScaleFactor(), $2.getUnit() ); }
;            


%%

namespace SketcherExpressions{
	void ExpressionParser::error(const location_type& l, const std::string& m){ }
}
