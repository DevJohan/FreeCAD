/* Parser for the FreeCAD  Units language           */
/* (c) 2013 Juergen Riegel  LGPL                    */

/* Bison options.  */
%defines
%language "C++"
%define api.namespace {QuantityParserExp}
%define api.value.type variant
%define api.token.constructor
%define parser_class_name {QuantityParser}

%code requires {
	#include "PreCompiled.h"
	#ifndef _PreComp_
	# include <sstream>
	#endif

    #define BaseExport
    #include <Base/Unit.h>
    #include <Base/Quantity.h>
    #include <Mod/Sketcher/App/expressions/ScaledUnit.h>
    using Base::Unit;
    using Base::Quantity;
    using SketcherExpressions::ScaledUnit;

	namespace QuantityParserExp {
		class QuantityParserContext;
	}		
}

%param { QuantityParserContext& context }

/* Represents the many different ways we can access our data */
%{
//        #define YYSTYPE Quantity
//        #define yyparse Quantity_yyparse
//        #define yyerror Quantity_yyerror
        #ifndef  DOUBLE_MAX
        # define DOUBLE_MAX 1.7976931348623157E+308    /* max decimal value of a "double"*/
        #endif
        #ifndef  DOUBLE_MIN
        # define DOUBLE_MIN 2.2250738585072014E-308    /* min decimal value of a "double"*/
        #endif
%}

%code {
	#include <limits>
	template <typename T>
    bool is_integer(T value){
    	const T result = abs(std::fmod(value,1.0)); 
    	return ( result > .5 ? 1.0 - result : result ) < std::numeric_limits<T>::epsilon();
    }
	 
	#include <Mod/Sketcher/App/expressions/QuantityParserContextExp.h>
	namespace QuantityParserExp {
		YY_DECL;
	}
}

/* Bison declarations.  */
%token <ScaledUnit> UNIT "unit";
%token <double> NUM "floating point number";
%token <std::string> IDENTIFIER "identifier";
%token <char> OPERATOR "operator char";
%token MINUSSIGN;
%token ACOS ASIN ATAN ATAN2 COS EXP ABS MOD LOG LOG10 POW SIN SINH TAN TANH SQRT;
%left MINUSSIGN '+'
%left '*' '/'
%left NEG     /* negation--unary minus */
%right '^'    /* exponentiation */

%{
        #include <cmath>
               
%}

%type <Quantity> expr
%type <double> num
%type <ScaledUnit> unit
%type <Quantity> quantity

%start input

%%

    input:                                  { context.result = Quantity(DOUBLE_MIN); /* empty input */ }
            |  expr                         { context.result = $1;            }
            |  num                          { context.result = $1;            }
            |  unit                         { context.result = Quantity(0,$1.getUnit());     }
            |  quantity                     { context.result = $1     ;            }
            |  quantity quantity            { context.result = $1 + $2;            }
 ;   
     expr:    num                			{ $$ = $1;                             }
             | expr '+' expr        		{ $$ = $1.getValue() + $3.getValue();  }
             | expr MINUSSIGN expr          { $$ = $1.getValue() - $3.getValue();  }
             | expr '*' expr       			{ $$ = $1.getValue() * $3.getValue();  }
             | expr '/' expr       			{ $$ = $1.getValue() / $3.getValue();  }
             | MINUSSIGN expr  %prec NEG    { $$ = -$2.getValue();        	       }
             | expr '^' expr       			{ $$ = pow ($1.getValue(), $3.getValue());}
             | '(' expr ')'        			{ $$ = $2;         	}
             | ACOS  '(' expr ')'  			{ $$ = acos($3.getValue());   	}
             | ASIN  '(' expr ')'  			{ $$ = asin($3.getValue());   	}
             | ATAN  '(' expr ')'  			{ $$ = atan($3.getValue());   	}
             | ATAN  '(' expr '\\' expr ')'	{ $$ = atan2($3.getValue(),$5.getValue());}
             | ABS  '(' expr ')'   			{ $$ = fabs($3.getValue());   	}
             | EXP  '(' expr ')'   			{ $$ = exp($3.getValue());    	}
             | LOG  '(' expr ')'			{ $$ = log($3.getValue());     }
             | LOG10  '(' expr ')'			{ $$ = log10($3.getValue());   }
             | SIN  '(' expr ')'   			{ $$ = sin($3.getValue());     }
             | SINH '(' expr ')'   			{ $$ = sinh($3.getValue());    }
             | TAN  '(' expr ')'   			{ $$ = tan($3.getValue());     }
             | TANH  '(' expr ')'   		{ $$ = tanh($3.getValue());    }
             | SQRT  '(' expr ')'   		{ $$ = sqrt($3.getValue());    }
             | COS  '(' expr ')'   			{ $$ = cos($3.getValue());    }
;            
     
     num:      NUM                			{ $$ = $1;                             }
             | num '+' num        			{ $$ = $1 + $3;  }
             | num MINUSSIGN num      		{ $$ = $1 - $3;  }
             | num '*' num        			{ $$ = $1 * $3;  }
             | num '/' num        			{ $$ = $1 / $3;  }
             | MINUSSIGN num  %prec NEG     { $$ = -$2;        	       }
             | num '^' num        			{ $$ = pow( $1 , $3 );}
             | '(' num ')'        			{ $$ = $2;         	}
             | ACOS  '(' num ')'  			{ $$ = acos( $3 );   	}
             | ASIN  '(' num ')'  			{ $$ = asin( $3 );   	}
             | ATAN  '(' num ')'  			{ $$ = atan( $3 );   	}
             | ATAN  '(' num '\\' num ')'	{ $$ = atan2( $3, $5 );}
             | ABS  '(' num ')'   			{ $$ = fabs( $3 );   	}
             | EXP  '(' num ')'   			{ $$ = exp( $3 );    	}
             | LOG  '(' num ')'				{ $$ = log( $3 );     }
             | LOG10  '(' num ')'			{ $$ = log10( $3 );   }
             | SIN  '(' num ')'   			{ $$ = sin( $3 );     }
             | SINH '(' num ')'   			{ $$ = sinh( $3 );    }
             | TAN  '(' num ')'   			{ $$ = tan( $3 );     }
             | TANH  '(' num ')'   			{ $$ = tanh( $3 );    }
             | SQRT  '(' num ')'   			{ $$ = sqrt( $3 );    }
             | COS  '(' num ')'   			{ $$ = cos( $3 );    }
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
