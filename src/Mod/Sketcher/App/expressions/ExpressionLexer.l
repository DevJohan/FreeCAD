%top{
/* Lexer for the FreeCAD  Units language   */
/* (c) 2013 Juergen Riegel  LGPL           */


/* This disables inclusion of unistd.h, which is not available under Visual C++
 * on Win32. The C++ scanner uses STL streams instead. */
#define YY_NO_UNISTD_H
/* Prevent lexer header file from being included in lexer compilation context */
#ifndef yyHEADER_H
#define yyHEADER_H
#endif

}

/*** Flex Declarations and Options ***/

/* the manual says "somewhat more optimized" */
%option batch
%option never-interactive 
%option reentrant


/* no support for include files is planned */
%option noyywrap nounput 

/* Needed inclusions */
%{
	#include <sstream>
	#include "ExpressionParser.hpp"
	#include <Mod/Sketcher/App/expressions/ExpressionParserContext.h>
	using SketcherExpressions::ExpressionParser;
	using SketcherExpressions::ExpressionParserContext;
	
	#define YY_USER_ACTION set_location( location, yyleng );
	
%}


DIGIT    [0-9]
IDRULE   [a-zA-z][a-zA-Z_0-9]*
EXPO     [eE][-+]?[0-9]+

%x C_COMMENT


%% /*** Filter language Part ***/

%{
	// Code run each time yylex is called.
	reset_location(location);
%}

"\["            { BEGIN(C_COMMENT); }
<C_COMMENT>"\]" { BEGIN(INITIAL); }
<C_COMMENT>.    { ;}


[ \t]       				reset_location(location);
[\n]+       				reset_location(location);

[+()=/\\*^]  				{ return ExpressionParser::make_OPERATOR_CHAR( *yytext ,location ); }

"-"                         return ExpressionParser::make_MINUSSIGN( location );
"\xe2\x88\x92"              return ExpressionParser::make_MINUSSIGN( location );

"nm"           				return ExpressionParser::make_UNIT( ScaledUnit::NanoMetre        ,location ); // nano meter
"um"           				return ExpressionParser::make_UNIT( ScaledUnit::MicroMetre       ,location ); // micro meter
"\xC2\xB5m"    				return ExpressionParser::make_UNIT( ScaledUnit::MicroMetre       ,location ); // micro meter    (greek micro in UTF8)
"mm"           				return ExpressionParser::make_UNIT( ScaledUnit::MilliMetre       ,location ); // milli meter    (internal standard length)
"cm"           				return ExpressionParser::make_UNIT( ScaledUnit::CentiMetre       ,location ); // centi meter
"dm"           				return ExpressionParser::make_UNIT( ScaledUnit::DeciMetre        ,location ); // deci meter
"m"            				return ExpressionParser::make_UNIT( ScaledUnit::Metre            ,location ); // metre
"km"           				return ExpressionParser::make_UNIT( ScaledUnit::KiloMetre        ,location ); // kilo meter

"l"            				return ExpressionParser::make_UNIT( ScaledUnit::Liter            ,location ); // Liter      dm^3

"ug"           				return ExpressionParser::make_UNIT( ScaledUnit::MicroGram        ,location ); // micro gram
"\xC2\xB5g"    				return ExpressionParser::make_UNIT( ScaledUnit::MicroGram        ,location ); // micro gram
"mg"           				return ExpressionParser::make_UNIT( ScaledUnit::MilliGram        ,location ); // milli gram
"g"            				return ExpressionParser::make_UNIT( ScaledUnit::Gram             ,location ); // gram
"kg"           				return ExpressionParser::make_UNIT( ScaledUnit::KiloGram         ,location ); // kilo gram      (internal standard for mass) 
"t"            				return ExpressionParser::make_UNIT( ScaledUnit::Ton              ,location ); // Metric Tonne

"s"            				return ExpressionParser::make_UNIT( ScaledUnit::Second           ,location ); // second         (internal standard time)
"min"          				return ExpressionParser::make_UNIT( ScaledUnit::Minute           ,location ); // minute
"h"            				return ExpressionParser::make_UNIT( ScaledUnit::Hour             ,location ); // hour  

"A"            				return ExpressionParser::make_UNIT( ScaledUnit::Ampere           ,location ); // Ampere         (internal standard electric current)
"mA"           				return ExpressionParser::make_UNIT( ScaledUnit::MilliAmpere      ,location ); // milli Ampere         
"kA"           				return ExpressionParser::make_UNIT( ScaledUnit::KiloAmpere       ,location ); // kilo Ampere         
"MA"           				return ExpressionParser::make_UNIT( ScaledUnit::MegaAmpere       ,location ); // Mega Ampere         

"K"            				return ExpressionParser::make_UNIT( ScaledUnit::Kelvin           ,location ); // Kelvin         (internal standard thermodynamic temperature)
"mK"           				return ExpressionParser::make_UNIT( ScaledUnit::MilliKelvin      ,location ); // Kelvin         
"\xC2\xB5K"    				return ExpressionParser::make_UNIT( ScaledUnit::MicroKelvin      ,location ); // Kelvin         
"uK"           				return ExpressionParser::make_UNIT( ScaledUnit::MicroKelvin      ,location ); // Kelvin         

"mol"          				return ExpressionParser::make_UNIT( ScaledUnit::Mole             ,location ); // Mole           (internal standard amount of substance)        

"cd"           				return ExpressionParser::make_UNIT( ScaledUnit::Candela          ,location ); // Candela        (internal standard luminous intensity)        

"in"           				return ExpressionParser::make_UNIT( ScaledUnit::Inch             ,location ); // inch
"\""           				return ExpressionParser::make_UNIT( ScaledUnit::Inch             ,location ); // inch
"ft"           				return ExpressionParser::make_UNIT( ScaledUnit::Foot             ,location ); // foot
"'"            				return ExpressionParser::make_UNIT( ScaledUnit::Foot             ,location ); // foot
"thou"         				return ExpressionParser::make_UNIT( ScaledUnit::Thou             ,location ); // thou (in/1000)
"mil"          				return ExpressionParser::make_UNIT( ScaledUnit::Thou             ,location ); // mil  (the thou in US)
"yd"           				return ExpressionParser::make_UNIT( ScaledUnit::Yard             ,location ); // yard
"mi"           				return ExpressionParser::make_UNIT( ScaledUnit::Mile             ,location ); // mile



"lb"           				return ExpressionParser::make_UNIT( ScaledUnit::Pound            ,location ); // pound
"lbm"          				return ExpressionParser::make_UNIT( ScaledUnit::Pound            ,location ); // pound 
"oz"           				return ExpressionParser::make_UNIT( ScaledUnit::Ounce            ,location ); // ounce
"st"           				return ExpressionParser::make_UNIT( ScaledUnit::Stone            ,location ); // Stone
"cwt"          				return ExpressionParser::make_UNIT( ScaledUnit::Hundredweights   ,location ); // hundredweights

"lbf"          				return ExpressionParser::make_UNIT( ScaledUnit::PoundForce       ,location ); // pound

"N"            				return ExpressionParser::make_UNIT( ScaledUnit::Newton           ,location ); // Newton (kg*m/s^2)
"kN"           				return ExpressionParser::make_UNIT( ScaledUnit::KiloNewton       ,location ); // Newton 
"MN"           				return ExpressionParser::make_UNIT( ScaledUnit::MegaNewton       ,location ); // Newton 
"mN"           				return ExpressionParser::make_UNIT( ScaledUnit::MilliNewton      ,location ); // Newton 

"Pa"           				return ExpressionParser::make_UNIT( ScaledUnit::Pascal           ,location ); // Pascal (kg/m*s^2 or N/m^2) 
"kPa"          				return ExpressionParser::make_UNIT( ScaledUnit::KiloPascal       ,location ); // Pascal 
"MPa"          				return ExpressionParser::make_UNIT( ScaledUnit::MegaPascal       ,location ); // Pascal  
"GPa"          				return ExpressionParser::make_UNIT( ScaledUnit::GigaPascal       ,location ); // Pascal  

"Torr"         				return ExpressionParser::make_UNIT( ScaledUnit::Torr             ,location ); // portion of Pascal ( 101325/760 ) 
"mTorr"        				return ExpressionParser::make_UNIT( ScaledUnit::mTorr            ,location ); //  
"uTorr"        				return ExpressionParser::make_UNIT( ScaledUnit::yTorr            ,location ); //  
"\xC2\xB5Torr" 				return ExpressionParser::make_UNIT( ScaledUnit::yTorr            ,location ); //  

"psi"          				return ExpressionParser::make_UNIT( ScaledUnit::PSI              ,location ); // pounds/in^2 
"ksi"          				return ExpressionParser::make_UNIT( ScaledUnit::KSI              ,location ); // 1000 x pounds/in^2 

"W"            				return ExpressionParser::make_UNIT( ScaledUnit::Watt             ,location ); // Watt (kg*m^2/s^3) 
"VA"           				return ExpressionParser::make_UNIT( ScaledUnit::VoltAmpere       ,location ); // VoltAmpere (kg*m^2/s^3) 

"J"            				return ExpressionParser::make_UNIT( ScaledUnit::Joule            ,location ); // Joule (kg*m^2/s^2) 
"Nm"           				return ExpressionParser::make_UNIT( ScaledUnit::NewtonMeter      ,location ); // N*m = Joule 
"VAs"          				return ExpressionParser::make_UNIT( ScaledUnit::VoltAmpereSecond ,location ); // V*A*s = Joule 
"CV"           				return ExpressionParser::make_UNIT( ScaledUnit::WattSecond       ,location ); //  
"Ws"           				return ExpressionParser::make_UNIT( ScaledUnit::WattSecond       ,location ); // W*s = Joule 

"\xC2\xB0"     				return ExpressionParser::make_UNIT( ScaledUnit::Degree           ,location ); // degree         (internal standard angle)
"deg"          				return ExpressionParser::make_UNIT( ScaledUnit::Degree           ,location ); // degree         (internal standard angle)
"rad"          				return ExpressionParser::make_UNIT( ScaledUnit::Radian           ,location ); // radian         
"gon"          				return ExpressionParser::make_UNIT( ScaledUnit::Gon              ,location ); // gon         

{DIGIT}+"."{DIGIT}*{EXPO}? {  return ExpressionParser::make_NUM( atof(yytext) ,location);  }
"."{DIGIT}+{EXPO}?         {  return ExpressionParser::make_NUM( atof(yytext) ,location);  }


"pi"                   		{ return ExpressionParser::make_NUM( MATH_CONSTANT_PI , location); } // constant pi         
"e"                    		{ return ExpressionParser::make_NUM( MATH_CONSTANT_E  , location); } // constant e         
                     
"acos"                 		return ExpressionParser::make_ACOS(location);
"asin"                 		return ExpressionParser::make_ASIN(location);
"atan"                 		return ExpressionParser::make_ATAN(location);
"cos"                  		return ExpressionParser::make_COS(location);
"exp"                  		return ExpressionParser::make_EXP(location);
"abs"                  		return ExpressionParser::make_ABS(location);
"log"                  		return ExpressionParser::make_LOG(location);
"log10"                		return ExpressionParser::make_LOG10(location);
"sin"                  		return ExpressionParser::make_SIN(location);
"sinh"                 		return ExpressionParser::make_SINH(location);
"tan"                  		return ExpressionParser::make_TAN(location);
"tanh"                 		return ExpressionParser::make_TANH(location);
"sqrt"                 		return ExpressionParser::make_SQRT(location);

{IDRULE}					return ExpressionParser::make_IDENTIFIER( yytext ,location );

.                      		{
								std::ostringstream oss; oss << "Unknown symbol "<< yytext;
								throw ExpressionParser::syntax_error( location, oss.str() );
							}
<<EOF>>    					return ExpressionParser::make_END(location);

%%