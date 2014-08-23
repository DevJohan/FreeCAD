%{
/* Lexer for the FreeCAD  Units language   */
/* (c) 2013 Juergen Riegel  LGPL           */


/* This disables inclusion of unistd.h, which is not available under Visual C++
 * on Win32. The C++ scanner uses STL streams instead. */
#define YY_NO_UNISTD_H
	#include "QuantityParserExp.hpp"
	#include <Mod/Sketcher/App/expressions/QuantityParserContextExp.h>
	using QuantityParserExp::QuantityParser;
	using QuantityParserExp::QuantityParserContext;
	using namespace QuantityParserExp;
%}

/*** Flex Declarations and Options ***/


/* the manual says "somewhat more optimized" */
/*%option c++*/
%option batch
%option never-interactive 


/* no support for include files is planned */
%option noyywrap nounput 

DIGIT    [0-9]
IDRULE   [a-zA-z][a-zA-Z_0-9]*
EXPO     [eE][-+]?[0-9]+

%x C_COMMENT

%% /*** Filter language Part ***/

"\["            { BEGIN(C_COMMENT); }
<C_COMMENT>"\]" { BEGIN(INITIAL); }
<C_COMMENT>.    { ;}


[ \t]       				;
[\n]+       				;

[+()=/\\*^]  				{ return QuantityParser::make_OPERATOR( *yytext ); }

"-"                         return QuantityParser::make_MINUSSIGN();
"\xe2\x88\x92"              return QuantityParser::make_MINUSSIGN();

"nm"           				return QuantityParser::make_UNIT( ScaledUnit::NanoMetre         ); // nano meter
"um"           				return QuantityParser::make_UNIT( ScaledUnit::MicroMetre        ); // micro meter
"\xC2\xB5m"    				return QuantityParser::make_UNIT( ScaledUnit::MicroMetre        ); // micro meter    (greek micro in UTF8)
"mm"           				return QuantityParser::make_UNIT( ScaledUnit::MilliMetre        ); // milli meter    (internal standard length)
"cm"           				return QuantityParser::make_UNIT( ScaledUnit::CentiMetre        ); // centi meter
"dm"           				return QuantityParser::make_UNIT( ScaledUnit::DeciMetre         ); // deci meter
"m"            				return QuantityParser::make_UNIT( ScaledUnit::Metre             ); // metre
"km"           				return QuantityParser::make_UNIT( ScaledUnit::KiloMetre         ); // kilo meter

"l"            				return QuantityParser::make_UNIT( ScaledUnit::Liter             ); // Liter      dm^3

"ug"           				return QuantityParser::make_UNIT( ScaledUnit::MicroGram         ); // micro gram
"\xC2\xB5g"    				return QuantityParser::make_UNIT( ScaledUnit::MicroGram         ); // micro gram
"mg"           				return QuantityParser::make_UNIT( ScaledUnit::MilliGram         ); // milli gram
"g"            				return QuantityParser::make_UNIT( ScaledUnit::Gram              ); // gram
"kg"           				return QuantityParser::make_UNIT( ScaledUnit::KiloGram          ); // kilo gram      (internal standard for mass) 
"t"            				return QuantityParser::make_UNIT( ScaledUnit::Ton               ); // Metric Tonne

"s"            				return QuantityParser::make_UNIT( ScaledUnit::Second            ); // second         (internal standard time)
"min"          				return QuantityParser::make_UNIT( ScaledUnit::Minute            ); // minute
"h"            				return QuantityParser::make_UNIT( ScaledUnit::Hour              ); // hour  

"A"            				return QuantityParser::make_UNIT( ScaledUnit::Ampere            ); // Ampere         (internal standard electric current)
"mA"           				return QuantityParser::make_UNIT( ScaledUnit::MilliAmpere       ); // milli Ampere         
"kA"           				return QuantityParser::make_UNIT( ScaledUnit::KiloAmpere        ); // kilo Ampere         
"MA"           				return QuantityParser::make_UNIT( ScaledUnit::MegaAmpere        ); // Mega Ampere         

"K"            				return QuantityParser::make_UNIT( ScaledUnit::Kelvin            ); // Kelvin         (internal standard thermodynamic temperature)
"mK"           				return QuantityParser::make_UNIT( ScaledUnit::MilliKelvin       ); // Kelvin         
"\xC2\xB5K"    				return QuantityParser::make_UNIT( ScaledUnit::MicroKelvin       ); // Kelvin         
"uK"           				return QuantityParser::make_UNIT( ScaledUnit::MicroKelvin       ); // Kelvin         

"mol"          				return QuantityParser::make_UNIT( ScaledUnit::Mole              ); // Mole           (internal standard amount of substance)        

"cd"           				return QuantityParser::make_UNIT( ScaledUnit::Candela           ); // Candela        (internal standard luminous intensity)        

"in"           				return QuantityParser::make_UNIT( ScaledUnit::Inch              ); // inch
"\""           				return QuantityParser::make_UNIT( ScaledUnit::Inch              ); // inch
"ft"           				return QuantityParser::make_UNIT( ScaledUnit::Foot              ); // foot
"'"            				return QuantityParser::make_UNIT( ScaledUnit::Foot              ); // foot
"thou"         				return QuantityParser::make_UNIT( ScaledUnit::Thou              ); // thou (in/1000)
"mil"          				return QuantityParser::make_UNIT( ScaledUnit::Thou              ); // mil  (the thou in US)
"yd"           				return QuantityParser::make_UNIT( ScaledUnit::Yard              ); // yard
"mi"           				return QuantityParser::make_UNIT( ScaledUnit::Mile              ); // mile



"lb"           				return QuantityParser::make_UNIT( ScaledUnit::Pound             ); // pound
"lbm"          				return QuantityParser::make_UNIT( ScaledUnit::Pound             ); // pound 
"oz"           				return QuantityParser::make_UNIT( ScaledUnit::Ounce             ); // ounce
"st"           				return QuantityParser::make_UNIT( ScaledUnit::Stone             ); // Stone
"cwt"          				return QuantityParser::make_UNIT( ScaledUnit::Hundredweights    ); // hundredweights

"lbf"          				return QuantityParser::make_UNIT( ScaledUnit::PoundForce        ); // pound

"N"            				return QuantityParser::make_UNIT( ScaledUnit::Newton            ); // Newton (kg*m/s^2)
"kN"           				return QuantityParser::make_UNIT( ScaledUnit::KiloNewton        ); // Newton 
"MN"           				return QuantityParser::make_UNIT( ScaledUnit::MegaNewton        ); // Newton 
"mN"           				return QuantityParser::make_UNIT( ScaledUnit::MilliNewton       ); // Newton 

"Pa"           				return QuantityParser::make_UNIT( ScaledUnit::Pascal            ); // Pascal (kg/m*s^2 or N/m^2) 
"kPa"          				return QuantityParser::make_UNIT( ScaledUnit::KiloPascal        ); // Pascal 
"MPa"          				return QuantityParser::make_UNIT( ScaledUnit::MegaPascal        ); // Pascal  
"GPa"          				return QuantityParser::make_UNIT( ScaledUnit::GigaPascal        ); // Pascal  

"Torr"         				return QuantityParser::make_UNIT( ScaledUnit::Torr              ); // portion of Pascal ( 101325/760 ) 
"mTorr"        				return QuantityParser::make_UNIT( ScaledUnit::mTorr             ); //  
"uTorr"        				return QuantityParser::make_UNIT( ScaledUnit::yTorr             ); //  
"\xC2\xB5Torr" 				return QuantityParser::make_UNIT( ScaledUnit::yTorr             ); //  

"psi"          				return QuantityParser::make_UNIT( ScaledUnit::PSI               ); // pounds/in^2 
"ksi"          				return QuantityParser::make_UNIT( ScaledUnit::KSI               ); // 1000 x pounds/in^2 

"W"            				return QuantityParser::make_UNIT( ScaledUnit::Watt              ); // Watt (kg*m^2/s^3) 
"VA"           				return QuantityParser::make_UNIT( ScaledUnit::VoltAmpere        ); // VoltAmpere (kg*m^2/s^3) 

"J"            				return QuantityParser::make_UNIT( ScaledUnit::Joule             ); // Joule (kg*m^2/s^2) 
"Nm"           				return QuantityParser::make_UNIT( ScaledUnit::NewtonMeter       ); // N*m = Joule 
"VAs"          				return QuantityParser::make_UNIT( ScaledUnit::VoltAmpereSecond  ); // V*A*s = Joule 
"CV"           				return QuantityParser::make_UNIT( ScaledUnit::WattSecond        ); //  
"Ws"           				return QuantityParser::make_UNIT( ScaledUnit::WattSecond        ); // W*s = Joule 

"\xC2\xB0"     				return QuantityParser::make_UNIT( ScaledUnit::Degree            ); // degree         (internal standard angle)
"deg"          				return QuantityParser::make_UNIT( ScaledUnit::Degree            ); // degree         (internal standard angle)
"rad"          				return QuantityParser::make_UNIT( ScaledUnit::Radian            ); // radian         
"gon"          				return QuantityParser::make_UNIT( ScaledUnit::Gon               ); // gon         

{DIGIT}+"."{DIGIT}*{EXPO}? {  return QuantityParser::make_NUM( atof(yytext) );  }
"."{DIGIT}+{EXPO}?         {  return QuantityParser::make_NUM( atof(yytext) );  }


"pi"                   		{ return QuantityParser::make_NUM( M_PI ); } // constant pi         
"e"                    		{ return QuantityParser::make_NUM( M_E ); } // constant e         
                     
"acos"                 		return QuantityParser::make_ACOS();
"asin"                 		return QuantityParser::make_ASIN();
"atan"                 		return QuantityParser::make_ATAN();
"cos"                  		return QuantityParser::make_COS();
"exp"                  		return QuantityParser::make_EXP();
"abs"                  		return QuantityParser::make_ABS();
"mod"                  		return QuantityParser::make_MOD();
"log"                  		return QuantityParser::make_LOG();
"log10"                		return QuantityParser::make_LOG10();
"pow"                  		return QuantityParser::make_POW();
"sin"                  		return QuantityParser::make_SIN();
"sinh"                 		return QuantityParser::make_SINH();
"tan"                  		return QuantityParser::make_TAN();
"tanh"                 		return QuantityParser::make_TANH();
"sqrt"                 		return QuantityParser::make_SQRT();

{IDRULE}					return QuantityParser::make_IDENTIFIER(yytext);

.                      		return *yytext;

%%