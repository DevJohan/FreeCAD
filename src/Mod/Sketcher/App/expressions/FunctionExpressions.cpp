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

#include "FunctionExpressions.h"

namespace SketcherExpressions{

template <int arg_count, ExpressionType::Types type, ExpressionBaseTypes data_value_type>
FunctionExpression<arg_count,type,data_value_type>::FunctionExpression(){ }

template <ExpressionType::Types type, ExpressionBaseTypes data_value_type>
FunctionExpression<dynamic,type,data_value_type>::FunctionExpression(){ }

template <int arg_count, ExpressionType::Types type, ExpressionBaseTypes data_value_type>
FunctionExpression<arg_count,type,data_value_type>::~FunctionExpression(){ }

template <ExpressionType::Types type, ExpressionBaseTypes data_value_type>
FunctionExpression<dynamic,type,data_value_type>::~FunctionExpression(){ }


//SIN
template < ExpressionBaseTypes data_value_type>
SinExpression<data_value_type>::SinExpression(ExpressionReference argument){

}

template < ExpressionBaseTypes data_value_type>
SinExpression<data_value_type>::~SinExpression(){

}

//COS
template < ExpressionBaseTypes data_value_type>
CosExpression<data_value_type>::CosExpression(ExpressionReference argument){

}

template < ExpressionBaseTypes data_value_type>
CosExpression<data_value_type>::~CosExpression(){

}

//TAN
template < ExpressionBaseTypes data_value_type>
TanExpression<data_value_type>::TanExpression(ExpressionReference argument){

}

template < ExpressionBaseTypes data_value_type>
TanExpression<data_value_type>::~TanExpression(){

}

//ASIN
template < ExpressionBaseTypes data_value_type>
ArcsinExpression<data_value_type>::ArcsinExpression(ExpressionReference argument){

}

template < ExpressionBaseTypes data_value_type>
ArcsinExpression<data_value_type>::~ArcsinExpression(){

}

//ACOS
template < ExpressionBaseTypes data_value_type>
ArccosExpression<data_value_type>::ArccosExpression(ExpressionReference argument){
}
template < ExpressionBaseTypes data_value_type>
ArccosExpression<data_value_type>::~ArccosExpression(){

}

//ATAN
template < ExpressionBaseTypes data_value_type>
ArctanExpression<data_value_type>::ArctanExpression(ExpressionReference argument){

}
template < ExpressionBaseTypes data_value_type>
ArctanExpression<data_value_type>::~ArctanExpression(){

}

//ATAN2
template < ExpressionBaseTypes data_value_type>
Arctan2Expression<data_value_type>::Arctan2Expression(ExpressionReference y, ExpressionReference x){

}
template < ExpressionBaseTypes data_value_type>
Arctan2Expression<data_value_type>::~Arctan2Expression(){
}

//ABS
template < ExpressionBaseTypes data_value_type>
AbsExpression<data_value_type>::AbsExpression(ExpressionReference argument){
}
template < ExpressionBaseTypes data_value_type>
AbsExpression<data_value_type>::~AbsExpression(){
}

//POW
template < ExpressionBaseTypes data_value_type>
PowExpression<data_value_type>::PowExpression(ExpressionReference base, ExpressionReference argument){

}
template < ExpressionBaseTypes data_value_type>
PowExpression<data_value_type>::~PowExpression(){
}

//LOG
template < ExpressionBaseTypes data_value_type>
LogExpression<data_value_type>::LogExpression(ExpressionReference base, ExpressionReference argument){
}
template < ExpressionBaseTypes data_value_type>
LogExpression<data_value_type>::~LogExpression(){

}

//EXP
template < ExpressionBaseTypes data_value_type>
ExpExpression<data_value_type>::ExpExpression(ExpressionReference argument){

}
template < ExpressionBaseTypes data_value_type>
ExpExpression<data_value_type>::~ExpExpression(){
}

//LN
template < ExpressionBaseTypes data_value_type>
LnExpression<data_value_type>::LnExpression(ExpressionReference argument){
}

template < ExpressionBaseTypes data_value_type>
LnExpression<data_value_type>::~LnExpression(){
}


//SINH
template < ExpressionBaseTypes data_value_type>
SinhExpression<data_value_type>::SinhExpression(ExpressionReference argument){
}

template < ExpressionBaseTypes data_value_type>
SinhExpression<data_value_type>::~SinhExpression(){
}

//COSH
template < ExpressionBaseTypes data_value_type>
CoshExpression<data_value_type>::CoshExpression(ExpressionReference argument){

}

template < ExpressionBaseTypes data_value_type>
CoshExpression<data_value_type>::~CoshExpression(){

}

//TANH
template < ExpressionBaseTypes data_value_type>
TanhExpression<data_value_type>::TanhExpression(ExpressionReference argument){

}

template < ExpressionBaseTypes data_value_type>
TanhExpression<data_value_type>::~TanhExpression(){

}

//SQRT
template < ExpressionBaseTypes data_value_type>
SqrtExpression<data_value_type>::SqrtExpression(ExpressionReference argument){

}
template < ExpressionBaseTypes data_value_type>
SqrtExpression<data_value_type>::~SqrtExpression(){

}



//SIN
template class SinExpression<plain>;
template class SinExpression<quantity>;

//COS
template class CosExpression<plain>;
template class CosExpression<quantity>;

//TAN
template class TanExpression<plain>;
template class TanExpression<quantity>;

//ASIN
template class ArcsinExpression<plain>;
template class ArcsinExpression<quantity>;

//ACOS
template class ArccosExpression<plain>;
template class ArccosExpression<quantity>;

//ATAN
template class ArctanExpression<plain>;
template class ArctanExpression<quantity>;

//ATAN2
template class Arctan2Expression<plain>;
template class Arctan2Expression<quantity>;

//ABS
template class AbsExpression<plain>;
template class AbsExpression<quantity>;

//POW
template class PowExpression<plain>;
template class PowExpression<quantity>;

//LOG
template class LogExpression<plain>;
template class LogExpression<quantity>;

//EXP
template class ExpExpression<plain>;
template class ExpExpression<quantity>;

//LN
template class LnExpression<plain>;
template class LnExpression<quantity>;

//SINH
template class SinhExpression<plain>;
template class SinhExpression<quantity>;

//COSH
template class CoshExpression<plain>;
template class CoshExpression<quantity>;

//TANH
template class TanhExpression<plain>;
template class TanhExpression<quantity>;

//SQRT
template class SqrtExpression<plain>;
template class SqrtExpression<quantity>;

} // namespace SketcherExpressions
