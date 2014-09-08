/***************************************************************************
 *   Copyright (c) 2013 Juergen Riegel                                     *
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

#include "ScaledUnit.h"
#include <cmath>
#include <Base/Exception.h>
#include <Base/UnitsApi.h>

// suppress annoying warnings from generated source files
#ifdef _MSC_VER
# pragma warning(disable : 4003)
# pragma warning(disable : 4018)
# pragma warning(disable : 4065)
# pragma warning( disable : 4273 )
# pragma warning(disable : 4335) // disable MAC file format warning on VC
#endif

namespace SketcherExpressions{


ScaledUnit::ScaledUnit()
{
    this->_ScaleFactor = 1.0;
}

ScaledUnit::ScaledUnit(const ScaledUnit& that)
{
    *this = that ;
}

ScaledUnit::ScaledUnit(double Scale, const Unit& unit)
{
    this->_Unit = unit;
    this->_ScaleFactor = Scale;
}


double ScaledUnit::getScaleFactorAs(const ScaledUnit &q) const {

    return _ScaleFactor/q.getScaleFactor();
}



bool ScaledUnit::operator ==(const ScaledUnit& that) const
{
    return (this->_ScaleFactor == that._ScaleFactor) && (this->_Unit == that._Unit) ;
}

ScaledUnit ScaledUnit::operator *(const ScaledUnit &p) const
{
    return ScaledUnit(this->_ScaleFactor * p._ScaleFactor,this->_Unit * p._Unit);
}

ScaledUnit ScaledUnit::operator /(const ScaledUnit &p) const
{
    return ScaledUnit(this->_ScaleFactor / p._ScaleFactor,this->_Unit / p._Unit);
}

ScaledUnit ScaledUnit::pow(int power) const
{
    return ScaledUnit(
        std::pow(_ScaleFactor, power),
        this->_Unit.pow(power) );
}

ScaledUnit& ScaledUnit::operator = (const ScaledUnit &New)
{
    this->_ScaleFactor = New._ScaleFactor;
    this->_Unit = New._Unit;
    return *this;
}

/// true if it has a number without a unit
bool ScaledUnit::isDimensionless(void)const
{
    return isValid() && _Unit.isEmpty();
}
// true if it has a number and a valid unit
bool ScaledUnit::isScaledUnit(void)const
{
    return isValid() && !_Unit.isEmpty();
}
void ScaledUnit::setInvalid(void)
{
    _ScaleFactor = -1 ;
}

// === Predefined types =====================================================

ScaledUnit ScaledUnit::NanoMetre        (1.0e-6         ,Unit(1));
ScaledUnit ScaledUnit::MicroMetre       (1.0e-3         ,Unit(1));
ScaledUnit ScaledUnit::MilliMetre       (1.0            ,Unit(1));
ScaledUnit ScaledUnit::CentiMetre       (10.0           ,Unit(1));
ScaledUnit ScaledUnit::DeciMetre        (100.0          ,Unit(1));
ScaledUnit ScaledUnit::Metre            (1.0e3          ,Unit(1));
ScaledUnit ScaledUnit::KiloMetre        (1.0e6          ,Unit(1));

ScaledUnit ScaledUnit::Liter            (1.0e6          ,Unit(3));

ScaledUnit ScaledUnit::MicroGram        (1.0e-9         ,Unit(0,1));
ScaledUnit ScaledUnit::MilliGram        (1.0e-6         ,Unit(0,1));
ScaledUnit ScaledUnit::Gram             (1.0e-3         ,Unit(0,1));
ScaledUnit ScaledUnit::KiloGram         (1.0            ,Unit(0,1));
ScaledUnit ScaledUnit::Ton              (1.0e3          ,Unit(0,1));

ScaledUnit ScaledUnit::Second           (1.0            ,Unit(0,0,1));
ScaledUnit ScaledUnit::Minute           (60.0           ,Unit(0,0,1));
ScaledUnit ScaledUnit::Hour             (3600.0         ,Unit(0,0,1));

ScaledUnit ScaledUnit::Ampere           (1.0           ,Unit(0,0,0,1));
ScaledUnit ScaledUnit::MilliAmpere      (0.001         ,Unit(0,0,0,1));
ScaledUnit ScaledUnit::KiloAmpere       (1000.0        ,Unit(0,0,0,1));
ScaledUnit ScaledUnit::MegaAmpere       (1.0e6         ,Unit(0,0,0,1));

ScaledUnit ScaledUnit::Kelvin           (1.0           ,Unit(0,0,0,0,1));
ScaledUnit ScaledUnit::MilliKelvin      (0.001         ,Unit(0,0,0,0,1));
ScaledUnit ScaledUnit::MicroKelvin      (0.000001      ,Unit(0,0,0,0,1));

ScaledUnit ScaledUnit::Mole             (1.0           ,Unit(0,0,0,0,0,1));

ScaledUnit ScaledUnit::Candela          (1.0           ,Unit(0,0,0,0,0,0,1));

ScaledUnit ScaledUnit::Inch             (25.4          ,Unit(1));
ScaledUnit ScaledUnit::Foot             (304.8         ,Unit(1));
ScaledUnit ScaledUnit::Thou             (0.0254        ,Unit(1));
ScaledUnit ScaledUnit::Yard             (914.4         ,Unit(1));
ScaledUnit ScaledUnit::Mile             (1609344.0     ,Unit(1));

ScaledUnit ScaledUnit::Pound            (0.45359237    ,Unit(0,1));
ScaledUnit ScaledUnit::Ounce            (0.0283495231  ,Unit(0,1));
ScaledUnit ScaledUnit::Stone            (6.35029318    ,Unit(0,1));
ScaledUnit ScaledUnit::Hundredweights   (50.80234544   ,Unit(0,1));

ScaledUnit ScaledUnit::PoundForce       (224.81        ,Unit(1,1,-2));  // Newton  are ~= 0.22481 lbF

ScaledUnit ScaledUnit::Newton           (1000.0        ,Unit(1,1,-2));  // Newton (kg*m/s^2)
ScaledUnit ScaledUnit::KiloNewton       (1e+6          ,Unit(1,1,-2));
ScaledUnit ScaledUnit::MegaNewton       (1e+9          ,Unit(1,1,-2));
ScaledUnit ScaledUnit::MilliNewton      (1.0           ,Unit(1,1,-2));

ScaledUnit ScaledUnit::Pascal           (0.001         ,Unit(-1,1,-2)); // Pascal (kg/m*s^2 or N/m^2)
ScaledUnit ScaledUnit::KiloPascal       (1.00          ,Unit(-1,1,-2));
ScaledUnit ScaledUnit::MegaPascal       (1000.0        ,Unit(-1,1,-2));
ScaledUnit ScaledUnit::GigaPascal       (1e+6          ,Unit(-1,1,-2));

ScaledUnit ScaledUnit::Torr             (101.325/760.0 ,Unit(-1,1,-2)); // Torr is a defined fraction of Pascal (kg/m*s^2 or N/m^2)
ScaledUnit ScaledUnit::mTorr            (0.101325/760.0,Unit(-1,1,-2)); // Torr is a defined fraction of Pascal (kg/m*s^2 or N/m^2)
ScaledUnit ScaledUnit::yTorr            (0.000101325/760.0 ,Unit(-1,1,-2)); // Torr is a defined fraction of Pascal (kg/m*s^2 or N/m^2)

ScaledUnit ScaledUnit::PSI              (0.145038      ,Unit(-1,1,-2)); // pounds/in^2
ScaledUnit ScaledUnit::KSI              (145.038       ,Unit(-1,1,-2)); // 1000 x pounds/in^2

ScaledUnit ScaledUnit::Watt             (1e+6          ,Unit(2,1,-3));  // Watt (kg*m^2/s^3)
ScaledUnit ScaledUnit::VoltAmpere       (1e+6          ,Unit(2,1,-3));  // VoltAmpere (kg*m^2/s^3)

ScaledUnit ScaledUnit::Joule            (1e+6          ,Unit(2,1,-2));  // Joule (kg*m^2/s^2)
ScaledUnit ScaledUnit::NewtonMeter      (1e+6          ,Unit(2,1,-2));  // Joule (kg*m^2/s^2)
ScaledUnit ScaledUnit::VoltAmpereSecond (1e+6          ,Unit(2,1,-2));  // Joule (kg*m^2/s^2)
ScaledUnit ScaledUnit::WattSecond       (1e+6          ,Unit(2,1,-2));  // Joule (kg*m^2/s^2)

ScaledUnit ScaledUnit::KMH              (277.778       ,Unit(1,0,-1));  // km/h
ScaledUnit ScaledUnit::MPH              (447.04        ,Unit(1,0,-1));  // Mile/h

ScaledUnit ScaledUnit::Degree           (1.0           ,Unit(0,0,0,0,0,0,0,1)); // degree         (internal standard angle)
ScaledUnit ScaledUnit::Radian           (180/M_PI      ,Unit(0,0,0,0,0,0,0,1)); // radian
ScaledUnit ScaledUnit::Gon              (360.0/400.0   ,Unit(0,0,0,0,0,0,0,1)); // gon



} // namespace SketcherExpressions

