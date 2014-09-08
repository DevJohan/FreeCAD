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


#ifndef SKETCHER_SCALEDUNIT_H_
#define SKETCHER_SCALEDUNIT_H_



#include <Base/Unit.h>

#ifndef  DOUBLE_MAX
# define DOUBLE_MAX 1.7976931348623157E+308    /* max decimal value of a "double"*/
#endif
#ifndef  DOUBLE_MIN
# define DOUBLE_MIN 2.2250738585072014E-308    /* min decimal value of a "double"*/
#endif

namespace SketcherExpressions {
using Base::Unit;

/**
 * The ScaledUnit class.
 */
class SketcherExport ScaledUnit {
public:
    /// default constructor
    ScaledUnit(void);
    ScaledUnit(const ScaledUnit&);
    ScaledUnit(double Value, const Unit& unit=Unit());
    /// Destruction
    ~ScaledUnit () {}

    /** Operators. */
    //@{
    ScaledUnit operator *(const ScaledUnit &p) const;
    ScaledUnit operator /(const ScaledUnit &p) const;
    bool operator ==(const ScaledUnit&) const;
    ScaledUnit& operator =(const ScaledUnit&);
    ScaledUnit pow(int) const;
    //@}

    /// returns the unit of the quantity
	const Unit & getUnit(void) const{return _Unit;}
    /// set the unit of the quantity
    void setUnit(const Unit &un){_Unit = un;}
    /// get the Value of the quantity
	double getScaleFactor(void) const{return _ScaleFactor;}
    /// set the value of the quantity
    void setScaleFactor(double val){_ScaleFactor = val;}
    /** get the ScaleFactor in a special unit given as ScaledUnit.
      * One can use one of the predefined ScaledUnits in this class
      */
    double getScaleFactorAs(const ScaledUnit &)const;


    /// true if it has a positive scale factor without a unit
    bool isDimensionless(void)const;
    /// true if it has a positive scale factor and a valid unit
    bool isScaledUnit(void)const;
    /// true if it has a positive scale factor with or without a unit
    bool isValid(void)const{return _ScaleFactor > 0.0;}
    /// sets the quantity invalid
    void setInvalid(void);


    /** Predefined Unit types. */
    //@{
	static ScaledUnit NanoMetre;
	static ScaledUnit MicroMetre;
	static ScaledUnit CentiMetre;
	static ScaledUnit DeciMetre;
	static ScaledUnit Metre;
	static ScaledUnit MilliMetre;
	static ScaledUnit KiloMetre;

	static ScaledUnit Liter;

	static ScaledUnit MicroGram;
	static ScaledUnit MilliGram;
	static ScaledUnit Gram;
	static ScaledUnit KiloGram;
	static ScaledUnit Ton;

	static ScaledUnit Second;
	static ScaledUnit Minute;
	static ScaledUnit Hour;

	static ScaledUnit Ampere;
	static ScaledUnit MilliAmpere;
	static ScaledUnit KiloAmpere;
	static ScaledUnit MegaAmpere;

	static ScaledUnit Kelvin;
	static ScaledUnit MilliKelvin;
	static ScaledUnit MicroKelvin;

	static ScaledUnit Mole;

	static ScaledUnit Candela;

	static ScaledUnit Inch;
	static ScaledUnit Foot;
	static ScaledUnit Thou;
	static ScaledUnit Yard;

	static ScaledUnit Pound;
	static ScaledUnit Ounce;
	static ScaledUnit Stone;
	static ScaledUnit Hundredweights;
	static ScaledUnit Mile;

    static ScaledUnit PoundForce;

	static ScaledUnit Newton;
	static ScaledUnit KiloNewton;
	static ScaledUnit MegaNewton;
	static ScaledUnit MilliNewton;

	static ScaledUnit Pascal;
	static ScaledUnit KiloPascal;
	static ScaledUnit MegaPascal;
	static ScaledUnit GigaPascal;

  	static ScaledUnit Torr;
  	static ScaledUnit mTorr;
  	static ScaledUnit yTorr;

	static ScaledUnit PSI;
	static ScaledUnit KSI;

	static ScaledUnit Watt;
	static ScaledUnit VoltAmpere;

	static ScaledUnit Joule;
	static ScaledUnit NewtonMeter;
	static ScaledUnit VoltAmpereSecond;
	static ScaledUnit WattSecond;

    static ScaledUnit KMH;
	static ScaledUnit MPH;

	static ScaledUnit Degree;
	static ScaledUnit Radian;
	static ScaledUnit Gon;


    //@}


protected:
    double _ScaleFactor;
    Unit   _Unit;

};

} // namespace SketcherExpressions



#endif /* SKETCHER_SCALEDUNIT_H_ */
