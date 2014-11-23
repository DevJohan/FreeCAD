/***************************************************************************
 *   Copyright (c) J�rgen Riegel          (juergen.riegel@web.de) 2010     *
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

#ifndef SKETCHER_SKETCHSOLVER_H
#define SKETCHER_SKETCHSOLVER_H

#include <App/PropertyStandard.h>
#include <App/PropertyFile.h>
#include <Mod/Part/App/Geometry.h>
#include <Mod/Part/App/TopoShape.h>
#include "Constraint.h"

#include <Base/Persistence.h>

namespace Sketcher
{

class SketcherExport SketchSolver :public Base::Persistence
{
    TYPESYSTEM_HEADER();

public:
    SketchSolver();
    virtual ~SketchSolver() = 0;

    // from base class
//    virtual unsigned int getMemSize(void) const;
//    virtual void Save(Base::Writer &/*writer*/) const;
//    virtual void Restore(Base::XMLReader &/*reader*/);

    virtual int64_t getGeometryVersion() const = 0;
    virtual double getSolveTime() const = 0;
    virtual int getLatestAlgorithm() const = 0;

    /// solve the actual set up sketch
    virtual int solve(void) = 0;
    /// delete all geometry and constraints, leave an empty sketch
    virtual void clear(void) = 0;
    /** set the sketch up with geoms and constraints
      * 
      * returns the degree of freedom of a sketch and calculates a list of
      * conflicting constraints
      *
      * 0 degrees of freedom correspond to a fully constrained sketch
      * -1 degrees of freedom correspond to an over-constrained sketch
      * positive degrees of freedom correspond to an under-constrained sketch
      *
      * an over-constrained sketch will always contain conflicting constraints
      * a fully constrained or under-constrained sketch may contain conflicting
      * constraints or may not
      */

    virtual int setUpSketch(const int64_t new_geometry_version, const std::vector<Part::Geometry *> &GeoList, const std::vector<Constraint *> &ConstraintList,
                    int extGeoCount=0) = 0;
    /// return the actual geometry of the sketch a TopoShape
    virtual Part::TopoShape toShape(void) const = 0;
    /// add unspecified geometry
    virtual int addGeometry(const Part::Geometry *geo, bool fixed=false) = 0;
    /// add unspecified geometry
    virtual int addGeometry(const std::vector<Part::Geometry *> &geo, bool fixed=false) = 0;
    /// returns the actual geometry
    virtual std::vector<Part::Geometry *> extractGeometry(bool withConstrucionElements=true,
                                                  bool withExternalElements=false) const = 0;
    /// get the geometry as python objects
    virtual Py::Tuple getPyGeometry(void) const = 0;

    /// retrieves the index of a point
    virtual int getPointId(int geoId, PointPos pos) const = 0;
    /// retrieves a point
    virtual Base::Vector3d getPoint(int geoId, PointPos pos) const = 0;

    virtual bool hasConflicts(void) const = 0;
    virtual const std::vector<int> &getConflicting(void) const = 0;
    virtual bool hasRedundancies(void) const = 0;
    virtual const std::vector<int> &getRedundant(void) const = 0;

    /** set the datum of a distance or angle constraint to a certain value and solve
      * This can cause the solving to fail!
      */
    virtual int setDatum(int constrId, double value) = 0;

    /** initializes a point (or curve) drag by setting the current
      * sketch status as a reference
      */
    virtual int initMove(int geoId, PointPos pos, bool fine=true) = 0;

    /** move this point (or curve) to a new location and solve.
      * This will introduce some additional weak constraints expressing
      * a condition for satisfying the new point location!
      * The relative flag permits moving relatively to the current position
      */
    virtual int movePoint(int geoId, PointPos pos, Base::Vector3d toPoint, bool relative=false) = 0;

    /// add dedicated geometry
    //@{
    /// add a point
    virtual int addPoint(const Part::GeomPoint &point, bool fixed=false) = 0;
    /// add an infinite line
    virtual int addLine(const Part::GeomLineSegment &line, bool fixed=false) = 0;
    /// add a line segment
    virtual int addLineSegment(const Part::GeomLineSegment &lineSegment, bool fixed=false) = 0;
    /// add a arc (circle segment)
    virtual int addArc(const Part::GeomArcOfCircle &circleSegment, bool fixed=false) = 0;
    /// add a circle
    virtual int addCircle(const Part::GeomCircle &circle, bool fixed=false) = 0;
    /// add a ellipse
    virtual int addEllipse(const Part::GeomEllipse &ellipse, bool fixed=false) = 0;
    //@}


    /// constraints
    //@{
    /// add all constraints in the list
    virtual int addConstraints(const std::vector<Constraint *> &ConstraintList) = 0;
    /// add one constraint to the sketch
    virtual int addConstraint(const Constraint *constraint) = 0;
    /// add a fixed coordinate constraint to a point
    virtual int addCoordinateXConstraint(int geoId, PointPos pos, double value) = 0;
    virtual int addCoordinateYConstraint(int geoId, PointPos pos, double value) = 0;
    /// add a horizontal distance constraint to two points or line ends
    virtual int addDistanceXConstraint(int geoId, double value) = 0;
    virtual int addDistanceXConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2, double value) = 0;
    /// add a vertical distance constraint to two points or line ends
    virtual int addDistanceYConstraint(int geoId, double value) = 0;
    virtual int addDistanceYConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2, double value) = 0;
    /// add a horizontal constraint to a geometry
    virtual int addHorizontalConstraint(int geoId) = 0;
    virtual int addHorizontalConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2) = 0;
    /// add a vertical constraint to a geometry
    virtual int addVerticalConstraint(int geoId) = 0;
    virtual int addVerticalConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2) = 0;
    /// add a coincident constraint to two points of two geometries
    virtual int addPointCoincidentConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2) = 0;
    /// add a length or distance constraint
    virtual int addDistanceConstraint(int geoId1, double value) = 0;
    virtual int addDistanceConstraint(int geoId1, int geoId2, double value) = 0;
    virtual int addDistanceConstraint(int geoId1, PointPos pos1, int geoId2, double value) = 0;
    virtual int addDistanceConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2, double value) = 0;
    /// add a parallel constraint between two lines
    virtual int addParallelConstraint(int geoId1, int geoId2) = 0;
    /// add a perpendicular constraint between two lines
    virtual int addPerpendicularConstraint(int geoId1, int geoId2) = 0;
    virtual int addPerpendicularConstraint(int geoId1, PointPos pos1, int geoId2) = 0;
    virtual int addPerpendicularConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2) = 0;
    /// add a tangency constraint between two geometries
    virtual int addTangentConstraint(int geoId1, int geoId2) = 0;
    virtual int addTangentConstraint(int geoId1, PointPos pos1, int geoId2) = 0;
    virtual int addTangentConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2) = 0;
    /// add a radius constraint on a circle or an arc
    virtual int addRadiusConstraint(int geoId, double value) = 0;
    /// add an angle constraint on a line or between two lines
    virtual int addAngleConstraint(int geoId, double value) = 0;
    virtual int addAngleConstraint(int geoId1, int geoId2, double value) = 0;
    virtual int addAngleConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2, double value) = 0;
    /// add an equal length or radius constraints between two lines or between circles and arcs
    virtual int addEqualConstraint(int geoId1, int geoId2) = 0;
    /// add a point on line constraint
    virtual int addPointOnObjectConstraint(int geoId1, PointPos pos1, int geoId2) = 0;
    /// add a symmetric constraint between two points with respect to a line
    virtual int addSymmetricConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2, int geoId3) = 0;
    /// add a symmetric constraint between three points, the last point is in the middle of the first two
    virtual int addSymmetricConstraint(int geoId1, PointPos pos1, int geoId2, PointPos pos2, int geoId3, PointPos pos3) = 0;
    //@}
};

} //namespace Part


#endif // SKETCHER_SKETCHSOLVER_H
