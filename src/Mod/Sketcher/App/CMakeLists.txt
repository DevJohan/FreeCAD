cmake_minimum_required(VERSION 2.8)
find_package(BISON)
find_package(FLEX)

if(MSVC)
    add_definitions(-DFCAppSketcher -DHAVE_ACOSH -DHAVE_ASINH -DHAVE_ATANH)
else(MSVC)
    add_definitions(-DHAVE_LIMITS_H -DHAVE_CONFIG_H)
endif(MSVC)

include_directories(
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_BINARY_DIR}/..
    ${Boost_INCLUDE_DIRS}
    ${OCC_INCLUDE_DIR}
    ${QT_QTCORE_INCLUDE_DIR}
    ${ZLIB_INCLUDE_DIR}
    ${PYTHON_INCLUDE_PATH}
    ${XERCESC_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)
link_directories(${OCC_LIBRARY_DIR})

set(Sketcher_LIBS
    Part
    FreeCADApp
)

generate_from_xml(SketchObjectSFPy)
generate_from_xml(SketchObjectPy)
generate_from_xml(ConstraintPy)
generate_from_xml(SketchPy)
generate_from_xml(SketchSolverPy)


SET(Properties_SRCS
    PropertyConstraintList.cpp
    PropertyConstraintList.h
)
SOURCE_GROUP("Properties" FILES ${Properties_SRCS})


SET(Features_SRCS
    SketchObjectSF.cpp
    SketchObjectSF.h
    SketchObject.cpp
    SketchObject.h
)
SOURCE_GROUP("Features" FILES ${Features_SRCS})

SET(Datatypes_SRCS
    Constraint.cpp
    Constraint.h
    Sketch.cpp
    Sketch.h
    Sketch_exp.cpp
    Sketch_exp.h
    SketchSolver.cpp
    SketchSolver.h
)
SOURCE_GROUP("Datatypes" FILES ${Datatypes_SRCS})

SET(Python_SRCS
    SketchObjectSFPy.xml
    SketchObjectSFPyImp.cpp
    SketchObjectPy.xml
    SketchObjectPyImp.cpp
    ConstraintPyImp.cpp
    ConstraintPy.xml
    SketchPy.xml
    SketchPyImp.cpp
    SketchSolverPy.xml
    SketchSolverPyImp.cpp
)
SOURCE_GROUP("Python" FILES ${Python_SRCS})


SET(FreeGCS_SRCS
    freegcs/qp_eq.h
    freegcs/GCS.cpp
    freegcs/GCS.h
    freegcs/Util.h
    freegcs/Geo.h
    freegcs/Constraints.cpp
    freegcs/Constraints.h
    freegcs/SubSystem.cpp
    freegcs/SubSystem.h
    freegcs/qp_eq.cpp
    freegcs/qp_eq.h
)
SOURCE_GROUP("FreeGCS" FILES ${FreeGCS_SRCS})

SET(FreeGCS_EXP_SRCS
    freegcs_exp/GCS_exp.cpp
    freegcs_exp/GCS_exp.h
    freegcs_exp/Util_exp.h
    freegcs_exp/Geo_exp.h
    freegcs_exp/Constraints_exp.cpp
    freegcs_exp/Constraints_exp.h
    freegcs_exp/ConstraintsInfo_exp.cpp
    freegcs_exp/ConstraintsInfo_exp.h
    freegcs_exp/SubSystem_exp.cpp
    freegcs_exp/SubSystem_exp.h
    freegcs_exp/qp_eq_exp.cpp
    freegcs_exp/qp_eq_exp.h
)
SOURCE_GROUP("FreeGCS_EXP" FILES ${FreeGCS_EXP_SRCS})


BISON_TARGET( ExpressionParser 
	${CMAKE_CURRENT_SOURCE_DIR}/expressions/ExpressionParser.y 
	${CMAKE_CURRENT_BINARY_DIR}/expressions/ExpressionParser.cpp )
FLEX_TARGET( ExpressionLexer  
	${CMAKE_CURRENT_SOURCE_DIR}/expressions/ExpressionLexer.l  
	${CMAKE_CURRENT_BINARY_DIR}/expressions/ExpressionLexer.cpp
	COMPILE_FLAGS --header-file=${CMAKE_CURRENT_BINARY_DIR}/expressions/ExpressionLexer.h )
ADD_FLEX_BISON_DEPENDENCY(ExpressionLexer ExpressionParser)
SET(ExpressionParser_SRCS
	${BISON_ExpressionParser_OUTPUTS} 
	${FLEX_ExpressionLexer_OUTPUTS}
)
fc_check_target_dirs( ${ExpressionParser_SRCS} )

SET(Expressions_SRCS
	${ExpressionParser_SRCS}
	expressions/ScaledUnit.cpp
	expressions/ScaledUnit.h
	expressions/ExpressionParserContext.cpp
	expressions/ExpressionParserContext.h
	expressions/ExpressionBase.cpp
	expressions/ExpressionBase.h
	expressions/FunctionExpressions.cpp
	expressions/FunctionExpressions.h
	expressions/ExpressionContext.cpp
	expressions/ExpressionContext.h
	expressions/SumExpression.cpp
	expressions/SumExpression.h
	expressions/ProductExpression.cpp
	expressions/ProductExpression.h
)
SOURCE_GROUP("Expressions" FILES ${Expressions_SRCS})

SET(SketchModule_SRCS
    AppSketcher.cpp
    AppSketcherPy.cpp
    PreCompiled.cpp
    PreCompiled.h
)
SOURCE_GROUP("Module" FILES ${SketchModule_SRCS})

SET(Sketcher_SRCS
    ${Features_SRCS}
    ${FreeGCS_SRCS}
    ${FreeGCS_EXP_SRCS}
    ${SketchModule_SRCS}
    ${Python_SRCS}
    ${Properties_SRCS}
    ${Datatypes_SRCS}
    ${Expressions_SRCS}
)

SET(Sketcher_Scripts
    Init.py
    SketcherExample.py
	TestSketcherApp.py
	Profiles.py
	ProfileLib/Hexagon.py
	ProfileLib/RegularPolygon.py
	ProfileLib/__init__.py
)

add_library(Sketcher SHARED ${Sketcher_SRCS})
target_link_libraries(Sketcher ${Sketcher_LIBS})


fc_target_copy_resource(Sketcher 
    ${CMAKE_SOURCE_DIR}/src/Mod/Sketcher
    ${CMAKE_BINARY_DIR}/Mod/Sketcher
    ${Sketcher_Scripts})

SET_BIN_DIR(Sketcher Sketcher /Mod/Sketcher)
SET_PYTHON_PREFIX_SUFFIX(Sketcher)

INSTALL(TARGETS Sketcher DESTINATION ${CMAKE_INSTALL_LIBDIR})
