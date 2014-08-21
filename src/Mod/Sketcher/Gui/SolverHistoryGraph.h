/*
 * SolverHistoryGraph.h
 *
 *  Created on: 16 jul 2014
 *      Author: johan
 */

#ifndef SOLVERHISTORYGRAPH_H_
#define SOLVERHISTORYGRAPH_H_

#include <qwidget.h>
#include <boost/signals2.hpp>

namespace SketcherGui{
#if __cplusplus < 201103L
#ifndef nullptr
#define nullptr 0
#endif
#endif

template<class T>
struct rgb_struct{
	T values[3];
	rgb_struct(){}
	rgb_struct(T r, T g, T b){values[0] = r; values[1] = g; values[2] = b; }
	rgb_struct(T gray){ values[0] = gray; values[1] = gray; values[2] = gray; }
	template<class S> rgb_struct<T>& operator =(const rgb_struct<S>& operand){
		values[0] = operand.values[0]; values[1] = operand.values[1]; values[2] = operand.values[2]; return *this;}
	template<class S> rgb_struct<T>& operator +=(const rgb_struct<S>& operand){
		values[0] += operand.values[0]; values[1] += operand.values[1]; values[2] += operand.values[2];	return *this; }
	template<class S> rgb_struct<T>& operator -=(const rgb_struct<S>& operand){
		values[0] -= operand.values[0]; values[1] -= operand.values[1]; values[2] -= operand.values[2]; return *this; }
	template<class S> rgb_struct<T>& operator *=(const S& operand){
		values[0] *= operand; values[1] *= operand; values[2] *= operand; return *this; }
	bool operator ==(const rgb_struct<T>& other){
		return values[0] == other.values[0] && values[1] == other.values[1] && values[2] == other.values[2]; }
	bool operator !=(const rgb_struct<T>& other){
		return values[0] != other.values[0] || values[1] != other.values[1] || values[2] != other.values[2]; }
};

template<class T, class S>
rgb_struct<T> operator +(const rgb_struct<T>& value, const rgb_struct<S>& operand)  {
	return rgb_struct<T>( value.values[0] - operand.values[0], value.values[1] - operand.values[1], value.values[2] - operand.values[2] ); }
template<class T, class S>
rgb_struct<T> operator -(const rgb_struct<T>& value, const rgb_struct<S>& operand)  {
	return rgb_struct<T>( value.values[0] - operand.values[0], value.values[1] - operand.values[1], value.values[2] - operand.values[2] ); }
template<class T>
rgb_struct<T> operator *(double d, const rgb_struct<T>& operand)  {
	return rgb_struct<T>( d * operand.values[0], d * operand.values[1], d * operand.values[2] ); }

typedef rgb_struct<uchar> rgb8;

template < typename T >
class SImage{
	static size_t ceil_to_multiple_of_16( size_t num ){ return num > 0 ? ((( num - 1 )>>4) + 1 ) << 4 : 0; }
public:
	SImage(int width, int height):
		_width(width),
		_height(height),
		_byte_stride( ceil_to_multiple_of_16( sizeof( T[ width ] ))),
		_data( new uchar[ _byte_stride * width ])
		{}
	~SImage(){ if( _data != nullptr ) delete[] _data; }
	T* row_data( int y ){ return reinterpret_cast<T*>(_data + y * _byte_stride);}
	const T* row_data( int y ) const { return reinterpret_cast<const T*>(_data + y * _byte_stride);}
	T& operator()(int y, int x){ return row_data(y)[ x ];}
	const T& operator()(int y, int x) const { return row_data(y)[ x ];}
	void fill( const T& value ){
		for( int x = 0; x < _width;  x++ ){ row_data(0)[x]=value; }
		for( int y = 1; y < _height; y++ ){ memcpy( row_data(y), row_data(0), sizeof( T[_width]) ); }}
	int width(){ return _width;}
	int height(){ return _height; }
	int byte_stride(){return _byte_stride;}
	uchar* data(){return _data;}
private:
	int _width;
	int _height;
	int _byte_stride;
	uchar* _data;
};

typedef boost::signals2::connection Connection;

class SolverHistoryGraph : public QWidget {
	Q_OBJECT
public:
	SolverHistoryGraph(QWidget* parent);
	virtual ~SolverHistoryGraph();

	void statusUpdate(		double solveTime, int solver, bool succeeded );
	void statusUpdate_exp(	double solveTime, int solver, bool succeeded );
protected:
	void paintEvent(QPaintEvent* event);
private:
	int value2line( double value );
	SImage<rgb8> imageData;
	QImage graph;
	int current_update;
	int current_update_exp;
};

}

#endif /* SOLVERHISTORYGRAPH_H_ */
