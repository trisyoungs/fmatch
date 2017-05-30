/*
	*** Column-Major (OpenGL-friendly) 4x4 Matrix class
	*** matrix.h
	Copyright T. Youngs 2011

	This file is part of FMatch3.

	FMatch3 is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	FMatch3 is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with FMatch3.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FMATCH3_MATRIX_H
#define FMATCH3_MATRIX_H

#include "vector3.h"

// Column-major 4x4 matrix
class Matrix
{
	public:
	// Constructor
	Matrix();
	
	private:
	// Matrix
	double matrix_[16];


	/*
	// Operators
	*/
	public:
	Matrix operator*(const Matrix &B) const;
	Matrix operator*(const double a) const;
	Matrix operator+(const Matrix &B) const;
	Matrix operator-(const Matrix &B) const;
	Vec3<double> operator*(const Vec3<double> &v) const;
	Matrix &operator*=(const Matrix &B);
	double &operator[](int);


	/*
	// Basic Set/Get
	*/
	public:
	// Reset the matrix to the identity
	void setIdentity();
	// Prints the matrix to stdout
	void print() const;
	// Set the zero matrix
	void zero();
	// Return matrix array
	double *matrix();
	// Return transpose of current matrix
	Matrix &transpose();
	// Calculate determinant
	double determinant();
	// Invert matrix
	void invert();


	/*
	// Column Operations
	*/
	public:
	// Copy column contents to supplied Vec3
	Vec3<double> columnAsVec3(int col);
	// Set specified row from supplied triplet of values
	void setRow(int row, double x, double y, double z);
	// Set specified row from supplied values
	void setRow(int row, double x, double y, double z, double w);
	// Set specified column from supplied values
	void setColumn(int col, double a, double b, double c, double d);
	// Set specified column from supplied Vec3
	void setColumn(int col, Vec3<double> vec, double w);
	// Adjust specified column from supplied values
	void adjustColumn(int col, double a, double b, double c, double d);
	// Adjust specified column from supplied Vec3
	void adjustColumn(int col, Vec3<double> vec, double w);
	// Calculate column magnitude
	double columnMagnitude(int column);
	// Multiply single column by single value
	void columnMultiply(int col, double d);
	// Multiply first three columns by values insupplied vector
	void columnMultiply(Vec3<double> vec);
	// Normalise specified column to 1
	void columnNormalise(int column);
	// Orthogonalise rotation matrix column w.r.t. one (or two) other columns)
	void orthogonaliseColumn(int targetcol, int orthcol1, int orthocol2 = -1);


	/*
	// Misc
	*/
	public:
	// Transform coordinates supplied and return as Vec3<double>
	Vec3<double> transform(double x, double y, double z) const;
	// Transform coordinates supplied and return as Vec3<double>
	Vec3<double> transform(Vec3<double> vec) const;
	// Apply rotational part of matrix to supplied Vec3
	Vec3<double> rotateVector(Vec3<double> &v) const;
	// Apply rotational part of matrix to supplied vector coordinates
	Vec3<double> rotateVector(double x, double y, double z) const;
	// Multiply against other matrix, but only rotational part, keeping translation/scaling intact
	void multiplyRotation(Matrix B);
	// Remove translation and scaling parts, leaving rotation only
	void removeTranslationAndScaling();
	// Copy translation and scaling parts from specified matrix
	void copyTranslationAndScaling(Matrix &source);
	// Create 3x3 matrix from cyclic permutations of supplied vector
	void cyclicPermute(Vec3<double> v);	
};

#endif

