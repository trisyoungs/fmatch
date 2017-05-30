/*
	*** Column-Major (OpenGL-friendly) 4x4 Matrix class
	*** matrix.cpp
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

#include "matrix.h"

// Constructor
Matrix::Matrix()
{
	setIdentity();
}

/*
// Operators
*/

// Matrix multiply (operator *) (return new matrix)
Matrix Matrix::operator*(const Matrix &B) const
{
	// [ row(A|this).column(B) ]
	Matrix AB;
	AB.matrix_[0] = matrix_[0]*B.matrix_[0] + matrix_[4]*B.matrix_[1] + matrix_[8]*B.matrix_[2] + matrix_[12]*B.matrix_[3];
	AB.matrix_[1] = matrix_[1]*B.matrix_[0] + matrix_[5]*B.matrix_[1] + matrix_[9]*B.matrix_[2] + matrix_[13]*B.matrix_[3];
	AB.matrix_[2] = matrix_[2]*B.matrix_[0] + matrix_[6]*B.matrix_[1] + matrix_[10]*B.matrix_[2] + matrix_[14]*B.matrix_[3];
	AB.matrix_[3] = matrix_[3]*B.matrix_[0] + matrix_[7]*B.matrix_[1] + matrix_[11]*B.matrix_[2] + matrix_[15]*B.matrix_[3];

	AB.matrix_[4] = matrix_[0]*B.matrix_[4] + matrix_[4]*B.matrix_[5] + matrix_[8]*B.matrix_[6] + matrix_[12]*B.matrix_[7];
	AB.matrix_[5] = matrix_[1]*B.matrix_[4] + matrix_[5]*B.matrix_[5] + matrix_[9]*B.matrix_[6] + matrix_[13]*B.matrix_[7];
	AB.matrix_[6] = matrix_[2]*B.matrix_[4] + matrix_[6]*B.matrix_[5] + matrix_[10]*B.matrix_[6] + matrix_[14]*B.matrix_[7];
	AB.matrix_[7] = matrix_[3]*B.matrix_[4] + matrix_[7]*B.matrix_[5] + matrix_[11]*B.matrix_[6] + matrix_[15]*B.matrix_[7];

	AB.matrix_[8] = matrix_[0]*B.matrix_[8] + matrix_[4]*B.matrix_[9] + matrix_[8]*B.matrix_[10] + matrix_[12]*B.matrix_[11];
	AB.matrix_[9] = matrix_[1]*B.matrix_[8] + matrix_[5]*B.matrix_[9] + matrix_[9]*B.matrix_[10] + matrix_[13]*B.matrix_[11];
	AB.matrix_[10] = matrix_[2]*B.matrix_[8] + matrix_[6]*B.matrix_[9] + matrix_[10]*B.matrix_[10] + matrix_[14]*B.matrix_[11];
	AB.matrix_[11] = matrix_[3]*B.matrix_[8] + matrix_[7]*B.matrix_[9] + matrix_[11]*B.matrix_[10] + matrix_[15]*B.matrix_[11];

	AB.matrix_[12] = matrix_[0]*B.matrix_[12] + matrix_[4]*B.matrix_[13] + matrix_[8]*B.matrix_[14] + matrix_[12]*B.matrix_[15];
	AB.matrix_[13] = matrix_[1]*B.matrix_[12] + matrix_[5]*B.matrix_[13] + matrix_[9]*B.matrix_[14] + matrix_[13]*B.matrix_[15];
	AB.matrix_[14] = matrix_[2]*B.matrix_[12] + matrix_[6]*B.matrix_[13] + matrix_[10]*B.matrix_[14] + matrix_[14]*B.matrix_[15];
	AB.matrix_[15] = matrix_[3]*B.matrix_[12] + matrix_[7]*B.matrix_[13] + matrix_[11]*B.matrix_[14] + matrix_[15]*B.matrix_[15];
	return AB;
}

Matrix Matrix::operator*(const double a) const
{
	Matrix AB;
	for (int n=0; n<16; ++n) AB.matrix_[n] = matrix_[n] * a;
	return AB;
}

Matrix Matrix::operator+(const Matrix &B) const
{
	Matrix A;
	for (int n=0; n<16; ++n) A[n] = matrix_[n] + B.matrix_[n];
	return A;
}

Matrix Matrix::operator-(const Matrix &B) const
{
	Matrix A;
	for (int n=0; n<16; ++n) A[n] = matrix_[n] - B.matrix_[n];
	return A;
}

Vec3<double> Matrix::operator*(const Vec3<double> &v) const
{
	Vec3<double> result;
	result.x = v.x*matrix_[0] + v.y*matrix_[4] + v.z*matrix_[8] + matrix_[12];
	result.y = v.x*matrix_[1] + v.y*matrix_[5] + v.z*matrix_[9] + matrix_[13];
	result.z = v.x*matrix_[2] + v.y*matrix_[6] + v.z*matrix_[10] + matrix_[14];
	return result;
}

// Matrix multiply (operator *=)
Matrix &Matrix::operator*=(const Matrix &B)
{
	// [ row(A|this).column(B) ]
	Matrix AB;
	AB.matrix_[0] = matrix_[0]*B.matrix_[0] + matrix_[4]*B.matrix_[1] + matrix_[8]*B.matrix_[2] + matrix_[12]*B.matrix_[3];
	AB.matrix_[1] = matrix_[1]*B.matrix_[0] + matrix_[5]*B.matrix_[1] + matrix_[9]*B.matrix_[2] + matrix_[13]*B.matrix_[3];
	AB.matrix_[2] = matrix_[2]*B.matrix_[0] + matrix_[6]*B.matrix_[1] + matrix_[10]*B.matrix_[2] + matrix_[14]*B.matrix_[3];
	AB.matrix_[3] = matrix_[3]*B.matrix_[0] + matrix_[7]*B.matrix_[1] + matrix_[11]*B.matrix_[2] + matrix_[15]*B.matrix_[3];

	AB.matrix_[4] = matrix_[0]*B.matrix_[4] + matrix_[4]*B.matrix_[5] + matrix_[8]*B.matrix_[6] + matrix_[12]*B.matrix_[7];
	AB.matrix_[5] = matrix_[1]*B.matrix_[4] + matrix_[5]*B.matrix_[5] + matrix_[9]*B.matrix_[6] + matrix_[13]*B.matrix_[7];
	AB.matrix_[6] = matrix_[2]*B.matrix_[4] + matrix_[6]*B.matrix_[5] + matrix_[10]*B.matrix_[6] + matrix_[14]*B.matrix_[7];
	AB.matrix_[7] = matrix_[3]*B.matrix_[4] + matrix_[7]*B.matrix_[5] + matrix_[11]*B.matrix_[6] + matrix_[15]*B.matrix_[7];

	AB.matrix_[8] = matrix_[0]*B.matrix_[8] + matrix_[4]*B.matrix_[9] + matrix_[8]*B.matrix_[10] + matrix_[12]*B.matrix_[11];
	AB.matrix_[9] = matrix_[1]*B.matrix_[8] + matrix_[5]*B.matrix_[9] + matrix_[9]*B.matrix_[10] + matrix_[13]*B.matrix_[11];
	AB.matrix_[10] = matrix_[2]*B.matrix_[8] + matrix_[6]*B.matrix_[9] + matrix_[10]*B.matrix_[10] + matrix_[14]*B.matrix_[11];
	AB.matrix_[11] = matrix_[3]*B.matrix_[8] + matrix_[7]*B.matrix_[9] + matrix_[11]*B.matrix_[10] + matrix_[15]*B.matrix_[11];

	AB.matrix_[12] = matrix_[0]*B.matrix_[12] + matrix_[4]*B.matrix_[13] + matrix_[8]*B.matrix_[14] + matrix_[12]*B.matrix_[15];
	AB.matrix_[13] = matrix_[1]*B.matrix_[12] + matrix_[5]*B.matrix_[13] + matrix_[9]*B.matrix_[14] + matrix_[13]*B.matrix_[15];
	AB.matrix_[14] = matrix_[2]*B.matrix_[12] + matrix_[6]*B.matrix_[13] + matrix_[10]*B.matrix_[14] + matrix_[14]*B.matrix_[15];
	AB.matrix_[15] = matrix_[3]*B.matrix_[12] + matrix_[7]*B.matrix_[13] + matrix_[11]*B.matrix_[14] + matrix_[15]*B.matrix_[15];
	*this = AB;
	return *this;
}

double &Matrix::operator[](int index)
{
	return matrix_[index];
}

/*
// Basic Set/Get
*/

// Reset to the identity matrix
void Matrix::setIdentity()
{
	matrix_[0] = 1.0;
	matrix_[1] = 0.0;
	matrix_[2] = 0.0;
	matrix_[3] = 0.0;
	matrix_[4] = 0.0;
	matrix_[5] = 1.0;
	matrix_[6] = 0.0;
	matrix_[7] = 0.0;
	matrix_[8] = 0.0;
	matrix_[9] = 0.0;
	matrix_[10] = 1.0;
	matrix_[11] = 0.0;
	matrix_[12] = 0.0;
	matrix_[13] = 0.0;
	matrix_[14] = 0.0;
	matrix_[15] = 1.0;
}

// Print matrix
void Matrix::print() const
{
	printf("CMaj   [0123]    [4567]    [8901] Translate\n");
	printf("        %8.4f %8.4f %8.4f %8.4f\n", matrix_[0], matrix_[4], matrix_[8], matrix_[12]);
	printf("        %8.4f %8.4f %8.4f %8.4f\n", matrix_[1], matrix_[5], matrix_[9], matrix_[13]);
	printf("        %8.4f %8.4f %8.4f %8.4f\n", matrix_[2], matrix_[6], matrix_[10], matrix_[14]);
	printf("Scale   %8.4f %8.4f %8.4f %8.4f\n", matrix_[3], matrix_[7], matrix_[11], matrix_[15]);
}

// Set zero matrix
void Matrix::zero()
{
	matrix_[0] = 0.0;
	matrix_[1] = 0.0;
	matrix_[2] = 0.0;
	matrix_[3] = 0.0;
	matrix_[4] = 0.0;
	matrix_[5] = 0.0;
	matrix_[6] = 0.0;
	matrix_[7] = 0.0;
	matrix_[8] = 0.0;
	matrix_[9] = 0.0;
	matrix_[10] = 0.0;
	matrix_[11] = 0.0;
	matrix_[12] = 0.0;
	matrix_[13] = 0.0;
	matrix_[14] = 0.0;
	matrix_[15] = 0.0;
}


// Return matrix array
double *Matrix::matrix()
{
	return matrix_;
}

// Return transpose of current matrix
Matrix &Matrix::transpose()
{
	static Matrix A;
	A.matrix_[0] = matrix_[0];
	A.matrix_[1] = matrix_[4];
	A.matrix_[2] = matrix_[8];
	A.matrix_[3] = matrix_[12];
	A.matrix_[4] = matrix_[1];
	A.matrix_[5] = matrix_[5];
	A.matrix_[6] = matrix_[9];
	A.matrix_[7] = matrix_[13];
	A.matrix_[8] = matrix_[2];
	A.matrix_[9] = matrix_[6];
	A.matrix_[10] = matrix_[10];
	A.matrix_[11] = matrix_[14];
	A.matrix_[12] = matrix_[3];
	A.matrix_[13] = matrix_[7];
	A.matrix_[14] = matrix_[11];
	A.matrix_[15] = matrix_[15];
	return A;
}

// Calculate determinant
double Matrix::determinant()
{

	double a = matrix_[0] * (matrix_[5]*(matrix_[10]*matrix_[15]-matrix_[11]*matrix_[14]) - matrix_[9]*(matrix_[6]*matrix_[15]-matrix_[7]*matrix_[14]) + matrix_[13]*(matrix_[6]*matrix_[11]-matrix_[7]*matrix_[10]) );
	double b = matrix_[4] * (matrix_[1]*(matrix_[10]*matrix_[15]-matrix_[11]*matrix_[14]) - matrix_[9]*(matrix_[2]*matrix_[15]-matrix_[3]*matrix_[14]) + matrix_[13]*(matrix_[2]*matrix_[11]-matrix_[3]*matrix_[10]) );
	double c = matrix_[8] * (matrix_[1]*(matrix_[6]*matrix_[15]-matrix_[7]*matrix_[14]) - matrix_[5]*(matrix_[2]*matrix_[15]-matrix_[3]*matrix_[14]) + matrix_[13]*(matrix_[2]*matrix_[7]-matrix_[3]*matrix_[6]) );
	double d = matrix_[12] * (matrix_[1]*(matrix_[6]*matrix_[11]-matrix_[7]*matrix_[10]) - matrix_[5]*(matrix_[2]*matrix_[11]-matrix_[3]*matrix_[10]) + matrix_[9]*(matrix_[2]*matrix_[7]-matrix_[3]*matrix_[6]) );
	return (a-b+c-d);
}

// Invert matrix
void Matrix::invert()
{
	// Gauss-Jordan Inversion
	// Invert the supplied matrix using Gauss-Jordan elimination
	int pivotrows[4], pivotcols[4], pivotrow = 0, pivotcol = 0;
	bool pivoted[4];
	int row, col, n, m;
	double large, element;
	for (n=0; n<4; ++n)
	{
		pivotrows[n] = 0;
		pivotcols[n] = 0;
		pivoted[n] = FALSE;
	}
	// Loop over columns to be reduced
	for (n=0; n<4; ++n)
	{
		// Locate suitable pivot element - find largest value in the matrix A
		large = 0.0;
		for (row=0; row<4; ++row)
		{
			// Only search this row if it has not previously contained a pivot element
			if (pivoted[row]) continue;
			for (col=0; col<4; ++col)
			{
				// Similarly, only look at the column element if the column hasn't been pivoted yet.
				if (pivoted[col]) continue;
				// Check the size of the element...
				element = fabs(matrix_[row*4+col]);
				if (element > large)
				{
					large = element;
					pivotrow = row;
					pivotcol = col;
				}
			}
		}
		
		// Mark the pivot row/column as changed
		pivoted[pivotcol] = TRUE;
		pivotrows[n] = pivotrow;
		pivotcols[n] = pivotcol;
		
		// Exchange rows to put pivot element on the diagonal
		if (pivotrow != pivotcol)
		{
			for (m=0; m<4; ++m)
			{
				element = matrix_[pivotrow*4+m];
				matrix_[pivotrow*4+m] = matrix_[pivotcol*4+m];
				matrix_[pivotcol*4+m] = element;
			}
		}
		
		// Now ready to divide through row elements.
		element = 1.0 / matrix_[pivotcol*4+pivotcol];
		matrix_[pivotcol*4+pivotcol] = 1.0;
		for (m=0; m<4; ++m) matrix_[pivotcol*4+m] *= element;
		
		// Divide through other rows by the relevant multiple of the pivot row
		for (row=0; row<4; ++row)
		{
			if (row == pivotcol) continue;
			element = matrix_[row*4 + pivotcol];
			matrix_[row*4 + pivotcol] = 0.0;
			for (m=0; m<4; ++m) matrix_[row*4+m] = matrix_[row*4+m] - matrix_[pivotcol*4+m] * element;
		}
	}
	// Rearrange columns to undo row exchanges performed earlier
	for (n=3; n>=0; --n)
	{
		if (pivotrows[n] != pivotcols[n]) for (m=0; m<4; ++m)
		{
			element = matrix_[m*4+pivotrows[n]];
			matrix_[m*4+pivotrows[n]] = matrix_[m*4+pivotcols[n]];
			matrix_[m*4+pivotcols[n]] = element;
		}
	}
}

/*
// Column Operations
*/

// Copy column contents to supplied Vec3
Vec3<double> Matrix::columnAsVec3(int col)
{
	Vec3<double> vec(matrix_[col*4], matrix_[col*4+1], matrix_[col*4+2]);
	return vec;
}

// Set specified row from supplied triplet of values
void Matrix::setRow(int row, double x, double y, double z)
{
	matrix_[row] = x;
	matrix_[4+row] = y;
	matrix_[8+row] = z;
}

// Set specified row from supplied values
void Matrix::setRow(int row, double x, double y, double z, double w)
{
	matrix_[row] = x;
	matrix_[4+row] = y;
	matrix_[8+row] = z;
	matrix_[12+row] = w;
}

// Set specified column from supplied values
void Matrix::setColumn(int col, double a, double b, double c, double d)
{
	matrix_[col*4] = a;
	matrix_[col*4+1] = b;
	matrix_[col*4+2] = c;
	matrix_[col*4+3] = d;
}

// Set specified column from supplied Vec3
void Matrix::setColumn(int col, Vec3<double> vec, double w)
{
	matrix_[col*4] = vec.x;
	matrix_[col*4+1] = vec.y;
	matrix_[col*4+2] = vec.z;
	matrix_[col*4+3] = w;
}

// Adjust specified column from supplied values
void Matrix::adjustColumn(int col, double a, double b, double c, double d)
{
	matrix_[col*4] += a;
	matrix_[col*4+1] += b;
	matrix_[col*4+2] += c;
	matrix_[col*4+3] += d;
}

// Adjust specified column from supplied Vec3
void Matrix::adjustColumn(int col, Vec3<double> vec, double w)
{
	matrix_[col*4] += vec.x;
	matrix_[col*4+1] += vec.y;
	matrix_[col*4+2] += vec.z;
	matrix_[col*4+3] += w;
}

// Calculate column magnitude
double Matrix::columnMagnitude(int column)
{
	double mag = 0.0;
	for (int n=column*4; n<column*4+4; ++n) mag += (matrix_[n] * matrix_[n]);
	return sqrt(mag);
}

// Multiply column by single value
void Matrix::columnMultiply(int col, double d)
{
	matrix_[col*4] *= d;
	matrix_[col*4+1] *= d;
	matrix_[col*4+2] *= d;
	matrix_[col*4+3] *= d;
}

// Multiply first three columns by values insupplied vector
void Matrix::columnMultiply(Vec3<double> vec)
{
	columnMultiply(0, vec.x);
	columnMultiply(1, vec.y);
	columnMultiply(2, vec.z);
}

// Normalise specified column to 1
void Matrix::columnNormalise(int col)
{
	double mag = 1.0/sqrt(matrix_[col*4]*matrix_[col*4] + matrix_[col*4+1]*matrix_[col*4+1] + matrix_[col*4+2]*matrix_[col*4+2] + matrix_[col*4+3]*matrix_[col*4+3]);
	matrix_[col*4] *= mag;
	matrix_[col*4+1] *= mag;
	matrix_[col*4+2] *= mag;
	matrix_[col*4+3] *= mag;
}

// Orthogonalise rotation matrix column w.r.t. one (or two) other columns)
void Matrix::orthogonaliseColumn(int targetcol, int orthocol1, int orthocol2)
{
	// Grab target column
	Vec3<double> v = columnAsVec3(targetcol);
	// Orthogonalising w.r.t one or two other vectors?
	if (orthocol2 == -1)
	{
		Vec3<double> source = columnAsVec3(orthocol1);
		double sourcemag = source.magnitude();
		double dpovermagsq = v.dp(source) / (sourcemag * sourcemag);
		v.x -= dpovermagsq * source.x;
		v.y -= dpovermagsq * source.y;
		v.z -= dpovermagsq * source.z;
	}
	else
	{
		// This routine actually generates the orthogonal vector via the cross-product
		// We also calculate the scalar resolute (dp) to ensure the new vector points in the same direction
		Vec3<double> source1 = columnAsVec3(orthocol1), source2 = columnAsVec3(orthocol2);
		Vec3<double> newvec = source1 * source2;
		newvec.normalise();
		double dp = newvec.dp(v);
		if (dp < 0.0) newvec *= -1.0;
		v = newvec;
	}
	setColumn(targetcol, v, matrix_[targetcol*4+3]);
}

/*
// Misc
*/

// Transform coordinates supplied and return as Vec3<double>
Vec3<double> Matrix::transform(double x, double y, double z) const
{
	Vec3<double> result;
	result.x = x*matrix_[0] + y*matrix_[4] + z*matrix_[8] + matrix_[12];
	result.y = x*matrix_[1] + y*matrix_[5] + z*matrix_[9] + matrix_[13];
	result.z = x*matrix_[2] + y*matrix_[6] + z*matrix_[10] + matrix_[14];
	return result;
}

// Transform coordinates supplied and return as Vec3<double>
Vec3<double> Matrix::transform(Vec3<double> vec) const
{
	Vec3<double> result;
	result.x = vec.x*matrix_[0] + vec.y*matrix_[4] + vec.z*matrix_[8] + matrix_[12];
	result.y = vec.x*matrix_[1] + vec.y*matrix_[5] + vec.z*matrix_[9] + matrix_[13];
	result.z = vec.x*matrix_[2] + vec.y*matrix_[6] + vec.z*matrix_[10] + matrix_[14];
	return result;
}

// Multiply against other matrix, but only rotational part, keeping translation/scaling intact
void Matrix::multiplyRotation(Matrix B)
{
	Matrix AB;
	AB.matrix_[0] = matrix_[0]*B.matrix_[0] + matrix_[4]*B.matrix_[1] + matrix_[8]*B.matrix_[2];
	AB.matrix_[1] = matrix_[1]*B.matrix_[0] + matrix_[5]*B.matrix_[1] + matrix_[9]*B.matrix_[2];
	AB.matrix_[2] = matrix_[2]*B.matrix_[0] + matrix_[6]*B.matrix_[1] + matrix_[10]*B.matrix_[2];
	
	AB.matrix_[4] = matrix_[0]*B.matrix_[4] + matrix_[4]*B.matrix_[5] + matrix_[8]*B.matrix_[6];
	AB.matrix_[5] = matrix_[1]*B.matrix_[4] + matrix_[5]*B.matrix_[5] + matrix_[9]*B.matrix_[6];
	AB.matrix_[6] = matrix_[2]*B.matrix_[4] + matrix_[6]*B.matrix_[5] + matrix_[10]*B.matrix_[6];
	
	AB.matrix_[8] = matrix_[0]*B.matrix_[8] + matrix_[4]*B.matrix_[9] + matrix_[8]*B.matrix_[10];
	AB.matrix_[9] = matrix_[1]*B.matrix_[8] + matrix_[5]*B.matrix_[9] + matrix_[9]*B.matrix_[10];
	AB.matrix_[10] = matrix_[2]*B.matrix_[8] + matrix_[6]*B.matrix_[9] + matrix_[10]*B.matrix_[10];
	
	matrix_[0] = AB.matrix_[0];
	matrix_[1] = AB.matrix_[1];
	matrix_[2] = AB.matrix_[2];
	matrix_[4] = AB.matrix_[4];
	matrix_[5] = AB.matrix_[5];
	matrix_[6] = AB.matrix_[6];
	matrix_[8] = AB.matrix_[8];
	matrix_[9] = AB.matrix_[9];
	matrix_[10] = AB.matrix_[10];
}

// Apply rotational part of matrix to supplied vector
Vec3<double> Matrix::rotateVector(Vec3<double> &v) const
{
	Vec3<double> result;
	result.x = v.x*matrix_[0] + v.y*matrix_[4] + v.z*matrix_[8];
	result.y = v.x*matrix_[1] + v.y*matrix_[5] + v.z*matrix_[9];
	result.z = v.x*matrix_[2] + v.y*matrix_[6] + v.z*matrix_[10];
	return result;
}

// Apply rotational part of matrix to supplied vector coordinates
Vec3<double> Matrix::rotateVector(double x, double y, double z) const
{
	Vec3<double> result;
	result.x = x*matrix_[0] + y*matrix_[4] + z*matrix_[8];
	result.y = x*matrix_[1] + y*matrix_[5] + z*matrix_[9];
	result.z = x*matrix_[2] + y*matrix_[6] + z*matrix_[10];
	return result;
}

// Remove translation and scaling parts, leaving rotation only
void Matrix::removeTranslationAndScaling()
{
	matrix_[3] = 0.0;
	matrix_[7] = 0.0;
	matrix_[11] = 0.0;
	matrix_[15] = 1.0;
	matrix_[12] = 0.0;
	matrix_[13] = 0.0;
	matrix_[14] = 0.0;
}

// Copy translation and scaling parts from specified matrix
void Matrix::copyTranslationAndScaling(Matrix &source)
{
	matrix_[3] = source.matrix_[3];
	matrix_[7] = source.matrix_[7];
	matrix_[11] = source.matrix_[11];
	matrix_[15] = source.matrix_[15];
	matrix_[12] = source.matrix_[12];
	matrix_[13] = source.matrix_[13];
	matrix_[14] = source.matrix_[14];
}


// Returns a unit vector in the specified direction
Vec3<double> unit_vector(int n)
{
	Vec3<double> result;
	result.zero();
	result.set(n,1.0);
	return result;
}

// Create 3x3 matrix from cyclic permutations of supplied vector
void Matrix::cyclicPermute(Vec3<double> v)
{
	Vec3<double> temp;
	for (int n=0; n<3; ++n)
	{
		temp = unit_vector((n+1)%3) * v.get((n+2)%3) - unit_vector((n+2)%3) * v.get((n+1)%3);
		setColumn(n,temp.x,temp.y,temp.z,0.0);
	}
	setColumn(3, 0.0, 0.0, 0.0, 1.0);
}
