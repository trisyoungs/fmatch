/*
	*** Set of Alpha
	*** alpha.h
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

#ifndef FM3_ALPHA_H
#define FM3_ALPHA_H

#include "list.h"

// Forward Declarations
class Parameter;
class Forcefield;
class TargetSystem;

// Set of Variables to be Optimised (Alpha)
class Alpha
{
	public:
	// Constructor / Destructor
	Alpha(int size = 0);
	~Alpha();
	// Copy constructor
        Alpha(const Alpha&);
	// List pointers
	Alpha *prev, *next;


	/*
	// Data Array
	*/
	private:
	// Number of alpha values in Alpha
	int nAlpha_;
	// Data array
	double *alpha_;
	// Pointers to original parameters
	Parameter **sourceParameters_;
	// Length scales for individual parameters (to be multiplied by basic length scale)
	double *lengthScales_;
	// Fit flags for parameters
	bool *fit_;
	// Flag whether current cost is valid (i.e. parameters have not been changed)
	bool costIsValid_;
	// Cost function at this Alpha
	double cost_;

	private:
	// Initialise to specified size, freeing old arrays if necessary
	void initialise(int size);
	// Copy information from supplied Alpha
	void copyData(const Alpha &v);

	public:
	// Clear all arrays and data
	void clear();
	// Copy data from supplied Alpha
	void copy(Alpha *sorce);
	// Set from copy of supplied parameter list
	void set(List<Parameter> *parameters);
	// Set from copy of supplied Alpha (optionally with single value adjusted)
	void set(const Alpha &v, int nudgeId = -1, double nudgeFactor = 0.0);
	// Set single value/parameter in Alpha
	void set(int n, double value, Parameter *param, double lengthScale, bool fit);
	// Multiply value in Alpha by supplied value
	void multiplyValue(int n, double d);
	// Calculate and return cost function at this Alpha
	double cost(TargetSystem* system, Forcefield* ff);
	// Reset parameter length scales
	void resetLengthScales();
	// Calculate parameter length scales
	void calculateLengthScales(TargetSystem *system, Forcefield *ff);
	// Return length scale for specified parameter
	double lengthScale(int n);
	// Return number of parameters in array
	int nAlpha() const;
	// Return specified alpha
	double alpha(int n);
	// Set value of specified alpha
	void setAlpha(int n, double value);
	// Multiply value of specified alpha
	void multiplyAlpha(int n, double factor);
	// Returnu specified source parameter
	Parameter *sourceParameter(int n);
	// Return pointers to original parameters
	Parameter **sourceParameters();
	// Reset all fit flags to TRUE
	void fitAll();
	// Set individual fit flag
	void setFit(int n, bool status);
	// Return fit flag for parameter specified
	bool fit(int n);
	// Poke current values into associated Parameters
	bool poke();
	// Retrieve current values from associated Parameters
	bool peek();
	// Print vertex alpha contents
	void print();
	// Return new alpha array containing only those parameters marked as fittable
	Alpha fittable();


	/*
	// Operators / Functions
	*/
	public:
	double operator[](int index);
	void zero();
	void operator=(const Alpha& v);
	void operator+=(const Alpha&);
	void operator-=(const Alpha&);
	void operator/=(const double);
	Alpha operator+(const Alpha&) const;
	Alpha operator-(const Alpha&) const;
	Alpha operator*(const double) const;
};

#endif
