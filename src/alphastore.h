/*
	*** Alpha Store
	*** alphastore.h
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

#ifndef FM3_ALPHASTORE_H
#define FM3_ALPHASTORE_H

#include "parameter.h"
#include "list.h"
#include "alpha.h"

// Forward Declarations
class Forcefield;
class TargetSystem;

// Alpha Statistics
class AlphaStore
{
	public:
	// Constructor / Destructor
	AlphaStore();
	~AlphaStore();
	
	private:
	// List of all sets of alpha stored, in ascending cost order
	List<Alpha> alpha_;
	// Size of arrays (number of alpha values)
	int nAlpha_;
	// Running minima, Maxima, Averages, and Standard Deviations for alpha
	double *minima_, *maxima_, *mean_, *sd_;
	
	public:
	// Clear all data from store
	void clear();
	// Add supplied alpha into store
	void add(Alpha *vertex, TargetSystem *targetSystem, Forcefield *ff);
	// Return number of sets of alpha currently in store
	int nData();
	// Print statistics
	void print();
	// Copy average values into supplied Alpha
	void copyAverages(Alpha *alpha);
	// Predict new set of alpha from stored vertex information
	Alpha predict(TargetSystem* targetSystem, Forcefield* ff);
};

#endif
