/*
	*** AtomPair Definition
	*** atompair.h
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

#ifndef FM3_ATOMPAIR_H
#define FM3_ATOMPAIR_H

#include "vector3.h"

// Forward Declarations
class Parameter;
class Atom;

// AtomPair Definition
class AtomPair
{
	public:
	// Constructor / Destructor
	AtomPair();
	~AtomPair();


	/*
	// Basic Data
	*/
	private:
	// Global ids of involved atoms
	int i_, j_;
	// Pointers to involved Atoms
	Atom *atomI_, *atomJ_;
	// Minimum image vector between atom pair
	Vec3<double> vij_;
	// Minimum image distance between atom pair
	double rij_;
	// Scaling factors for interactions between these two atoms
	double elecScale_, vdwScale_;

	public:
	// Set information for atom pair
	void set(int i, int j, Atom *atomi, Atom *atomj, Vec3<double> vij, double rij, double escale, double vscale);
	// Return id of atom i
	int i();
	// Return id of atom j
	int j();
	// Return pointer to atom i
	Atom *atomI();
	// Return pointer to atom j
	Atom *atomJ();
	// Return type of atom i
	int typeI();
	// Return type of atom j
	int typeJ();
	// Return charge on atom i
	double chargeI();
	// Return charge on atom j
	double chargeJ();
	// Return minimum image vector between atoms
	Vec3<double> vij();
	// Return minimum image distance between atoms
	double rij();
	// Return electrostatic scaling factor for interactinos between these two atoms
	double elecScale();
	// Return van der Waals scaling factor for interactinos between these two atoms
	double vdwScale();
};

#endif
