/*
	*** AtomPair Definition
	*** atompair.cpp
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

#include "atompair.h"
#include "atom.h"
#include <stdio.h>

// Constructor
AtomPair::AtomPair()
{
	// Private variables
	i_ = -1;
	j_ = -1;
	atomI_ = NULL;
	atomJ_ = NULL;
}

// Destructor
AtomPair::~AtomPair()
{
}

// Set information for atom pair
void AtomPair::set(int i, int j, Atom *atomi, Atom *atomj, Vec3<double> vij, double rij, double escale, double vscale)
{
	i_ = i;
	j_ = j;
	atomI_ = atomi;
	atomJ_ = atomj;
	vij_ = vij;
	rij_ = rij;
	elecScale_ = escale;
	vdwScale_ = vscale;
}

// Return id of atom i
int AtomPair::i()
{
	return i_;
}

// Return id of atom j
int AtomPair::j()
{
	return j_;
}

// Return pointer to atom i
Atom *AtomPair::atomI()
{
	return atomI_;
}

// Return pointer to atom j
Atom *AtomPair::atomJ()
{
	return atomJ_;
}

// Return type of atom i
int AtomPair::typeI()
{
	if (atomI_ == NULL) printf("BUG: Type of atom i requested, but atomI_ pointer is NULL.\n");
	else return atomI_->typeId();
	return -1;
}

// Return type of atom j
int AtomPair::typeJ()
{
	if (atomJ_ == NULL) printf("BUG: Type of atom j requested, but atomJ_ pointer is NULL.\n");
	else return atomJ_->typeId();
	return -1;
}

// Return charge on atom i
double AtomPair::chargeI()
{
	if (atomI_ == NULL) printf("BUG: Charge of atom i requested, but atomI_ pointer is NULL.\n");
	else return atomI_->charge();
	return 0.0;
}

// Return charge on atom j
double AtomPair::chargeJ()
{
	if (atomJ_ == NULL) printf("BUG: Charge of atom j requested, but atomJ_ pointer is NULL.\n");
	else return atomJ_->charge();
	return 0.0;
}

// Return minimum image vector between atoms
Vec3<double> AtomPair::vij()
{
	return vij_;
}

// Return minimum image distance between atoms
double AtomPair::rij()
{
	return rij_;
}

// Return electrostatic scaling factor for interactinos between these two atoms
double AtomPair::elecScale()
{
	return elecScale_;
}

// Return van der Waals scaling factor for interactinos between these two atoms
double AtomPair::vdwScale()
{
	return vdwScale_;
}
