/*
	*** Intramolecular Interaction Definition
	*** intradef.cpp
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

#include "intradef.h"
#include <stdio.h>

// Constructor
Intra::Intra()
{
	// Private variables
	for (int n=0; n<MAXATOMS; ++n) atoms_[n] = -1;
	interaction_ = NULL;
	
	// Public variables
	prev = NULL;
	next = NULL;
}

// Destructor
Intra::~Intra()
{
}

// Set involved atoms
void Intra::setAtoms(int i, int j, int k, int l)
{
	atoms_[0] = i;
	atoms_[1] = j;
	atoms_[2] = k;
	atoms_[3] = l;
}

// Set scale factors (if used)
void Intra::setScaleFactors(double escale, double vscale)
{
	elecScale_ = escale;
	vdwScale_ = vscale;
}

// Return specified atom index
int Intra::atom(int index)
{
	// Check range
	if (index < 0 || index >= MAXATOMS)
	{
		printf("Error: Index out of range for integer retrieval\n");
		return atoms_[0];
	}
	return atoms_[index];
}

// Return atom i
int Intra::i()
{
	return atoms_[0];
}

// Return atom j
int Intra::j()
{
	return atoms_[1];
}

// Return atom k
int Intra::k()
{
	return atoms_[2];
}

// Return atom l
int Intra::l()
{
	return atoms_[3];
}

// Return electrostatic scale factor
double Intra::elecScale()
{
	return elecScale_;
}

// Return van der Waals scale factor
double Intra::vdwScale()
{
	return vdwScale_;
}

// Set pointer to interaction definition
void Intra::setInteraction(Interaction *i)
{
	interaction_ = i;
}

// Return interaction pointer
Interaction *Intra::interaction()
{
	return interaction_;
}
