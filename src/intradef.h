/*
	*** Intramolecular Interaction Definition
	*** intradef.h
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

#ifndef FM3_INTRADEFINITION_H
#define FM3_INTRADEFINITION_H

#define MAXATOMS 4

#include "interaction.h"

// Intramolecular interaction definition
class Intra
{
	public:
	// Constructor / Destructor
	Intra();
	~Intra();
	// List pointers
	Intra *prev, *next;


	/*
	// Basic Data
	*/
	private:
	// Atom indices involved
	int atoms_[MAXATOMS];
	// Scale factors (for torsion interactions)
	double elecScale_, vdwScale_;
	// Pointer to interaction describing term
	Interaction *interaction_;

	public:
	// Set involved atoms
	void setAtoms(int i, int j, int k = -1, int l = -1);
	// Set scale factors (if used)
	void setScaleFactors(double escale, double vscale);
	// Return specified atom index
	int atom(int index);
	// Return atom i
	int i();
	// Return atom j
	int j();
	// Return atom k
	int k();
	// Return atom l
	int l();
	// Return electrostatic scale factor
	double elecScale();
	// Return van der Waals scale factor
	double vdwScale();
	// Set pointer to interaction definition
	void setInteraction(Interaction *i);
	// Return interaction pointer
	Interaction *interaction();
};

#endif
