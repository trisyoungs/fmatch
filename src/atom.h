/*
	*** Atom Definition
	*** atom.h
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

#ifndef FM3_ATOM_H
#define FM3_ATOM_H

#include "dnchar.h"

// Forward Declarations
class Parameter;

// Atom Definition
class Atom
{
	public:
	// Constructor / Destructor
	Atom();
	~Atom();
	// List pointers
	Atom *prev, *next;


	/*
	// Basic Data
	*/
	private:
	// Name (type) of atom
	Dnchar name_;
	// Integer type ID
	int typeId_;
	// Pointer to parameter containing charge for this atom
	Parameter *chargeParameter_;
	// Atomic mass (used in MD runs only)
	double mass_;
	
	public:
	// Set name of atom
	void setName(const char *s);
	// Return name of atom
	const char *name();
	// Set type id of atom
	void setTypeId(int id);
	// Return type ID
	int typeId();
	// Set charge parameter pointer
	void setChargeParameter(Parameter *param);
	// Return charge parameter
	Parameter *chargeParameter();
	// Return charge parameter valuel
	double charge();
	// Set mass of atom
	void setMass(double d);
	// Return mass of atom
	double mass();
};

#endif
