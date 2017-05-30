/*
	*** Atom Definition
	*** atom.cpp
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

#include "atom.h"
#include "parameter.h"
#include <stdio.h>

// Constructor
Atom::Atom()
{
	// Private variables
	typeId_ = -1;
	chargeParameter_ = NULL;
	mass_ = -1.0;
	
	// Public variables
	prev = NULL;
	next = NULL;
}

// Destructor
Atom::~Atom()
{
}

// Set name of atom
void Atom::setName(const char *s)
{
	name_ = s;
}

// Return name of atom
const char *Atom::name()
{
	return name_.get();
}

// Set type id of atom
void Atom::setTypeId(int id)
{
	typeId_ = id;
}

// Return type ID
int Atom::typeId()
{
	return typeId_;
}

// Set charge parameter pointer
void Atom::setChargeParameter(Parameter *param)
{
	chargeParameter_ = param;
}

// Return charge parameter
Parameter *Atom::chargeParameter()
{
	return chargeParameter_;
}

// Return charge parameter valuel
double Atom::charge()
{
	if (chargeParameter_ == NULL) printf("BUG: Charge j requested, but pointer in Atom is NULL.n");
	else return chargeParameter_->value();
	return 0.0;
}

// Set mass of atom
void Atom::setMass(double d)
{
	mass_ = d;
}

// Return mass of atom
double Atom::mass()
{
	return mass_;
}
