/*
	*** Van der Waals Definition / Calculation
	*** vdw.cpp
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

#include "vdw.h"

// Constructor
VdwInteraction::VdwInteraction()
{
	// Private variables
	
	// Public variables
	next = NULL;
	prev = NULL;
}

// Destructor
VdwInteraction::~VdwInteraction()
{
}


