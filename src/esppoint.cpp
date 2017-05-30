/*
	*** Electrostatic Potential Point Definition
	*** esppoint.cpp
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

#include "esppoint.h"

// Constructor
ESPPoint::ESPPoint()
{
	// Private variables
	weight_ = 0.0;
	energy_ = 0.0;

	// Public variables
	prev = NULL;
	next = NULL;
}

// Set data
void ESPPoint::set(Vec3<double> r, double energy, double weight)
{
	r_ = r;
	energy_ = energy;
	weight_ = weight;
}

// Return coordinates of point
Vec3<double> ESPPoint::r()
{
	return r_;
}

// Return energy at point
double ESPPoint::energy()
{
	return energy_;
}

// Return weight of contribution to SOSE
double ESPPoint::weight()
{
	return weight_;
}

