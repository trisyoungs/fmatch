/*
	*** Electrostatic Potential Point Definition
	*** esppoint.h
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

#ifndef FM3_ESPPOINT_H
#define FM3_ESPPOINT_H

#include "vector3.h"

// ESP Point Definition
class ESPPoint
{
	public:
	// Constructor
	ESPPoint();
	// List pointers
	ESPPoint *prev, *next;

	private:
	// Coordinates of point
	Vec3<double> r_;
	// Energy at point
	double energy_;
	// Weighting of point for SOSE
	double weight_;

	public:
	// Set data
	void set(Vec3<double> r, double energy, double weight = 1.0);
	// Return coordinates of point
	Vec3<double> r();
	// Return energy at point
	double energy();
	// Return weight of contribution to SOSE
	double weight();
};

#endif
