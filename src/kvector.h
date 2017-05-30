/*
	*** KVector Definition
	*** kvector.h
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

#ifndef FM3_KVECTOR_H
#define FM3_KVECTOR_H

#include "vector3.h"

// KVector Definition
class KVector
{
	/*
	// Basic Data
	*/
	private:
	// Integer multiples of reciprocal vectors
	int x_, y_, z_;
	// Reciprocal-space kvector
	Vec3<double> kVector_;
	// Squared magnitude of kVector
	double magSq_;
	// Scaling factors for contribution from this kvector
	double factor_;

	public:
	// Set information for atom pair
	void set(int kx, int ky, int kz, Vec3< double > kvec, double magsq, double factor);
	// Return integer multiple along X-vector
	int x();
	// Return integer multiple along Y-vector
	int y();
	// Return integer multiple along Z-vector
	int z();
	// Return reciprocal-space k-vvector
	Vec3<double> vector();
	// Return squared magnitude of k-vector
	double magSq();
	// Return scaling factor for this kvector
	double factor();
};

#endif
