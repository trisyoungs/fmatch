/*
	*** KVector Definition
	*** kvector.cpp
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

#include "kvector.h"

// Set information for atom pair
void KVector::set(int kx, int ky, int kz, Vec3<double> kvec, double magsq, double factor)
{
	x_ = kx;
	y_ = ky;
	z_ = kz;
	kVector_ = kvec;
	magSq_ = magsq;
	factor_ = factor;
}

// Return integer multiple along X-vector
int KVector::x()
{
	return x_;
}

// Return integer multiple along Y-vector
int KVector::y()
{
	return y_;
}

// Return integer multiple along Z-vector
int KVector::z()
{
	return z_;
}

// Return reciprocal-space k-vvector
Vec3<double> KVector::vector()
{
	return kVector_;
}

// Return squared magnitude of k-vector
double KVector::magSq()
{
	return magSq_;
}

// Return scaling factor for this kvector
double KVector::factor()
{
	return factor_;
}
