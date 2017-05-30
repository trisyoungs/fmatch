/*
	*** XYZ file Input
	*** io_xyz.cpp
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

#include "system.h"
#include "lineparser.h"
#include <string.h>

// Read in frames from XYZ file specified
bool TargetSystem::readConfigurationsXYZ(const char* filename, Matrix &cell, double forceMultiplier, int cellType, bool hasVelocities)
{
	LineParser parser(filename);
	
	if (!parser.isFileGood())
	{
		printf("Error opening configurations file '%s'.\n", filename);
		return FALSE;
	}
	
	// Continue reading from file until end
	int nread = 0;
	do
	{
		// First line is number of atoms
		if (parser.getArgsDelim() != 0) break;
		int natoms = parser.argi(0);
		if (natoms < 1)
		{
			printf("Error: XYZ file contains a frame (no. %i) with natoms given as %i.\n", nread+1, natoms); 
		}

		// Second line is title
		parser.readNextLine();
		Configuration *cfg = configurations_.add();
		cfg->setName(parser.line());
		cfg->initialiseAtomArrays(natoms);
		cfg->setCell(cell, cellType);

		// Next lines are atoms, coordinates and forces
		int success = 0;
		for (int n=0; n<natoms; ++n)
		{
			success = parser.getArgsDelim();
			if (success == 1)
			{
				printf("Error: Corrupt atom data (%i) in configuration %i in file '%s'.\n", n+1, nread+1, filename);
				return FALSE;
			}
			if (parser.nArgs() != (hasVelocities ? 10 : 7))
			{
				printf("Error: Atom data %i in configuration %i in file '%s' is missing force data or contains extraneous data.\n", n+1, nread+1, filename);
				return FALSE;
			}
			Vec3<double> r, f;
			r.set(parser.argd(1), parser.argd(2), parser.argd(3));
			if (hasVelocities) f.set(parser.argd(7), parser.argd(8), parser.argd(9));
			else f.set(parser.argd(4), parser.argd(5), parser.argd(6));
			if (!cfg->setAtomData(n, r, f*forceMultiplier, 1.0)) return FALSE;
// 			if ((success == -1) && (n == natoms-1)) break;
		}
		++nread;
	} while (!parser.eofOrBlank());
	parser.closeFile();
	printf("Read in %i configurations from file '%s'.\n", nread, filename);
	return TRUE;
}
