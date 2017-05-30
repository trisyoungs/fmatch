/*
	*** Read in frames for custom FMatch3 input
	*** io_fm.cpp
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

// Read in frames from FMatch3 input file specified
bool TargetSystem::readConfigurationsFM3(const char* filename, Matrix &cell, double forceMultiplier, int cellType)
{
	LineParser parser(filename);
	
	if (!parser.isFileGood())
	{
		printf("Error opening configurations file '%s'.\n", filename);
		return FALSE;
	}
	
	// Continue reading from file until end
	// Input is keyword driven
	Configuration *cfg;
	Matrix unitCell = cell;
	double temp;
	Vec3<double> u, v;
	int n, natoms, nread = 0;
	do
	{
		// Read keyword line (hopefully)
		if (parser.getArgsDelim() != 0) break;
		TargetSystem::FMatch3Keyword keywd = TargetSystem::fmatch3Keyword(parser.argc(0));
		switch (keywd)
		{
			// Specify the start of a new configuration
			case (TargetSystem::ConfigurationKeyword):
				cfg = configurations_.add();
				cfg->setName(parser.argc(1));
				cfg->setCell(unitCell, cellType);
				++nread;
				break;
			case (TargetSystem::AtomsKeyword):
				natoms = parser.argi(1);
				if (!cfg->initialiseAtomArrays(natoms)) return FALSE;
				for (n=0; n<natoms; ++n)
				{
					if (parser.getArgsDelim() != 0)
					{
						printf("Error reading atom %i from file for configuration '%s'\n", n+1, cfg->name());
						return FALSE;
					}
					// Check correct number of arguments was given
					if (parser.nArgs() < 4)
					{
					}
					else if (parser.nArgs() > 4) printf("Warning - extra data found for atom %i\n", n+1);
					if (strcmp(parser.argc(0),atomTypeName(n)) != 0) printf("Warning - Atom name %i (%s) does not match that inferred from the species definition (%s).\n", n+1, parser.argc(0), atomTypeName(n));
					// Set atom data
					v.set(parser.argd(1), parser.argd(2), parser.argd(3));
					cfg->setAtomCoordinates(n, v);
				}
				break;
			case (TargetSystem::AtomsAndForcesKeyword):
				natoms = parser.argi(1);
				if (!cfg->initialiseAtomArrays(natoms)) return FALSE;
				for (n=0; n<natoms; ++n)
				{
					if (parser.getArgsDelim() != 0)
					{
						printf("Error reading atom and forces %i from file for configuration '%s'\n", n+1, cfg->name());
						return FALSE;
					}
					// Check correct number of arguments was given
					if (parser.nArgs() < 7)
					{
					}
					else if (parser.nArgs() > 8) printf("Warning - extra data found after atom and forces %i\n", n+1);
					if (strcmp(parser.argc(0),atomTypeName(n)) != 0) printf("Warning - Atom name %i (%s) does not match that inferred from the species definition (%s).\n", n+1, parser.argc(0), atomTypeName(n));
					// Set atom data
					v.set(parser.argd(1), parser.argd(2), parser.argd(3));
					u.set(parser.argd(4), parser.argd(5), parser.argd(6));
					u *= forceMultiplier;
					cfg->setAtomData(n, v, u, parser.hasArg(4) ? parser.argd(4) : 1.0);
				}
			case (TargetSystem::CellKeyword):
				// Get cell type
				if (strcmp(parser.argc(1),"cubic"))
				{
					if (parser.nArgs() != 3)
					{
						printf("Error: 'cubic' cell requires the side length to be specified (and nothing else).\n");
						return FALSE;
					}
					unitCell.setIdentity();
					unitCell[0] = parser.argd(2);
					unitCell[5] = unitCell[0];
					unitCell[10] = unitCell[0];
					cfg->setCell(unitCell, 1);
				}
				else if (strcmp(parser.argc(1),"orthorhombic"))
				{
					if (parser.nArgs() != 5)
					{
						printf("Error: 'orthorhombic' cell requires the three side lengths to be specified (and nothing else).\n");
						return FALSE;
					}
					unitCell.setIdentity();
					unitCell[0] = parser.argd(2);
					unitCell[5] = parser.argd(3);
					unitCell[10] = parser.argd(4);
					cfg->setCell(unitCell, 1);
				}
				else if (strcmp(parser.argc(1),"parallelepiped"))
				{
					if (parser.nArgs() != 8)
					{
						printf("Error: 'parallelepiped' cell requires the three side lengths and angles to be specified (and nothing else).\n");
						return FALSE;
					}
					// Work in unit vectors. Assume that A lays along x-axis
					unitCell.setColumn(0,1.0,0.0,0.0,0.0);
					// Assume that B lays in the xy plane. Since A={1,0,0}, cos(gamma) equals 'x' of the B vector.
					temp = cos(parser.argd(7)/DEGRAD);
					unitCell.setColumn(1,temp,sqrt(1.0 - temp*temp),0.0,0.0);
					// The C vector can now be determined in parts.
					// It's x-component is equal to cos(beta) since {1,0,0}{x,y,z} = {1}{x} = cos(beta)
					unitCell.setColumn(2,cos(parser.argd(6)/DEGRAD),0.0,0.0,0.0);
					// The y-component can be determined by completing the dot product between the B and C vectors
					cell[9] = ( cos(parser.argd(5)/DEGRAD) - cell[4]*cell[8] ) / cell[5];
					// The z-component is simply the remainder of the unit vector...
					cell[10] = sqrt(1.0 - cell[8]*cell[8] - cell[9]*cell[9]);
					// Lastly, adjust these unit vectors to give the proper cell lengths
					unitCell.columnMultiply(0,parser.argd(2));
					unitCell.columnMultiply(1,parser.argd(3));
					unitCell.columnMultiply(2,parser.argd(4));
					unitCell.setColumn(3, 0.0, 0.0, 0.0, 1.0);
					cfg->setCell(unitCell, 3);
				}
				else
				{
					printf("Error: Unrecognised celltype '%s'.\n", parser.argc(0));
					return FALSE;
				}
				break;
// 			case (TargetSystem::EnergyKeyword):
// 			case (TargetSystem::EnergyPerAtomKeyword):
// 			case (TargetSystem::ESPKeyword):
			case (TargetSystem::ForcesKeyword):
				natoms = parser.argi(1);
				if (!cfg->initialiseAtomArrays(natoms)) return FALSE;
				for (n=0; n<natoms; ++n)
				{
					if (parser.getArgsDelim() != 0)
					{
						printf("Error reading atom forces  %i from file for configuration '%s'\n", n+1, cfg->name());
						return FALSE;
					}
					// Check correct number of arguments was given
					if (parser.nArgs() < 4)
					{
					}
					else if (parser.nArgs() > 5) printf("Warning - extra data found for after atom forces %i\n", n+1);
					if (strcmp(parser.argc(0),atomTypeName(n)) != 0) printf("Warning - Atom name %i (%s) does not match that inferred from the species definition (%s).\n", n+1, parser.argc(0), atomTypeName(n));
					// Set atom data
					v.set(parser.argd(1), parser.argd(2), parser.argd(3));
					v *= forceMultiplier;
					cfg->setAtomForces(n, v, parser.hasArg(4) ? parser.argd(4) : 1.0);
				}
			default:
				printf("Internal Error: Action for keyword is undefined.\n");
				return FALSE;
		}
	} while (!parser.eofOrBlank());
	parser.closeFile();
	printf("Read in %i configurations from file '%s'.\n", nread, filename);
	return TRUE;
}
