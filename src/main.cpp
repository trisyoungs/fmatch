/*
	*** FMatch3 Main
	*** main.cpp
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

#include <iostream>
#include "system.h"

int main(int argc, char *argv[])
{
	// Variables
	Alpha testParameters;
	LineParser parser;
	int nSteps;
	double temperature, timestep;
	TargetSystem::JobType jobType = TargetSystem::Optimise;

	// Print GPL license information
	printf("FMatch3 version Alpha, Copyright (C) 2011 T. Youngs.\n");
	printf("FMatch3 comes with ABSOLUTELY NO WARRANTY.\n");
	printf("This is free software, and you are welcome to redistribute it under certain conditions.\n");
	printf("For more details read the GPL at <http://www.gnu.org/copyleft/gpl.html>.\n\n");

	// Create main structure definition
	TargetSystem targetSystem;

	// Parse command-line options
	int n=1;
	Dnchar inputFile;
	while (n < argc)
	{
		// Does argument begin with a '-'?
		if (argv[n][0] == '-') switch (argv[n][1])
		{
			// Dummy run - check setup, print strategy, and quit
			case ('d'):
				jobType = TargetSystem::DummyRun;
				++n;
				break;
			// Perform forces test
			case ('f'):
				jobType = TargetSystem::TestForces;
				++n;
				break;
			// Display help
			case ('h'):
				printf("Usage:  fmatch3 <inputfile> [options]\n\nCommand-line switches are:\n");
				printf("\t-d\n\t\tCheck setup, print job strategy, and quit\n");
				printf("\t-f\n\t\tCalculate force SOSE for each configuration stack and quit\n");
				printf("\t-h\n\t\tDisplay this message\n");
				printf("\t-m <nsteps> <temperature, K> <timestep, ps>\n\t\tRun molecular dynamics simulation\n");
				printf("\t-t\n\t\tTest supplied parameter values in context of current forcefield (use '*' to keep for original value)\n");
				return 0;
				break;
			// Run MD simulation
			case ('m'):
				jobType = TargetSystem::MolecularDynamics;
				nSteps = atoi(argv[++n]);
				temperature = atof(argv[++n]);
				timestep = atof(argv[++n]);
				targetSystem.md.set(nSteps, temperature, timestep);
				++n;
				break;
			// Test supplied parameters
			case ('t'):
				// Get parameters from supplied argument
				parser.getArgsDelim(argv[++n]);
				for (int i=0; i<parser.nArgs(); ++i)
				{
					// Check alpha is within range
					if (i < testParameters.nAlpha()) testParameters.setAlpha(i, parser.argd(i));
					else
					{
						printf("Error: Too many parameter values supplied - number of free parameters is %i.\n", testParameters.nAlpha());
						return -1;
					}
					jobType = TargetSystem::TestParameters;
					testParameters.poke();
				}
				++n;
				break;
			default:
				printf("Unrecognised command-line option '%c'. Run 'fmatch' without any arguments to see available options.\n", argv[n][1]);
				return -1;
		}
		else
		{
			// Must be the input filename, so open file and prepare system
			// Check that we haven't already set up
			if (targetSystem.setup())
			{
				printf("Error: System has already been set up. More than one input file supplied?\n");
				return -1;
			}
			if (!targetSystem.open(argv[n])) return -1;
			if (!targetSystem.setup()) return -1;
			// Grab alpha values here, in case some of the other CLI options need them
			testParameters = targetSystem.originalAlpha();
			++n;
		}
	}

	// Initialise random seed
	srand( targetSystem.seed() );

	// Perform the task requested
	int result = targetSystem.run(jobType);

	// Done.
	return result;
}

