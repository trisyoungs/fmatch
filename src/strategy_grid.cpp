/*
	*** Grid Search Minimise Step
	*** strategy_grid.cpp
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

#include "strategy.h"
#include "system.h"

// Constructor
GridMinimiseStep::GridMinimiseStep() : StrategyStep(StrategyStep::GridMinimise)
{
}

// Print step information
void GridMinimiseStep::print(int indent)
{
	int nSteps = parameters_[0]->asInteger();
	double stepSize = parameters_[1]->asDouble();
	msg(indent, "Grid Search (+/-steps per parameter = %i, fractional change per step = %5.2f)\n", nSteps, stepSize);
}

// Execute step
Alpha GridMinimiseStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int nSteps = parameters_[0]->asInteger();
	double stepSize = parameters_[1]->asDouble();
	int n, m, nCycles;
	
	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();

	nCycles = nSteps*2+1;
	for (n=1; n<nParams; ++n) nCycles *= nSteps*2+1;
	msg("Grid search over %i parameters with %i steps per parameter requires %i steps.\n", optAlpha.nAlpha(), nSteps*2+1, nCycles);

	// Setup arrays
	double deltas[nParams], newCost;
	int counts[nParams];
	for (n=0; n<nParams; ++n)
	{
		deltas[n] = optAlpha.alpha(n)*stepSize;
		counts[n] = -nSteps;
	}
	
	// Begin Grid Search
	Alpha bestAlpha, gridAlpha = optAlpha;
	bestAlpha = optAlpha;
	
	for (n=0; n<nCycles; ++n)
	{
// 		printf("LOOP INFO COUNTS ");
// 		for (m=0; m<nParams; ++m) printf("%4i ", counts[m]);
// 		printf("\n");
		// Change most quickly varying parameternSteps*2+1
		gridAlpha.setAlpha(0, optAlpha.alpha(0) + counts[0]*deltas[0]);
		
		// Evaluate cost of new parameters
		newCost = gridAlpha.cost(targetSystem, ff);
		if (newCost < bestAlpha.cost(targetSystem, ff)) bestAlpha = gridAlpha;
		targetSystem->updateBestAlpha(&bestAlpha);
		
		// Increase count parameter(s)
		++counts[0];
		if (counts[0] > nSteps)
		{
			// Have done a full cycle of parameter zero, so reset its counter and check the others
			counts[0] = -nSteps;
			for (m=1; m<nParams; ++m)
			{
				if (counts[m] == nSteps)
				{
					// Reset this parameter, and allow loop to move onto the next
					counts[m] = -nSteps;
					gridAlpha.setAlpha(m, optAlpha.alpha(m) + counts[m]*deltas[m]);
					
				}
				else
				{
					// Not done with this parameter's loop yet, so increase count and continue
					++counts[m];
					gridAlpha.setAlpha(m, optAlpha.alpha(m) + counts[m]*deltas[m]);
					break;
				}
			}
		}
	}

	// Poke best alpha back into travelling parameters
	bestAlpha.poke();
	inputAlpha.peek();

	return inputAlpha;
}
