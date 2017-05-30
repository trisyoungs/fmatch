/*
	*** Rattle Step
	*** strategy_rattle.cpp
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
#include "random.h"

// Constructor
RattleStep::RattleStep() : StrategyStep(StrategyStep::GridMinimise)
{
}

// Print step information
void RattleStep::print(int indent)
{
	int nSteps = parameters_[0]->asInteger();
	double maxVariance = parameters_[1]->asDouble();
	msg(indent, "Rattle parameters %i times by up to %5.2f\n", nSteps, maxVariance);
}

// Execute step
Alpha RattleStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int nSteps = parameters_[0]->asInteger();
	double maxVariance = parameters_[1]->asDouble();
	
	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();

	// Begin Rattle!
	Alpha bestAlpha = optAlpha, newAlpha = optAlpha;
	int n, step;
	double r, newCost;
	for (step=0; step<nSteps; ++step)
	{
// 		// Randomise parameters
// 		for (n=0; n<optAlpha.nAlpha(); ++n)
// 		{
// 			r = 1.0 + (2.0*Random::number() - 1.0)*maxVariance;
// 			newAlpha.setAlpha(n, bestAlpha.alpha(n)*r);
// 		}
		
		// TEST
		r = (2.0*Random::number() - 1.0)*maxVariance;
		newAlpha.setAlpha(0, bestAlpha.alpha(0)*(1.0+r));
		newAlpha.setAlpha(1, bestAlpha.alpha(1)*(1.0-r));

		// Evaluate cost of new parameters
		newCost = newAlpha.cost(targetSystem, ff);
		printf("Trial cost = %11.5e  ", newCost); newAlpha.print();
		if (newCost < bestAlpha.cost(targetSystem, ff))
		{
			bestAlpha = newAlpha;
			targetSystem->updateBestAlpha(&bestAlpha);
		}
	}

	// Poke best alpha back into travelling parameters
	bestAlpha.poke();
	inputAlpha.peek();

	return inputAlpha;
}
