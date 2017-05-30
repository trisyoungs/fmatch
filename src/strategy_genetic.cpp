/*
	*** Genetic Algorithm Step
	*** strategy_genetic.cpp
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
#include "genetic.h"
#include "system.h"

// Constructor
GeneticStep::GeneticStep() : StrategyStep(StrategyStep::GeneticMinimise)
{
}

// Print step information
void GeneticStep::print(int indent)
{
	int popSize = parameters_[0]->asInteger();
	int nGenerations = parameters_[1]->asInteger();
	double lengthScale = parameters_[2]->asDouble();
	double selectionThreshold = parameters_[3]->asDouble();
	double mutateRate = parameters_[4]->asDouble();
	msg("Genetic Algorithm Search (populationSize = %i, nGenerations = %i, baseLengthScale = %5.2f, selectionRatio = %5.2f, mutate% = %5.2f)\n", popSize, nGenerations, lengthScale, selectionThreshold, mutateRate);
}

// Execute step
Alpha GeneticStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int popSize = parameters_[0]->asInteger();
	int nGenerations = parameters_[1]->asInteger();
	double lengthScale = parameters_[2]->asDouble();
	double selectionThreshold = parameters_[3]->asDouble();
	double mutateRate = parameters_[4]->asDouble();

	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();
	
	// Initialise Genetic Algorithm
	if (targetSystem->equalise()) optAlpha.calculateLengthScales(targetSystem, ff);
	else optAlpha.resetLengthScales();
	Genetic primordialSoup(targetSystem, ff);
	primordialSoup.initialise(optAlpha, popSize, lengthScale);
	
	// Begin Genetic Algorithm
	Alpha bestAlpha = optAlpha, newAlpha;
	for (int cycle=0; cycle<nGenerations; ++cycle)
	{
		msg("Cycle %i\n",  cycle+1);
		newAlpha = primordialSoup.evolve(selectionThreshold, mutateRate);
		if (newAlpha.cost(targetSystem, ff) < bestAlpha.cost(targetSystem, ff)) bestAlpha = newAlpha;
		if (primordialSoup.converged(0.1, selectionThreshold))
		{
			msg("Gene pool has converged to a single set of parameters.\nPool will now bo refreshed...\n");
			primordialSoup.refresh();
		}
	}
	
	// Poke best alpha back into travelling parameters
	bestAlpha.poke();
	inputAlpha.peek();

	return bestAlpha;
}
