/*
	*** Hidden Feature Step
	*** strategy_hidden.cpp
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
#include "simplex.h"
#include "system.h"
#include "ff.h"

// Constructor
HiddenFeatureStep::HiddenFeatureStep() : StrategyStep(StrategyStep::HiddenFeature)
{
}

// Print step information
void HiddenFeatureStep::print(int indent)
{
	msg("Hidden Feature\n");
}

// Execute step
Alpha HiddenFeatureStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int popSize = 100; //parameters_[0]->asInteger();
	int nGenerations = 100; //parameters_[1]->asInteger();
	double variation = 0.1; //parameters_[2]->asDouble();
	double selectionThreshold = 0.5; //parameters_[3]->asDouble();
	double mutateRate = 0.02; //parameters_[5]->asDouble();
	
	// Initialise Gene Population
	Genetic primordialSoup(targetSystem, ff);
	primordialSoup.initialise(inputAlpha, popSize, variation);
	Alpha gene, newAlpha = inputAlpha;
	
	// Begin Genetic Algorithm
	for (int cycle=0; cycle<nGenerations; ++cycle)
	{
		msg("Cycle %i\n",  cycle+1);

		// Perform simplex annealing / minimisation on each gene in turn
		for (int n=0; n<popSize; ++n)
		{
			msg("Gene %i of %i\n", n+1, popSize);
			gene = primordialSoup.gene(n);
			msg("Starting cost of gene is %11.5e\n", gene.cost(targetSystem, ff));
			//gene.poke();
			//ff->printFreeParameters();
			if (targetSystem->equalise()) gene.calculateLengthScales(targetSystem, ff);
			else gene.resetLengthScales();
			Simplex simplex(targetSystem, ff);
			simplex.initialise(gene, 0.1);
			gene = simplex.minimise(50, 0.01, gene.cost(targetSystem, ff));
			msg("Final cost of gene is %11.5e\n", gene.cost(targetSystem, ff));
// 			gene.poke();
// 			ff->printFreeParameters();
 			gene.resetLengthScales();
			Simplex minSimplex(targetSystem, ff);
			minSimplex.initialise(gene, 0.10);
			gene = minSimplex.minimise(50, 0.01, 0.0);

			primordialSoup.setGene(n, gene);
// 			msg("Initial Simplex:\n");
// 			simplex.printVertexInformation();
		}
		
		// Evolve gene population
		newAlpha = primordialSoup.evolve(selectionThreshold, mutateRate);
		
		if (primordialSoup.converged(0.1, selectionThreshold))
		{
			msg("Gene pool has converged to a single set of parameters.\n");
			break;
		}
/*		
		// Has pool become stagnant?  (i.e. converged on a set of solutions)
		if (primordialSoup.converged(0.001, deathRate))
		{
			msg("Pool has converged on a set of parameters, and so will now be refreshed.\n");
			primordialSoup.refresh();
		}*/
	}
	return newAlpha;
}
