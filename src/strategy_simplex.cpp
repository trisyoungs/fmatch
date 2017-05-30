/*
	*** Simplex Minimisation Step
	*** strategy_simplex.cpp
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
#include "simplex.h"
#include "constants.h"
#include "ff.h"
#include "system.h"
#include "random.h"

/*
// Anneal
*/

// Constructor
SimplexAnnealStep::SimplexAnnealStep() : StrategyStep(StrategyStep::SimplexAnneal)
{
}

// Print step information
void SimplexAnnealStep::print(int indent)
{
	int maxLoops = parameters_[0]->asInteger();
	int nMoves = parameters_[1]->asInteger();
	double convergence = parameters_[2]->asDouble();
	double steplength = parameters_[3]->asDouble();
	msg(indent, "Simplex Simulated Annealing (max %i loops of %i moves, convergence = %12.5e, steplength = %f)\n", maxLoops, nMoves, convergence, steplength);
}

// Execute step
Alpha SimplexAnnealStep::execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem)
{
	int maxLoops = parameters_[0]->asInteger();
	int nMoves = parameters_[1]->asInteger();
	double convergence = parameters_[2]->asDouble();
	double steplength = parameters_[3]->asDouble();

	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();
	
	// Calculate and print initial cost
	msg("Initial cost = %11.5e\n", optAlpha.cost(targetSystem, ff));

	// Begin master Simplex Loop
	int loop = 1;
	bool betterPointFound = TRUE;
	Alpha newAlpha, bestAlpha = optAlpha;

	// Set temperature for this cycle
// 	double tStart = bestAlpha.cost(targetSystem, ff)*10.0, t;
	double tStart = 0.1, t;
	t = tStart;

	do
	{
		msg("Loop %5i of %5i\n", loop, maxLoops);

		msg("Annealing temperature for this loop is %12.5e\n", t);
		
		// Equalise length scales?
		if (targetSystem->equalise()) bestAlpha.calculateLengthScales(targetSystem, ff);
		else bestAlpha.resetLengthScales();
		
		// Initialise Simplex and determine temperature
		Simplex simplex(targetSystem, ff);
		simplex.initialise(bestAlpha, steplength);
		msg("Initial Simplex:\n");
		simplex.printVertexInformation();

		// Minimise
		newAlpha = simplex.minimise(nMoves, convergence, t);
		betterPointFound = simplex.betterPointFound();

		// Print loop information
		msg("Loop %i summary:\n", loop);
		
		// Find current extreme (cold) cost values
		if (newAlpha.cost(targetSystem, ff) < bestAlpha.cost(targetSystem, ff)) bestAlpha = newAlpha;
// 		msg("\tBest vertex in Simplex has cost of   %14.6e\n", newAlpha.cost(targetSystem, ff));
		msg("Best vertex found so far has cost of %14.6e\n", bestAlpha.cost(targetSystem, ff));
		
		// Print move information
		simplex.printMoveInformation();
		
		// Has the simplex converged on a set of parameters
		if (simplex.converged()) msg("Simplex converged to within threshold of %12.5e.\n\n", convergence);
		else msg("Simplex did not converge in %i moves - MSD = %12.5e, threshold = %12.5e.\n\n", nMoves, simplex.costMSD(), convergence);

		// Increase loop and decrease temperature
		++loop;
		t -= tStart/maxLoops;

	} while (loop <= maxLoops);

	// Poke best alpha back into travelling parameters
	if (bestAlpha.cost(targetSystem, ff) < inputAlpha.cost(targetSystem, ff))
	{
		bestAlpha.poke();
		inputAlpha.peek();
	}
	else msg("Simplex annealing failed to find better point than the current one.\n");

	return inputAlpha;
}

/*
// Anneal
*/

// Constructor
SimplexAnneal2Step::SimplexAnneal2Step() : StrategyStep(StrategyStep::SimplexAnneal2)
{
}

// Print step information
void SimplexAnneal2Step::print(int indent)
{
	int nCycles = parameters_[0]->asInteger();
	int movesPerCycle = parameters_[1]->asInteger();
	double tMultiplier = parameters_[2]->asDouble();
	double vertexTemperature = parameters_[3]->asDouble();
	double stepLength = parameters_[4]->asDouble();
	msg(indent, "Simplex Annealing method 2 (%i cycles of %i moves, tMult = %f, vTemp = %f, stepLength = %f)\n", nCycles, movesPerCycle, tMultiplier, vertexTemperature, stepLength);
}

// Execute step
Alpha SimplexAnneal2Step::execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem)
{
	int nCycles = parameters_[0]->asInteger();
	int movesPerCycle = parameters_[1]->asInteger();
	double tMultiplier = parameters_[2]->asDouble();
	double vertexTemperature = parameters_[3]->asDouble();
	double stepLength = parameters_[4]->asDouble();

	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();
	
	// Calculate and print initial cost
	msg("Initial cost = %11.5e\n", optAlpha.cost(targetSystem, ff));

	// Set temperature for this cycle
	double tStart = optAlpha.cost(targetSystem, ff)*tMultiplier;

	msg("Starting Simplex temperature is %12.5e\n", tStart);
	msg("Starting Vertex temperature is  %12.5e\n", vertexTemperature);
		
	// Equalise length scales?
	if (targetSystem->equalise()) optAlpha.calculateLengthScales(targetSystem, ff);
	else optAlpha.resetLengthScales();

	msg("***********************************\n");
	msg("Generating Initial Simplex Vertices\n");
	msg("***********************************\n");
	
	// Initialise Simplex
	Simplex simplex(targetSystem, ff);
	simplex.initialise(optAlpha, stepLength);

	msg("Initial Simplex:\n");
	simplex.printVertexInformation();

	// Minimise
	msg("****************************\n");
	msg("Performing Simplex Annealing\n");
	msg("****************************\n");
	Alpha newAlpha = simplex.minimiseTGAY(nCycles, movesPerCycle, 0.1, tStart, vertexTemperature);

	// Print move information
	simplex.printMoveInformation();

	// Poke best alpha back into travelling parameters
	if (newAlpha.cost(targetSystem, ff) < inputAlpha.cost(targetSystem, ff))
	{
		newAlpha.poke();
		inputAlpha.peek();
	}
	else msg("Simplex annealing2 failed to find better point than the current one.\n");

	return inputAlpha;
}

/*
// Minimise
*/

// Constructor
SimplexMinimiseStep::SimplexMinimiseStep() : StrategyStep(StrategyStep::SimplexMinimise)
{
}

// Print step information
void SimplexMinimiseStep::print(int indent)
{
	int maxLoops = parameters_[0]->asInteger();
	int nMoves = parameters_[1]->asInteger();
	double convergence = parameters_[2]->asDouble();
	double steplength = parameters_[3]->asDouble();
	msg(indent, "Downhill Simplex minimisation (max %i loops of %i moves, convergence = %12.5e, steplength = %f)\n", maxLoops, nMoves, convergence, steplength);
}

// Execute step
Alpha SimplexMinimiseStep::execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem)
{
	int maxLoops = parameters_[0]->asInteger();
	int nMoves = parameters_[1]->asInteger();
	double convergence = parameters_[2]->asDouble();
	double steplength = parameters_[3]->asDouble();

	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();

	// Calculate and print initial cost
	msg("Initial cost = %11.5e\n", optAlpha.cost(targetSystem, ff));

	// Begin master Simplex Loop
	int loop = 1;
	Alpha newAlpha, bestAlpha = optAlpha;
	do
	{
// 		msg("*********************\n");
		msg(" Loop %5i of %5i\n", loop, maxLoops);
// 		msg("*********************\n");

		// Initialise Simplex
		Simplex simplex(targetSystem, ff);
		simplex.initialise(bestAlpha, steplength);
		msg("Initial Simplex:\n");
		simplex.printVertexInformation();

		// Minimise
		newAlpha = simplex.minimise(nMoves, convergence, 0.0);

		// Print loop information
		msg("Loop %i summary:\n", loop);
		
		// Find current extreme (cold) cost values
		if (newAlpha.cost(targetSystem, ff) < bestAlpha.cost(targetSystem, ff)) bestAlpha = newAlpha;
		msg("Best vertex found so far has cost of %14.6e\n", bestAlpha.cost(targetSystem, ff));
		
		// Print move information
		simplex.printMoveInformation();
		
		// Has the simplex converged on a set of parameters
		if (simplex.converged()) msg("Simplex converged to within threshold of %12.5e.\n\n", convergence);
		else msg("Simplex did not converge in %i moves - MSD = %12.5e, threshold = %12.5e.\n\n", nMoves, simplex.costMSD(), convergence);

		// Increase loop
		++loop;

	} while (loop <= maxLoops);

	// Poke best alpha back into travelling parameters
	if (bestAlpha.cost(targetSystem, ff) < inputAlpha.cost(targetSystem, ff))
	{
		bestAlpha.poke();
		inputAlpha.peek();
	}
	else msg("Simplex annealing2 failed to find better point than the current one.\n");

	return inputAlpha;
}


/*
// SuperSimplex
*/

// Constructor
SuperSimplexStep::SuperSimplexStep() : StrategyStep(StrategyStep::SuperSimplex)
{
}

// Print step information
void SuperSimplexStep::print(int indent)
{
	int nMoves = parameters_[0]->asInteger();
	double convergence = parameters_[1]->asDouble();
	msg(indent, "Initial minimisation via Super Simplex (nMoves = %i, convergence = %f)\n", nMoves, convergence);
}

// Execute step
Alpha SuperSimplexStep::execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem)
{
	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();
	
	// Find starting point with Super-Simplex
	msg("Performing Super-Simplex to search for better parameters...\n");
	bool foundBetter = FALSE;
	Alpha tempAlpha, newAlpha, bestAlpha = optAlpha;
	int nMoves = parameters_[0]->asInteger();
	double convergence = parameters_[1]->asDouble();

	// Cycle over each alpha, applying a normal nudge to each. Then, construct a Simplex around this point and zero-T minimise
	for (int n=0; n<nParams+1; ++n)
	{
		tempAlpha.set(optAlpha, n == 0 ? -1 : n-1, 2.0);
		Simplex simplex(targetSystem, ff);
		simplex.initialise(tempAlpha, 0.5);
		newAlpha = simplex.minimise(10, 0.01, 0.0);
		simplex.printMoveInformation();
		if (newAlpha.cost(targetSystem, ff) < bestAlpha.cost(targetSystem, ff))
		{
			bestAlpha = newAlpha;
			foundBetter = TRUE;
		}
	}
	if (foundBetter)
	{
		msg("SuperSimplex found better point with cost = %11.5e\n", bestAlpha.cost(targetSystem, ff));
		msg("Free parameters after SuperSimplex:\n");
		ff->printFreeParameters();
		if (ff->nDependentParameters() > 0)
		{
			msg("Dependent parameters after SuperSimplex:\n");
			ff->printDependentParameters();
		}
		bestAlpha.poke();
		inputAlpha.peek();
	}
	else msg("SuperSimplex could not find better point than the starting one.\n");

	return inputAlpha;
}
