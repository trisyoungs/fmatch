/*
	*** Loop Step
	*** strategy_loop.cpp
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
#include "ff.h"
#include "system.h"

// Constructor
LoopStep::LoopStep() : StrategyStep(StrategyStep::Loop)
{
}

// Loop step information
void LoopStep::print(int indent)
{
	int nLoops = parameters_[0]->asInteger();
	msg(indent, "Loop for %i iterations:\n", nLoops);
	// Now print out loop data
	for (StrategyStep *step = strategy_.first(); step != NULL; step = step->next) step->print(indent+2);
}

// Execute step
Alpha LoopStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int nLoops = parameters_[0]->asInteger();
	
	Alpha currentAlpha = inputAlpha;
	
	for (int n=0; n<nLoops; ++n)
	{
		msg("Iteration %i of %i\n", n+1, nLoops);
		
		// Now loop over strategy steps
		for (StrategyStep *step = strategy_.first(); step != NULL; step = step->next)
		{
			// Make sure currentAlpha are in the forcefield
// 			currentAlpha.poke();

			// Print out next step information and execute
			step->print();
			currentAlpha = step->execute(currentAlpha, ff, targetSystem);

			// Abort?
			if (targetSystem->abort()) break;
			
			// Threshold reached?
			if (targetSystem->bestAlpha().cost(targetSystem, ff) < targetSystem->threshold()) break;
			
		}

		// Abort?
		if (targetSystem->abort()) break;
		
		// Threshold reached?
		if (targetSystem->bestAlpha().cost(targetSystem, ff) < targetSystem->threshold()) break;

	}
	return currentAlpha;
}

// Return strategy list
List<StrategyStep> *LoopStep::strategy()
{
	return &strategy_;
}