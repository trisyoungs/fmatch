/*
	*** Print Step
	*** strategy_print.cpp
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
#include "configuration.h"
#include "system.h"

// Constructor
PrintStep::PrintStep() : StrategyStep(StrategyStep::Print)
{
}

// Print step information
void PrintStep::print(int indent)
{
	int info = parameters_[0]->asInteger();
	
	if (info == 1) msg(indent, "Print summary of CURRENT parameter information.\n");
	else if (info == 2) msg(indent, "Print summary of BEST parameter information obtained so far.\n");
}

// Execute step
Alpha PrintStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int info = parameters_[0]->asInteger();

	if (info == 1)
	{
// 		// Put values back into forcefield, calculate dependencies and print
// 		msg("Best obtained parameters so far have cost of %14.6e:\n\n", inputAlpha.cost(targetSystem, ff));
// 		inputAlpha.poke();
// 		ff->calculateDependencies();
// 		msg("Free parameters:\n");
// 		ff->printFreeParameters();
// 		if (ff->nDependentParameters() > 0)
// 		{
// 			msg("Dependent parameters:\n");
// 			ff->printDependentParameters();
// 		}
		msg("Current (%11.5e) = ", inputAlpha.cost(targetSystem, ff));
		inputAlpha.print();
	}
	else if (info == 2)
	{
		msg("   Best (%11.5e) = ", targetSystem->bestAlpha().cost(targetSystem, ff));
		targetSystem->bestAlpha().print();
	}
	
	return inputAlpha;
}
