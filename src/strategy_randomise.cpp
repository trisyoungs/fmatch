/*
	*** Parameter Randomisation Step
	*** strategy_randomise.cpp
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
#include "random.h"

// Constructor
RandomisationStep::RandomisationStep() : StrategyStep(StrategyStep::Randomise)
{
}

// Print step information
void RandomisationStep::print(int indent)
{
	msg(indent, "Free parameters will be randomised by +/- %f%%\n", parameters_[0]->asDouble()*100);
}

// Execute step
Alpha RandomisationStep::execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem)
{
	// Randomise starting parameters (if randomFactor_ != 0.0)
	double factor = parameters_[0]->asDouble();
	
// 	msg("Randomising input forcefield by +/- %f%%...\n", factor * 100.0);
	for (int n=0; n<inputAlpha.nAlpha(); ++n)
	{
		if (inputAlpha.fit(n)) inputAlpha.multiplyValue(n, 1.0 + (2.0*factor*Random::number() - factor));
	}
	inputAlpha.poke();
	ff->calculateDependencies();
	msg("Free parameters after randomisation:\n");
	ff->printFreeParameters();
	if (ff->nDependentParameters() > 0)
	{
		msg("Dependent parameters after randomisation:\n");
		ff->printDependentParameters();
	}
	targetSystem->updateBestAlpha(&inputAlpha);
	return inputAlpha;
}
