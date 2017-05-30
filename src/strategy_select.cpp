/*
	*** Select Step
	*** strategy_select.cpp
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
#include <string.h>

// Constructor
SelectStep::SelectStep() : StrategyStep(StrategyStep::Select)
{
}

// Print step information
void SelectStep::print(int indent)
{
	Dnchar s;
	for (StepParameter *p = parameters_.first(); p != NULL; p = p->next)
	{
		if (p->asInteger() < 0) s.strcatf("%i ", abs(p->asInteger()));
		else s.strcatf("%s ", Parameter::parameterStrength( (Parameter::ParameterStrength) p->asInteger() ));
		s.strcat("  ");
	}
	msg(indent, "Select parameter (sub)set for optimisation: %s\n", s.get());
}

// Execute step
Alpha SelectStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int i;
	
	// Get total number of all alpha
	int nAlpha = inputAlpha.nAlpha();
	bool optimiseParameter[nAlpha];
	bool optimiseType[Parameter::nParameterStrengths];
	for (i=0; i<nAlpha; ++i) optimiseParameter[i] = FALSE;
	for (i=0; i<Parameter::nParameterStrengths; ++i) optimiseType[i] = FALSE;
	
	// Set optimisation flags
	for (StepParameter *p = parameters_.first(); p != NULL; p = p->next)
	{
		if (p->asInteger() < 0)
		{
			i = abs(p->asInteger());
			// Check parameter ID is within range....
			if (i <= nAlpha) optimiseParameter[i-1] = TRUE;
			else
			{
				msg("Error: Parameter index %i is out of range (there are %i adjustable parameters in total)\n", i, nAlpha);
				targetSystem->flagAbort();
				return inputAlpha;
			}
		}
		else optimiseType[p->asInteger()] = TRUE;
	}

	// Count parameters while setting flags
	int count = 0;
	Parameter **sourceParams = inputAlpha.sourceParameters();
	for (i = 0; i<nAlpha; ++i)
	{
		if (optimiseType[Parameter::Any] || optimiseParameter[i] || optimiseType[sourceParams[i]->strength()])
		{
			++count;
			inputAlpha.setFit(i, TRUE);
		}
		else inputAlpha.setFit(i, FALSE);
	}
	if (count == 0)
	{
		msg("Error: Selection gives zero optimisable parameters.\n");
		targetSystem->flagAbort();
		return inputAlpha;
	}
	msg("Selected %i parameters.\n", count);
	
	return inputAlpha;
}
