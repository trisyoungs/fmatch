/*
	*** Optimisation Strategy Definition
	*** strategy.cpp
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

#include "constants.h"
#include "strategy.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

/*
// Step Parameter
*/

// Constructor
StepParameter::StepParameter()
{
	// Private variables
	valueI_ = 0;
	valueD_ = 0.0;
	isDouble_ = FALSE;
	
	// Public variables
	prev = NULL;
	next = NULL;
}

// Destructor
StepParameter::~StepParameter()
{
}

// Set data (integer)
void StepParameter::set(int i)
{
	valueI_ = i;
	isDouble_ = FALSE;
}

// Set data (double)
void StepParameter::set(double d)
{
	valueD_ = d;
	isDouble_ = TRUE;
}

// Return as integer
int StepParameter::asInteger()
{
	return (isDouble_ ? int(valueD_) : valueI_);
}

// Return as double
double StepParameter::asDouble()
{
	return (isDouble_ ? valueD_ : double(valueI_));
}

/*
// Strategy Step
*/

// Constructor
StrategyStep::StrategyStep(StrategyStep::StepFunction sf)
{
	// Private variables
	function_ = sf;
	
	// Public variables
	prev = NULL;
	next = NULL;
}

// Destructor
StrategyStep::~StrategyStep()
{
}

const char *StepFunctionNames[StrategyStep::nStepFunctions] = { "[CLEAR]", "_ENDLOOP", "[GENES]", "[GRID]", "[LOOP] ", "[PRDCT]", "[PRINT]", "[RAND] ", "[RATTL]", "[SD]   ", "[SELCT]", "[HOT]  ", "[HOT2] ", "[COLD] ", "[STORE]", "[SUPER]", "[?????]" };
const char *StepFunctionKeywords[StrategyStep::nStepFunctions] = { "clear", "endloop", "genetic", "grid", "loop", "predict", "print", "randomise", "rattle", "sd", "select", "anneal", "anneal2", "downhill", "store", "super", "HIDDENFEATURE" };
StrategyStep::StepFunction StrategyStep::stepFunction(const char *s)
{
	for (int ik=0; ik < StrategyStep::nStepFunctions; ++ik) if (strcmp(s, StepFunctionKeywords[ik]) == 0) return (StrategyStep::StepFunction) ik;
	return StrategyStep::nStepFunctions;
}
const char *StrategyStep::stepFunction(StrategyStep::StepFunction sf)
{
	return StepFunctionKeywords[sf];
}

// Print step-tagged output
void StrategyStep::msg(const char *fmt, ...)
{
	static char msgs[8096];
	va_list arguments;
	msgs[0] = '\0';
	// Parse the argument list (...) and internally write the output string into msgs[]
	va_start(arguments,fmt);
	vsprintf(msgs,fmt,arguments);
	va_end(arguments);
	printf("%s\t%s", StepFunctionNames[function_], msgs);
}

// Print indented step-tagged output
void StrategyStep::msg(int indent, const char *fmt, ...)
{
	static char msgs[8096], tabs[128];
	va_list arguments;
	msgs[0] = '\0';
	// Construct tab list
	int n;
	for (n=0; n<indent; ++n) tabs[n] = '\t';
	tabs[n] = '\0';
	
	// Parse the argument list (...) and internally write the output string into msgs[]
	va_start(arguments,fmt);
	vsprintf(msgs,fmt,arguments);
	va_end(arguments);
	printf("%s%s\t%s", tabs, StepFunctionNames[function_], msgs);
}

// Add parameter to step (integer)
void StrategyStep::addParameter(int i)
{
	StepParameter *sp = parameters_.add();
	sp->set(i);
}

// Add parameter to step (double)
void StrategyStep::addParameter(double d)
{
	StepParameter *sp = parameters_.add();
	sp->set(d);
}
