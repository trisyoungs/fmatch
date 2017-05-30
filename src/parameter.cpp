/*
	*** Parameter Definition
	*** parameter.cpp
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

#include "parameter.h"
#include "constants.h"
#include "interaction.h"
#include "dnchar.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// Constructor
Parameter::Parameter()
{
	// Private variables
	parent_ = NULL;
	type_ = Parameter::FreeParameter;
	strength_ = Parameter::Strong;
	lowerLimit_ = 0.0;
	upperLimit_ = 0.0;
	hasLowerLimit_ = FALSE;
	hasUpperLimit_ = FALSE;
	changeLog_ = 1;

	// Public variables
	next = 0;
	prev = 0;
}

// Destructor
Parameter::~Parameter()
{
}

const char *ParameterStrengthKeywords[Parameter::nParameterStrengths] = { "strong", "moderate", "weak", "all" };
Parameter::ParameterStrength Parameter::parameterStrength(const char *s)
{
	for (int ps=0; ps < Parameter::nParameterStrengths; ++ps) if (strcmp(s, ParameterStrengthKeywords[ps]) == 0) return (Parameter::ParameterStrength) ps;
	return Parameter::nParameterStrengths;
}
const char *Parameter::parameterStrength(Parameter::ParameterStrength ps)
{
	return ParameterStrengthKeywords[ps];
}

// Set name, value, and variability of this parameter
void Parameter::set(const char *s, double value)
{
	// Check first character of parameter name - if '@' then value is degrees and should be converted to radians
	if (s[0] == '@')
	{
		name_ = ++s;
		value_ = value/DEGRAD;
		printf("Note: Input value '%f' for parameter '%s' converted from degrees to radians as requested (now '%f').\n", value, name_.get(), value_);
	}
	else
	{
		name_ = s;
		value_ = value;
	}
}

// Return name of this parameter
const char *Parameter::name()
{
	return name_.get();
}

// Set parameter type
void Parameter::setType(ParameterType type, ParameterStrength strength)
{
	type_ = type;
	strength_ = strength;
}

// Return parameter type
Parameter::ParameterType Parameter::type()
{
	return type_;
}

// Return parameter strength
Parameter::ParameterStrength Parameter::strength()
{
	return strength_;
}

// Set current value of the parameter
bool Parameter::setValue(double d)
{
	if (type_ == Parameter::FixedParameter)
	{
		printf("BUG: This parameter ('%s') is fixed and cannot have its value set.\n", name_.get());
		return FALSE;
	}
	value_ = d;
	++changeLog_;
	return TRUE;
}

// Return current value of the parameter
double Parameter::value()
{
	return value_;
}

// Set parent interaction
void Parameter::setParent(Interaction *parent)
{
	parent_ = parent;
}

// Return parent interaction
Interaction *Parameter::parent()
{
	return parent_;
}

// Return parent interaction's functional form
const char *Parameter::parentForm()
{
	if (parent_ == NULL) printf("BUG: Interaction parent is not set in parameter '%s'\n", name_.get());
	else
	{
		if (parent_->form() == Interaction::NoForm) printf("BUG: Interaction parent has no defined form.\n");
		else return parent_->formText();
	}
	return "???";
}

// Return parent interaction's involved atoms
const char *Parameter::parentAtoms()
{
	if (parent_ == NULL) printf("BUG: Interaction parent is not set in parameter '%s'\n", name_.get());
	else return parent_->atomsText();
	return "???";
}

// Set lower limit for parameter
void Parameter::setLowerLimit(double minValue)
{
	lowerLimit_ = minValue;
	hasLowerLimit_ = TRUE;
}

// Set upper limit for parameter
void Parameter::setUpperLimit(double maxValue)
{
	upperLimit_ = maxValue;
	hasUpperLimit_ = TRUE;
}

// Return whether parameter has any associated limits
bool Parameter::hasLimits()
{
	return (hasLowerLimit_ || hasUpperLimit_);
}

// Return whether parameter has a lower limit
bool Parameter::hasLowerLimit()
{
	return hasLowerLimit_;
}

// Return whether parameter has an upper limit
bool Parameter::hasUpperLimit()
{
	return hasUpperLimit_;
}

// Return whether current parameter value is within accepted limits
bool Parameter::withinLimits()
{
	if (hasLowerLimit_ && (value_ < lowerLimit_)) return FALSE;
	if (hasUpperLimit_ && (value_ > upperLimit_)) return FALSE;
	return TRUE;
}

// Return the limit penalty of the parameter, based on how much it has strayed outside its allowed range
double Parameter::limitPenalty()
{
	double delta = 0.0;
	if (hasLowerLimit_ && (value_ < lowerLimit_)) delta = fabs(value_ - lowerLimit_);
	else if (hasUpperLimit_ && (value_ > upperLimit_)) delta = fabs(value_ - upperLimit_);
	// If both limits are defined, use the range within it to scale the penalty
	if (hasLowerLimit_ && hasUpperLimit_) return 100.0*(delta/(upperLimit_-lowerLimit_));
	else return 100.0*delta;
}

// Return lower limit as text
const char *Parameter::lowerLimitAsText()
{
	if (hasLowerLimit_)
	{
		static Dnchar result;
		result.sprintf("%f",lowerLimit_);
		return result;
	}
	return "<none>";
}

// Return upper limit as text
const char *Parameter::upperLimitAsText()
{
	if (hasUpperLimit_)
	{
		static Dnchar result;
		result.sprintf("%f",upperLimit_);
		return result;
	}
	return "<none>";
}

	// Return changelog value
unsigned int Parameter::changeLog()
{
	return changeLog_;
}
