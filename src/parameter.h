/*
	*** Interaction Parameter
	*** parameter.h
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

#ifndef FM3_PARAMETER_H
#define FM3_PARAMETER_H

#include "dnchar.h"

// Forward declarations
class Interaction;

// Interaction Parameter
class Parameter
{
	public:
	// Constructor / Destructor
	Parameter();
	~Parameter();
	// List pointers
	Parameter *prev, *next;
	// Type of parameter
	enum ParameterType { FixedParameter, FreeParameter, DependentParameter };
	// Role of parameter
	enum ParameterStrength { Strong, Moderate, Weak, Any, nParameterStrengths };
	static ParameterStrength parameterStrength(const char *s);
	static const char *parameterStrength(Parameter::ParameterStrength ps);


	/*
	// Basic Data
	*/
	private:
	// Name of this parameter
	Dnchar name_;
	// Type of the parameter
	ParameterType type_;
	// Perceived strength of the parameter
	ParameterStrength strength_;
	// Current value of the parameter
	double value_;
	// Parent interaction
	Interaction *parent_;
	// Value limits, beyond which penalty occurs
	double lowerLimit_, upperLimit_;
	// Flags specifying whether limit(s) have been set
	bool hasLowerLimit_, hasUpperLimit_;
	// Integer log of parameter changes
	unsigned int changeLog_;

	public:
	// Set name and value of this parameter
	void set(const char *s, double value);
	// Return name of this parameter
	const char *name();
	// Set parameter type
	void setType(ParameterType type, ParameterStrength strength);
	// Return parameter type
	ParameterType type();
	// Return parameter strength
	ParameterStrength strength();
	// Set current value of the parameter
	bool setValue(double d);
	// Return current value of the parameter
	double value();
	// Set parent interaction
	void setParent(Interaction *parent);
	// Return parent interaction
	Interaction *parent();
	// Return parent interaction's functional form
	const char *parentForm();
	// Return parent interaction's involved atoms
	const char *parentAtoms();
	// Set lower limit for parameter
	void setLowerLimit(double minValue);
	// Set upper limit for parameter
	void setUpperLimit(double maxValue);
	// Return whether parameter has any associated limits
	bool hasLimits();
	// Return whether parameter has a lower limit
	bool hasLowerLimit();
	// Return whether parameter has an upper limit
	bool hasUpperLimit();
	// Return whether current parameter value is within accepted limits
	bool withinLimits();
	// Return the limit penalty of the parameter, based on how much it has strayed outside its allowed range
	double limitPenalty();
	// Return lower limit as text
	const char *lowerLimitAsText();
	// Return upper limit as text
	const char *upperLimitAsText();
	// Return changelog value
	unsigned int changeLog();
};

#endif
