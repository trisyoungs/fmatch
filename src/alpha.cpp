/*
	*** Set of Alpha
	*** alpha.cpp
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

#include "alpha.h"
#include "constants.h"
#include "ff.h"
#include "system.h"

// Constructor / Destructor
Alpha::Alpha(int size)
{
	// Private variables
	nAlpha_ = 0;
	alpha_ = NULL;
	lengthScales_ = NULL;
	sourceParameters_ = NULL;
	fit_ = NULL;
	if (size > 0) initialise(size);
	costIsValid_ = FALSE;
	cost_ = -1.0;
	
	// Public variables
	prev = NULL;
	next = NULL;
}

Alpha::~Alpha()
{
}

// Copy constructor
Alpha::Alpha(const Alpha &v)
{
	// Private variables
	nAlpha_ = 0;
	alpha_ = NULL;
	lengthScales_ = NULL;
	sourceParameters_ = NULL;
	fit_ = NULL;
	initialise(v.nAlpha_);
	copyData(v);
	costIsValid_ = FALSE;
	
	// Public variables
	prev = NULL;
	next = NULL;
}

// Initialise to specified size
void Alpha::initialise(int size)
{
	if (alpha_ != NULL) delete[] alpha_;
	if (sourceParameters_ != NULL) delete[] sourceParameters_;
	if (lengthScales_ != NULL) delete[] lengthScales_;
	if (fit_ != NULL) delete[] fit_;

	nAlpha_ = size;
	alpha_ = new double[nAlpha_];
	lengthScales_ = new double[nAlpha_];
	sourceParameters_ = new Parameter*[nAlpha_];
	fit_ = new bool[nAlpha_];
	for (int n=0; n<nAlpha_; ++n)
	{
		alpha_[n] = 0.0;
		lengthScales_[n] = 0.0;
		sourceParameters_[n] = NULL;
		fit_[n] = TRUE;
	}
	costIsValid_ = FALSE;
}

// Copy information from supplied Alpha
void Alpha::copyData(const Alpha &v)
{
	if (alpha_ == NULL) initialise(v.nAlpha_);
	else if (nAlpha_ != v.nAlpha()) initialise(v.nAlpha_);
	for (int n=0; n<nAlpha_; ++n)
	{
		alpha_[n] = v.alpha_[n];
		lengthScales_[n] = v.lengthScales_[n];
		sourceParameters_[n] = v.sourceParameters_[n];
		fit_[n] = v.fit_[n];
	}
}

// Clear all arrays and data
void Alpha::clear()
{
	if (alpha_ != NULL) delete[] alpha_;
	if (sourceParameters_ != NULL) delete[] sourceParameters_;
	if (lengthScales_ != NULL) delete[] lengthScales_;
	if (fit_ != NULL) delete[] fit_;
	nAlpha_ = 0;
	alpha_ = NULL;
	lengthScales_ = NULL;
	sourceParameters_ = NULL;
	fit_ = NULL;
	costIsValid_ = FALSE;
	cost_ = -1.0;
}

// Copy data from supplied Alpha
void Alpha::copy(Alpha *source)
{
	if (alpha_ == NULL) initialise(source->nAlpha_);
	else if (nAlpha_ != source->nAlpha()) initialise(source->nAlpha_);
	for (int n=0; n<nAlpha_; ++n)
	{
		alpha_[n] = source->alpha_[n];
		lengthScales_[n] = source->lengthScales_[n];
		sourceParameters_[n] = source->sourceParameters_[n];
		fit_[n] = source->fit_[n];
	}
}

// Set from copy of supplied parameter list
void Alpha::set(List<Parameter> *parameters)
{
	if (alpha_ == NULL) initialise(parameters->nItems());
	else if (nAlpha_ != parameters->nItems()) initialise(parameters->nItems());
	// Copy parameter values
	for (int n=0; n<nAlpha_; ++n)
	{
		alpha_[n] = (*parameters)[n]->value();
		sourceParameters_[n] = (*parameters)[n];
	}
	costIsValid_ = FALSE;
}

// Set from copy of supplied vertex (optionally with single value adjusted)
void Alpha::set(const Alpha &v, int nudgeId, double nudgeFactor)
{
	if (alpha_ == NULL) initialise(v.nAlpha());
	else if (nAlpha_ != v.nAlpha_) initialise(v.nAlpha());
	// Copy parameter values
	copyData(v);
	// Nudge alpha value if requested
	if (nudgeId != -1)
	{
		if (nudgeId >= nAlpha_) printf("BUG: Alpha index %i in Alpha is out of bounds (0 to %i).\n", nudgeId, nAlpha_-1);
		else alpha_[nudgeId] *= nudgeFactor;
	}
	costIsValid_ = FALSE;
}

// Set single value/parameter in vertex
void Alpha::set(int n, double value, Parameter* param, double lengthScale, bool fit)
{
	if (alpha_ == NULL) printf("BUG: Alpha has not been initialised before setting an individual vertex.\n");
	else if (n < 0 || n >= nAlpha_) printf("BUG: Index %i is out of range for setting a vertex vertex.\n", n);
	else
	{
		alpha_[n] = value;
		lengthScales_[n] = lengthScale;
		sourceParameters_[n] = param;
		fit_[n] = fit;
	}
}

// Multiply value in Alpha by supplied value
void Alpha::multiplyValue(int n, double d)
{
	if ((n < 0) || (n >= nAlpha_)) printf("BUG: Index to multiply (%i) in Alpha is out of range.\n", n);
	else alpha_[n] *= d;
}

// Calculate and return cost function at this vertex
double Alpha::cost(TargetSystem *system, Forcefield *ff)
{
	if (costIsValid_) return cost_;
	
	// Must calculate new cost. Begin by poking this vertex's parameters into the Forcefield...
	if (!poke()) return -1.0;
	
	// Calculate dependencies in forcefield
	if (!ff->calculateDependencies()) return -1.0;
	
	// Calculate forces and SOSE for current configuration stack
	system->stackForces();
	cost_ = system->stackSOSE();
	costIsValid_ = TRUE;
	return cost_;
}

// Reset parameter length scales
void Alpha::resetLengthScales()
{
	if (lengthScales_ == NULL) printf("BUG: No lengthScales_ array in Alpha to reset.\n");
	else for (int n=0; n<nAlpha_; ++n) lengthScales_[n] = 1.0;
}

// Calculate parameter length scales
void Alpha::calculateLengthScales(TargetSystem *system, Forcefield *ff)
{
	if (alpha_ == NULL)
	{
		printf("BUG: No lengthScales_ array in Alpha to reset.\n");
		return;
	}

	// Determine rough 'gradient' of parameters
	double originalValue, nudge = 0.01, baseRMS = cost(system, ff);
	int n, minimumValue = 0;
	for (n=0; n<nAlpha_; ++n)
	{
		originalValue = alpha_[n];
		// Calculate positive displacement
		setAlpha(n, originalValue * (1.0+nudge));
		lengthScales_[n] = fabs(cost(system, ff) - baseRMS);
		// Calculate negative displacement
		setAlpha(n, originalValue * (1.0-nudge));
		lengthScales_[n] += fabs(cost(system, ff) - baseRMS);

		setAlpha(n, originalValue);
		lengthScales_[n] = sqrt(0.5*lengthScales_[n]);
// 		printf("lenghscale %i = %f\n", n, lengthScales_[n]);
		if (lengthScales_[n] < lengthScales_[minimumValue]) minimumValue = n;
	}
	
	// Normalise with respect to lowest parameter
	for (n=0; n<nAlpha_; ++n) lengthScales_[n] = 1.0 / (lengthScales_[n] / lengthScales_[minimumValue]);
	
	// Check for stagnant length scales
	int nSmall = 0;
	for (n=0; n<nAlpha_; ++n)
	{
		if (lengthScales_[n] < 0.05)
		{
			++nSmall;
			lengthScales_[n] = 0.05;
		}
	}
	if (nSmall > (nAlpha_*0.5))
	{
		printf("[ALPHA]\tToo many parameters (%i) have length scales less than 5%% of the weakest. All scales reset to 1.0.\n", nSmall);
		resetLengthScales();
	}
	else
	{
		// Print information
		printf("\n[ALPHA]\tLength scales:");
		for (n=0; n<nAlpha_; ++n) printf(" %s %f%c", sourceParameters_[n]->name(), lengthScales_[n], n<nAlpha_-1 ? ',' : ' ');
		printf("\n");
	}
}

// Return length scale for specified parameter
double Alpha::lengthScale(int n)
{
	if (lengthScales_ == NULL) printf("BUG: No lengthScales_ array in Alpha to reset.\n");
	else return lengthScales_[n];
	return 0.0;
}

// Return number of parameters in array
int Alpha::nAlpha() const
{
	return nAlpha_;
}

// Return specified alpha
double Alpha::alpha(int n)
{
	return alpha_[n];
}

// Set value of specified alpha
void Alpha::setAlpha(int n, double value)
{
	alpha_[n] = value;
	costIsValid_ = FALSE;
}

// Multiply value of specified alpha
void Alpha::multiplyAlpha(int n, double factor)
{
	alpha_[n] *= factor;
	costIsValid_ = FALSE;
}
	
// Return pointers to original parameters
Parameter **Alpha::sourceParameters()
{
	return sourceParameters_;
}

// Returnu specified source parameter
Parameter *Alpha::sourceParameter(int n)
{
	return sourceParameters_[n];
}

// Reset all fit flags to TRUE
void Alpha::fitAll()
{
	for (int n=0; n<nAlpha_; ++n) fit_[n] = TRUE;
}

// Set individual fit flag
void Alpha::setFit(int n, bool status)
{
	fit_[n] = status;
}

// Return fit flag for parameter specified
bool Alpha::fit(int n)
{
	return fit_[n];
}

// Poke current values into associated Parameters
bool Alpha::poke()
{
	if (alpha_ == NULL) printf("BUG: No alpha values in Alpha to poke.\n");
	else for (int n=0; n<nAlpha_; ++n)
	{
		if (sourceParameters_[n] == NULL)
		{
			printf("BUG: Source parameter %i has no associated pointer.\n", n);
			return FALSE;
		}
		else if (!sourceParameters_[n]->setValue(alpha_[n])) return FALSE;
	}
	return TRUE;
}

// Retrieve current values from associated Parameters
bool Alpha::peek()
{
	if (alpha_ == NULL) printf("BUG: No alpha values in Alpha to peek.\n");
	else for (int n=0; n<nAlpha_; ++n)
	{
		if (sourceParameters_[n] == NULL)
		{
			printf("BUG: Source parameter %i has no associated pointer in peek.\n", n);
			return FALSE;
		}
		else alpha_[n] = sourceParameters_[n]->value();
	}
	return TRUE;
}

// Print vertex alpha contents
void Alpha::print()
{
	for (int n=0; n<nAlpha_; ++n)
	{
		if (fit_[n]) printf("%8.4e [O] (%s) ", alpha_[n], sourceParameters_[n] == NULL ? "???" : sourceParameters_[n]->name());
		else printf("%8.4e (%s) ", alpha_[n], sourceParameters_[n] == NULL ? "???" : sourceParameters_[n]->name());
	}
	printf("\n");
}

// Return new alpha array containing only those parameters marked as fittable
Alpha Alpha::fittable()
{
	// Count fittable parameters first...
	int n, nFit = 0, count;
	for (n=0; n<nAlpha_; ++n) if (fit_[n]) ++nFit;
	if (nFit == 0) printf("[ALPHA]\tWarning - Alpha contains no fittable parameters.\n");

	// Create new Alpha
	Alpha fittableAlpha(nFit);
	count = 0;
	for (n=0; n<nAlpha_; ++n)
	{
		if (fit_[n])
		{
			fittableAlpha.set(count, alpha_[n], sourceParameters_[n], lengthScales_[n], TRUE);
			++count;
		}
	}
	return fittableAlpha;
}

/*
// Operators / Functions
*/

void Alpha::zero()
{
	if (alpha_ == NULL) return;
	for (int n=0; n<nAlpha_; ++n) alpha_[n] = 0.0;
}


// Element access operator
double Alpha::operator[](int index)
{
	if (index < 0 || index >= nAlpha_) printf("BUG: Alpha array index %i is out of bounds.\n", index);
	else return alpha_[index];
	return 0.0;
}

void Alpha::operator=(const Alpha &v)
{
	// Check sizes
	if (alpha_ == NULL) initialise(v.nAlpha());
	else if (nAlpha_ != v.nAlpha_) initialise(v.nAlpha());
// 	{
// 		printf("BUG: Supplied Alpha is of a different size (%i vs %i) to that previously stored in this Alpha.\n", v.nAlpha_, nAlpha_);
// 		return;
// 	}
	// Copy data
	copyData(v);
	cost_ = v.cost_;
	costIsValid_ = v.costIsValid_;
}

void Alpha::operator+=(const Alpha &v)
{
	if (alpha_ == NULL)
	{
		initialise(v.nAlpha());
		// Can just copy data of other vertex...
		copyData(v);
	}
	else if (nAlpha_ != v.nAlpha_)
	{
		printf("BUG: Operator += Alpha is of a different size (%i vs %i) to that previously stored in this Alpha.\n", v.nAlpha_, nAlpha_);
		return;
	}
	else for (int n=0; n<nAlpha_; ++n)
	{
		alpha_[n] += v.alpha_[n];
// 		if (
	}
	costIsValid_ = FALSE;
}

void Alpha::operator-=(const Alpha &v)
{
	if (alpha_ == NULL)
	{
		printf("BUG: Operator -= Alpha is not initialised.\n");
		return;
	}
	else if (nAlpha_ != v.nAlpha_)
	{
		printf("BUG: Operator -= Alpha is of a different size (%i vs %i) to that previously stored in this Alpha.\n", v.nAlpha_, nAlpha_);
		return;
	}
	for (int n=0; n<nAlpha_; ++n) alpha_[n] -= v.alpha_[n];
	costIsValid_ = FALSE;
}

void Alpha::operator/=(const double d)
{
	if (alpha_ == NULL)
	{
		printf("BUG: Operator /= Alpha is not initialised.\n");
		return;
	}
	for (int n=0; n<nAlpha_; ++n) alpha_[n] /= d;
	costIsValid_ = FALSE;
}

Alpha Alpha::operator*(double d) const
{
	Alpha vertex(*this);
	if (alpha_ == NULL) printf("BUG: Attempted to operate (*) on an Alpha that hasn't been initialised yet.\n"); 
	else for (int n=0; n<nAlpha_; ++n) vertex.alpha_[n] *= d;
	return vertex;
}

Alpha Alpha::operator+(const Alpha &v) const
{
	Alpha vertex(nAlpha_);
	if (alpha_ == NULL) printf("BUG: Attempted to operate (+) on an Alpha that hasn't been initialised yet.\n"); 
	else
	{
		vertex.copyData(*this);
		for (int n=0; n<nAlpha_; ++n) vertex.alpha_[n] += v.alpha_[n];
	}
	return vertex;
}

Alpha Alpha::operator-(const Alpha &v) const
{
	Alpha vertex(nAlpha_);
	if (alpha_ == NULL) printf("BUG: Attempted to operate (-) on an Alpha that hasn't been initialised yet.\n"); 
	else
	{
		vertex.copyData(*this);
		for (int n=0; n<nAlpha_; ++n) vertex.alpha_[n] -= v.alpha_[n];
	}
	return vertex;
}
