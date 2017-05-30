/*
	*** Optimisation Strategy Definition
	*** strategy.h
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

#ifndef FM3_STRATEGY_H
#define FM3_STRATEGY_H

#include "alpha.h"

// Forward Declarations
class Forcefield;
class Configuration;
class TargetSystem;

// Step Parameter
class StepParameter
{
	public:
	// Constructor / Destructor
	StepParameter();
	~StepParameter();
	// List pointers
	StepParameter *prev, *next;
	
	private:
	// Integer value
	int valueI_;
	// Double value
	double valueD_;
	// Flag specifying whether an integer or a double is stored
	bool isDouble_;
	
	public:
	// Set data (integer)
	void set(int i);
	// Set data (double)
	void set(double d);
	// Return as integer
	int asInteger();
	// Return as double
	double asDouble();
};

// Optimisation Strategy Step
class StrategyStep
{
	public:
	// Step Functions
	enum StepFunction { Clear, EndLoop, GeneticMinimise, GridMinimise, Loop, Predict, Print, Randomise, Rattle, SD, Select, SimplexAnneal, SimplexAnneal2, SimplexMinimise, Store, SuperSimplex, HiddenFeature, nStepFunctions };
	static StepFunction stepFunction(const char *s);
	static const char *stepFunction(StepFunction sf);
	// Constructor / Destructor
	StrategyStep(StrategyStep::StepFunction sf);
	~StrategyStep();
	// List pointers
	StrategyStep *prev, *next;

	protected:
	// Step function
	StepFunction function_;
	// List of parameters to step
	List<StepParameter> parameters_;

	protected:
	// Print step-tagged output
	void msg(const char *fmt, ...);
	// Print indent = 0ed step-tagged output
	void msg(int indent, const char *fmt, ...);

	public:
	// Add parameter to step (integer)
	void addParameter(int i);
	// Add parameter to step (double)
	void addParameter(double d);
	// Print step information
	virtual void print(int indent = 0) = 0;
	// Execute step
	virtual Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem) = 0;
};

/*
// Classes Derived from StrategyStep
*/

// Genetic Algorithm
class GeneticStep : public StrategyStep
{
	public:
	// Constructor
	GeneticStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Loop
class LoopStep : public StrategyStep
{
	public:
	// Constructor
	LoopStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
	
	private:
	// List of steps to loop over
	List<StrategyStep> strategy_;
	
	public:
	// Return strategy list
	List<StrategyStep> *strategy();
};

// Print Current Parameters
class PrintStep : public StrategyStep
{
	public:
	// Constructor
	PrintStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Print Summary of Current Parameters
class SummaryStep : public StrategyStep
{
	public:
	// Constructor
	SummaryStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Select Parameters
class SelectStep : public StrategyStep
{
	public:
	// Constructor
	SelectStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Parameter Randomisation
class RandomisationStep : public StrategyStep
{
	public:
	// Constructor
	RandomisationStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Simplex Anneal
class SimplexAnnealStep : public StrategyStep
{
	public:
	// Constructor
	SimplexAnnealStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Simplex Anneal 2
class SimplexAnneal2Step : public StrategyStep
{
	public:
	// Constructor
	SimplexAnneal2Step();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Simplex Minimisation
class SimplexMinimiseStep : public StrategyStep
{
	public:
	// Constructor
	SimplexMinimiseStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Super Simplex
class SuperSimplexStep : public StrategyStep
{
	public:
	// Constructor
	SuperSimplexStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Hidden Feature
class HiddenFeatureStep : public StrategyStep
{
	public:
	// Constructor
	HiddenFeatureStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Grid Search Minimise
class GridMinimiseStep : public StrategyStep
{
	public:
	// Constructor
	GridMinimiseStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Steepest Descent
class SteepestDescentStep : public StrategyStep
{
	private:
	// Return cost function using parameters moved along specified gradient
	double costAt(Alpha alpha, double *gradient, double delta, TargetSystem *targetSystem, Forcefield *ff);
	// Bracket initial minimum
	void bracket(Alpha alpha, double *gradient, double &ax, double &bx, double &cx, TargetSystem *targetSystem, Forcefield *ff);
	// Perform Golden Search to find minimum
	double goldenSearch(Alpha alpha, double *gradient, double ax, double bx, double cx, double &finalCost, TargetSystem *targetSystem, Forcefield *ff);

	public:
	// Constructor
	SteepestDescentStep(StrategyStep::StepFunction sf = StrategyStep::SD);
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Predict Step
class PredictStep : public StrategyStep
{
	private:
	// Return cost function using parameters moved along specified gradient
	double costAt(Alpha alpha, double* a, double* b, double delta, TargetSystem* targetSystem, Forcefield* ff);
	// Bracket initial minimum
	void bracket(Alpha alpha, double *a, double *b, double &ax, double &bx, double &cx, TargetSystem *targetSystem, Forcefield *ff);
	// Perform Golden Search to find minimum
	double goldenSearch(Alpha alpha, double *a, double *b, double ax, double bx, double cx, double &finalCost, TargetSystem *targetSystem, Forcefield *ff);

	public:
	// Constructor
	PredictStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Clear Step
class ClearStep : public StrategyStep
{
	public:
	// Constructor
	ClearStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Store Step
class StoreStep : public StrategyStep
{
	public:
	// Constructor
	StoreStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

// Rattle Step
class RattleStep : public StrategyStep
{
	public:
	// Constructor
	RattleStep();
	// Print step information
	void print(int indent = 0);
	// Execute step
	Alpha execute(Alpha inputAlpha, Forcefield *ff, TargetSystem *targetSystem);
};

#endif
