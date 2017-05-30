/*
	*** Predict Step
	*** strategy_predict.cpp
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
#include "ff.h"
#include "simplex.h"

/*
// Clear Step
*/

// Constructor
ClearStep::ClearStep() : StrategyStep(StrategyStep::Clear)
{
}

// Print step information
void ClearStep::print(int indent)
{
	int id = parameters_[0]->asInteger();
	msg(indent, "Clear stored alpha list (store id %i).\n", id);
}

// Execute step
Alpha ClearStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int id = parameters_[0]->asInteger() - 1;
	targetSystem->clearAlphaList(id);
	return inputAlpha;
}

/*
// Store Step
*/

// Constructor
StoreStep::StoreStep() : StrategyStep(StrategyStep::Store)
{
}

// Print step information
void StoreStep::print(int indent)
{
	int id = parameters_[0]->asInteger();
	msg(indent, "Store current alpha in list (store id %i).\n", id);
}

// Execute step
Alpha StoreStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int id = parameters_[0]->asInteger() - 1;
	targetSystem->storeAlpha(id, &inputAlpha);
	return inputAlpha;
}

/*
// Predict Step
*/

// Constructor
PredictStep::PredictStep() : StrategyStep(StrategyStep::Predict)
{
}

// Return cost function using parameters moved along specified gradient
double PredictStep::costAt(Alpha alpha, double* a, double *b, double delta, TargetSystem* targetSystem, Forcefield* ff)
{
	// for (int n=0; n<alpha.nAlpha(); ++n) alpha.setAlpha(n, a[n] + b[n]*delta);
	for (int n=0; n<alpha.nAlpha(); ++n) alpha.setAlpha(n, alpha.alpha(n) + b[n]*delta);
	return alpha.cost(targetSystem, ff);
}

// Bracket initial minimum
void PredictStep::bracket(Alpha alpha, double *a, double *b, double &ax, double &bx, double &cx, TargetSystem *targetSystem, Forcefield *ff)
{
	// Bracket minimum and find new point
	double fa, fb, fc;
	double goldenRatio = 1.618304, magLimit = 100.0;

	fa = costAt(alpha, a, b, ax, targetSystem, ff);
	fb = costAt(alpha, a, b, bx, targetSystem, ff);
	// fa should be lower than fb, but check just in case
	if (fb > fa)
	{
		fc = fa;
		fa = fb;
		fb = fc;
		cx = ax;
		ax = bx;
		bx = cx;
	}
	// Set initial value of cx
	cx = bx + goldenRatio*(bx-ax);
	fc = costAt(alpha, a, b, cx, targetSystem, ff);
	// Iterate until we are converged
	double r, q, ux, ulim, x, fu;
	while (fb > fc)
	{
		r = (bx-ax)*(fb-fc);
		q = (bx-cx)*(fb-fa);
		x = max(fabs(q-r), 1.0e-10);
		ux = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*(x * (q-r < 0 ? -1.0 : 1.0)));
		ulim = bx + magLimit*(cx-bx);
		// Test new point
		if ((bx-ux)*(ux-cx) > 0.0)
		{
			// Parabolic minimum 'ux' is between points bx and cx...
			fu = costAt(alpha, a, b, ux, targetSystem, ff);
			if (fu < fc)
			{
				// Minimum is between bc and cx, so replace a with b and set b to u
				ax = bx;
				bx = ux;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb)
			{
				// Got a minimum between ax and ux
				cx = ux;
				fc = fu;
				return;
			}
			// Parabolic estimate no good....
			ux = cx+goldenRatio*(cx-bx);
			fu = costAt(alpha, a, b, ux, targetSystem, ff);
		}
		else if ((cx-ux)*(ux-ulim) > 0.0)
		{
			// Parabolic estimate is between cx and the imposed limit
			fu = costAt(alpha, a, b, ux, targetSystem, ff);
			if (fu < fc)
			{
				bx = cx;
				cx = ux;
				ux = ux*goldenRatio*(ux-cx);
				fb = fc;
				fc = fu;
				fu = costAt(alpha, a, b, ux, targetSystem, ff);
			}
		}
		else if ((ux-ulim)*(ulim-cx) >= 0.0)
		{
			// Limit parabolic value...
			ux = ulim;
			fu = costAt(alpha, a, b, ux, targetSystem, ff);
		}
		else
		{
			// Reject parabolic estimate, and use default ratio
			ux = cx + goldenRatio*(cx-bx);
			fu = costAt(alpha, a, b, ux, targetSystem, ff);
		}
		// Update points
		ax = bx;
		bx = cx;
		cx = ux;
		fa = fb;
		fb = fc;
		fc = fu;
	}
}

// Perform Golden Search to find minimum
double PredictStep::goldenSearch(Alpha alpha, double *a, double *b, double ax, double bx, double cx, double &finalCost, TargetSystem* targetSystem, Forcefield* ff)
{
	double ratio = 0.61803399, cRatio = 1.0-ratio, x0, x1, x2, x3, f1, f2, tolerance = 1.0e-5;
	// Sort out initial values, making x0-x1 the smaller segment
	x0 = ax;
	x3 = cx;
	if (fabs(cx-bx) > fabs(bx-ax))
	{
		x1 = bx; 
		x2 = bx+cRatio*(cx-bx);
	}
	else
	{
		x2 = bx; 
		x1 = bx+cRatio*(bx-ax);
	}
	// Get initial function values
	f1 = costAt(alpha, a, b, x1, targetSystem, ff);
	f2 = costAt(alpha, a, b, x2, targetSystem, ff);
	while (fabs(x3-x0) > tolerance*(fabs(x1)+fabs(x2)))
	{
		if (f2 < f1)
		{
			x0 = x1;
			x1 = x2;
			x2 = ratio*x2 + cRatio*x3;
			f1 = f2;
			f2 = costAt(alpha, a, b, x2, targetSystem, ff);
		}
		else
		{
			x3 = x2;
			x2 = x1;
			x1 = ratio*x1 + cRatio*x0;
			f2 = f1;
			f1 = costAt(alpha, a, b, x1, targetSystem, ff);
		}
	}
	// Done
	finalCost = f1 < f2 ? f1 : f2;
	return (f1 < f2 ? x1 : x2);
}

// Print step information
void PredictStep::print(int indent)
{
	int id = parameters_[0]->asInteger();
	int minSize = parameters_[1]->asInteger();
	int sampleSize = parameters_[2]->asInteger();
	msg(indent, "Predict new forcefield parameters from previous best (store id %i) using last %i sets from minimum of %i.\n", id, sampleSize, minSize);
}

// Execute step
Alpha PredictStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int id = parameters_[0]->asInteger() - 1;
	int minSize = parameters_[1]->asInteger();
	int sampleSize = parameters_[2]->asInteger();
	// Check current population in best alpha list
	if (targetSystem->alphaListSize(id) < minSize)
	{
		msg("Not enough alpha sets have been acquired to make a prediction at this point.\n");
		return inputAlpha;
	}

	// Get some data, including the last set of alpha which will be the best in the list
	Alpha bestAlpha, *alpha;
	double bestCost;
	bestAlpha.copy(targetSystem->alphaList(id));
	bestCost = bestAlpha.cost(targetSystem, ff);
	int nAlpha = bestAlpha.nAlpha();
	double a[nAlpha], b[nAlpha], chi2[nAlpha];

	// Determine 'gradient' of parameters by linear regression of best alpha list
	double sumx, sosx, sumy, avgx, dx, x, y;
	int n, i;
	for (n=0; n<nAlpha; ++n)
	{
		// Get sums of x and y data 
		sumx = 0.0;
		sumy = 0.0;
		alpha = targetSystem->alphaList(id);
		for (i = 0; i < sampleSize; ++i)
		{
			x = alpha->cost(targetSystem, ff);
			y = alpha->alpha(n);
			sumx += x;
			sumy += y;
			alpha = alpha->next;
		}

		// Determine average x, SOS diffs of x with average, and xdelta*y
		avgx = sumx / sampleSize;
		sosx = 0.0;
		b[n] = 0.0;
		alpha = targetSystem->alphaList(id);
		for (i = 0; i < sampleSize; ++i)
		{
			x = alpha->cost(targetSystem, ff);
			y = alpha->alpha(n);
			dx = x - avgx;
			sosx += dx*dx;
			b[n] += dx * y;
			alpha = alpha->next;
		}
		
		// Now solve for a and b in 'y = a + bx', and calculate chi2
		b[n] /= sosx;
		a[n] = (sumy-sumx*b[n])/sampleSize;
		chi2[n] = 0.0;
		alpha = targetSystem->alphaList(id);
		for (i = 0; i < sampleSize; ++i)
		{
			x = alpha->cost(targetSystem, ff);
			y = alpha->alpha(n);
			chi2[n] += (y - a[n] - b[n]*x)*(y - a[n] - b[n]*x);
			alpha = alpha->next;
		}

		msg("\tRegression for parameter %3i : y = %8.4e + %8.4ex,  chi2=%f (%s)\n", n+1, a[n], b[n], chi2[n], bestAlpha.sourceParameter(n)->name());
	}

	// Try to find better set of parameters
	Alpha newAlpha, newBest;

	// Extrapolate to zero
	newAlpha = bestAlpha;
	for (n=0; n<nAlpha; ++n) newAlpha.setAlpha(n, a[n]);
	msg("Extrapolation to zero gives cost of %11.5e\n", newAlpha.cost(targetSystem, ff));

	// Extrapolate to minimum - if we fail, try again but excluding the worst parameter each time
// 	for (i = 0; i<200; ++i)
// 	{
// 		for (n=0; n<nAlpha; ++n) newAlpha.setAlpha(n, a[n]+b[n]*(bestAlpha.cost(targetSystem,ff)+(bestAlpha.cost(targetSystem,ff)/100.0)*i));
// 		printf("trial step %i  line = %12.5e  actual cost = %11.5e\n", i, bestAlpha.cost(targetSystem,ff)+(bestAlpha.cost(targetSystem,ff)/100.0)*i, newAlpha.cost(targetSystem, ff));
// 	}
	
	double ax, bx, cx;
	for (i=0; i<nAlpha; ++i)
	{
		// Bracket initial minimum
		//ax = 0.0;
		ax = -bestAlpha.cost(targetSystem, ff);
		bx = bestAlpha.cost(targetSystem, ff);
		cx = 0.0;
		bracket(bestAlpha, a, b, ax, bx, cx, targetSystem, ff);
		msg(" -- brackets for minimum (extrapolating %i parameters) are %f and %f\n", nAlpha-i, ax, cx);
		
		// Perform Golden Search to get minimum
		double finalCost;
		double step = goldenSearch(bestAlpha, a, b, ax, bx, cx, finalCost, targetSystem, ff);
		//for (n=0; n<nAlpha; ++n) newAlpha.setAlpha(n, a[n] + b[n]*step);
		for (n=0; n<nAlpha; ++n) newAlpha.setAlpha(n, bestAlpha.alpha(n) + b[n]*step);
		msg(" -- minimum is at %12.5e giving cost of %11.5e\n", step, newAlpha.cost(targetSystem, ff));

		if (targetSystem->updateBestAlpha(&newAlpha))
		{
			//TEST Simplex this point...
			// Initialise Simplex and determine temperature
			Simplex simplex(targetSystem, ff);
			simplex.initialise(newAlpha, 0.05);
			// Minimise
			Alpha smplxAlpha = simplex.minimise(100, 0.01, 0.0);
			if (smplxAlpha.cost(targetSystem, ff) < newAlpha.cost(targetSystem, ff))
			{
				msg(" -- Simplex around predicted point found better minimum.\n");
				return smplxAlpha;
			}
			else msg(" -- accepting extrapolation of %i parameters to %11.5e.\n", nAlpha-i, step);
			return newAlpha;
		}

		// Couldn't find a better point, so knock out the worst parameter
		int worst = 0;
		for (n=0; n<nAlpha; ++n) if (chi2[n] > chi2[worst]) worst = n;
		b[worst] = 0.0;
		// TEST Kill stored data
		printf("Prediction failed badly - clearing current parameter store.\n");
		targetSystem->clearAlphaList(id);
		break;
	}	

	// Step along the gradient
// 	for (n=0; n<bestAlpha.nAlpha(); ++n) bestAlpha.multiplyAlpha(n, (1.0-b[n]*step));
// 	targetSystem->updateBestAlpha(&bestAlpha);
	
// 	if (foundBetter) return newBest;
	return inputAlpha;
}
