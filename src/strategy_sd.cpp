/*
	*** Steepest Descent Step
	*** strategy_sd.cpp
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

// Constructor
SteepestDescentStep::SteepestDescentStep(StrategyStep::StepFunction sf) : StrategyStep(sf)
{
}

// Return cost function using parameters moved along specified gradient
double SteepestDescentStep::costAt(Alpha alpha, double* gradient, double delta, TargetSystem* targetSystem, Forcefield* ff)
{
	for (int n=0; n<alpha.nAlpha(); ++n) alpha.multiplyAlpha(n, 1.0-(gradient[n]*delta));
	return alpha.cost(targetSystem, ff);
}

// Bracket initial minimum
void SteepestDescentStep::bracket(Alpha alpha, double *gradient, double &ax, double &bx, double &cx, TargetSystem *targetSystem, Forcefield *ff)
{
	// Bracket minimum and find new point
	double fa, fb, fc;
	double goldenRatio = 1.618304, magLimit = 100.0;

	fa = costAt(alpha, gradient, ax, targetSystem, ff);
	fb = costAt(alpha, gradient, bx, targetSystem, ff);
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
	fc = costAt(alpha, gradient, cx, targetSystem, ff);
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
			fu = costAt(alpha, gradient, ux, targetSystem, ff);
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
			fu = costAt(alpha, gradient, ux, targetSystem, ff);
		}
		else if ((cx-ux)*(ux-ulim) > 0.0)
		{
			// Parabolic estimate is between cx and the imposed limit
			fu = costAt(alpha, gradient, ux, targetSystem, ff);
			if (fu < fc)
			{
				bx = cx;
				cx = ux;
				ux = ux*goldenRatio*(ux-cx);
				fb = fc;
				fc = fu;
				fu = costAt(alpha, gradient, ux, targetSystem, ff);
			}
		}
		else if ((ux-ulim)*(ulim-cx) >= 0.0)
		{
			// Limit parabolic value...
			ux = ulim;
			fu = costAt(alpha, gradient, ux, targetSystem, ff);
		}
		else
		{
			// Reject parabolic estimate, and use default ratio
			ux = cx + goldenRatio*(cx-bx);
			fu = costAt(alpha, gradient, ux, targetSystem, ff);
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
double SteepestDescentStep::goldenSearch(Alpha alpha, double *gradient, double ax, double bx, double cx, double &finalCost, TargetSystem* targetSystem, Forcefield* ff)
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
	f1 = costAt(alpha, gradient, x1, targetSystem, ff);
	f2 = costAt(alpha, gradient, x2, targetSystem, ff);
	int ncycles = 0;
	while (fabs(x3-x0) > tolerance*(fabs(x1)+fabs(x2)))
	{
		if (f2 < f1)
		{
			x0 = x1;
			x1 = x2;
			x2 = ratio*x2 + cRatio*x3;
			f1 = f2;
			f2 = costAt(alpha, gradient, x2, targetSystem, ff);
		}
		else
		{
			x3 = x2;
			x2 = x1;
			x1 = ratio*x1 + cRatio*x0;
			f2 = f1;
			f1 = costAt(alpha, gradient, x1, targetSystem, ff);
		}
		// Safety check - terminate after 100 cycles regardless
		++ncycles;
		if (ncycles == 100)
		{
			msg("Golden Search failed to converge - terminating after 100 iterations.\n");
			break;
		}
	}
	// Done
	finalCost = f1 < f2 ? f1 : f2;
	return (f1 < f2 ? x1 : x2);
}

// Print step information
void SteepestDescentStep::print(int indent)
{
	int maxCycles = parameters_[0]->asInteger();
	double convergence = parameters_[1]->asDouble();
	msg(indent, "Steepest Descent minimisation (maxCycles = %i, convergence = %8.3e)\n", maxCycles, convergence);
}

// Execute step
Alpha SteepestDescentStep::execute(Alpha inputAlpha, Forcefield* ff, TargetSystem *targetSystem)
{
	int maxCycles = parameters_[0]->asInteger();
	double convergence = parameters_[1]->asDouble();
	
	// Set current forcefield and grab optimisable parameters
	inputAlpha.poke();
	Alpha optAlpha = inputAlpha.fittable();
	int nParams = optAlpha.nAlpha();
	
	Alpha bestAlpha = optAlpha;
	double gradient[bestAlpha.nAlpha()], startCost, thisCost, lastCost;
	startCost = bestAlpha.cost(targetSystem, ff);
	lastCost = startCost;
	
	int cycle = 1;
	do
	{
		// Determine 'gradient' of parameters
		double originalValue, nudge = 0.001;
// 		msg("Cycle %i of %i (max) : current cost is %11.5e\n", cycle, maxCycles, lastCost);
		double left, right;
		int n, minimumValue = 0;
		for (n=0; n<bestAlpha.nAlpha(); ++n)
		{
			originalValue = bestAlpha.alpha(n);
			// Calculate positive displacement
			bestAlpha.setAlpha(n, originalValue + originalValue*nudge);
			right = bestAlpha.cost(targetSystem, ff) - lastCost;
			// Calculate negative displacement
			bestAlpha.setAlpha(n, originalValue - originalValue*nudge);
			left = bestAlpha.cost(targetSystem, ff) - lastCost;

			bestAlpha.setAlpha(n, originalValue);
			
			gradient[n] = (right-left)/(2.0*nudge);
// 			printf("Gradient of parameter %i = %f\n", n, gradient[n]);
		}
		
		// Normalise with respect to largest absolute gradient
		double maximumValue = 0.0, maxN;
		for (n=0; n<bestAlpha.nAlpha(); ++n) if (fabs(gradient[n]) > maximumValue)
		{
			maximumValue = fabs(gradient[n]);
			maxN = n;
		}
// 		msg("\t\tMaximal parameter gradient is %11.5e (%s)\n", maximumValue, bestAlpha.sourceParameter(maxN)->name());
		for (n=0; n<bestAlpha.nAlpha(); ++n)
		{
			gradient[n] /= maximumValue;
// 			msg("\t\tScaled gradient of parameter %i = %f\n", n, gradient[n]);
		}

		// Bracket initial minimum
		double ax = 0.0, bx = 1.0, cx;
// 		printf("Bracketing...\n");
		bracket(bestAlpha, gradient, ax, bx, cx, targetSystem, ff);
		
		// Perform Golden Search to get minimum
// 		printf("Golden search...\n");
		double step = goldenSearch(bestAlpha, gradient, ax, bx, cx, thisCost, targetSystem, ff);
// 		printf("  --> Step is %f\n", step);
		
		// Step along the gradient
		for (n=0; n<bestAlpha.nAlpha(); ++n) bestAlpha.multiplyAlpha(n, (1.0-gradient[n]*step));
		targetSystem->updateBestAlpha(&bestAlpha);
		
		// Increase cycle number and check for convergence
		++cycle;
		if (cycle > maxCycles) msg("Limit of %i cycles reached.\n", maxCycles);
		if (fabs(thisCost-lastCost) < convergence)
		{
			msg("Convergence limit achieved (step cost delta = %12.5e)\n", fabs(thisCost-lastCost));
			break;
		}
		lastCost = thisCost;
	} while (cycle <= maxCycles);
	msg("Starting cost was %11.5e, final cost is %11.5e (delta = %12.5e)\n", startCost, thisCost, lastCost - startCost);
	
	// Poke best alpha back into travelling parameters
	if (bestAlpha.cost(targetSystem, ff) < inputAlpha.cost(targetSystem, ff))
	{
		bestAlpha.poke();
		inputAlpha.peek();
	}
	else msg("Steepest descent failed to find better point than the current one.\n");
	
	return inputAlpha;
}
