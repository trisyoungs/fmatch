/*
	*** Alpha Store
	*** alphastore.cpp
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

#include "alphastore.h"
#include "constants.h"
#include "ff.h"
#include "system.h"

// Constructor
AlphaStore::AlphaStore()
{
	// Private variables
	nAlpha_ = -1;
	minima_ = NULL;
	maxima_ = NULL;
	mean_ = NULL;
	sd_ = NULL;
}

// Destructor
AlphaStore::~AlphaStore()
{
	clear();
}

// Clear all data from store
void AlphaStore::clear()
{
	alpha_.clear();
	if (minima_ != NULL) delete[] minima_;
	if (maxima_ != NULL) delete[] maxima_;
	if (mean_ != NULL) delete[] mean_;
	if (sd_ != NULL) delete[] sd_;
	minima_ = NULL;
	maxima_ = NULL;
	mean_ = NULL;
	sd_ = NULL;
}

// Add supplied alpha into store
void AlphaStore::add(Alpha *source, TargetSystem *targetSystem, Forcefield *ff)
{
	// Insert in cost order....
	Alpha *newAlpha, *listPos;

	// First addition to list? Then also create arrays...
	if (alpha_.nItems() == 0)
	{
		nAlpha_ = source->nAlpha();
		if (nAlpha_ < 1) return;
		minima_ = new double[nAlpha_];
		maxima_ = new double[nAlpha_];
		mean_ = new double[nAlpha_];
		sd_ = new double[nAlpha_];
		for (int n=0; n<nAlpha_; ++n)
		{
			minima_[n] = source->alpha(n);
			maxima_[n] = source->alpha(n);
			mean_[n] = source->alpha(n);
			sd_[n] = 0.0;
		}
		newAlpha = alpha_.add();
	}
	else
	{
		for (listPos = alpha_.first(); listPos != NULL; listPos = listPos->next)
		{
			if (source->cost(targetSystem, ff) < listPos->cost(targetSystem, ff))
			{
				listPos = listPos->prev;
				break;
			}
		}
		newAlpha = alpha_.insert(listPos);

		// Accumulate values
		double newMean, a;
		int nAccumulated = alpha_.nItems() - 1;
		for (int n=0; n<nAlpha_; ++n)
		{
			a = source->alpha(n);
			if (a < minima_[n]) minima_[n] = a;
			if (a > maxima_[n]) maxima_[n] = a;
			// Update mean and SD
			newMean = (mean_[n]*nAccumulated + a) / (nAccumulated+1);
			sd_[n] = (nAccumulated*sd_[n] + (a-mean_[n])*(a-newMean)) / (nAccumulated+1);
			mean_[n] = newMean;
		}

	}
	newAlpha->copy(source);
}

// Return number of sets of alpha currently in store
int AlphaStore::nData()
{
	return alpha_.nItems();
}

// Predict new set of alpha from stored vertex information
Alpha AlphaStore::predict(TargetSystem *targetSystem, Forcefield *ff)
{
	// Simple TEST - take average parameters from list
// 	printf("Taking stats from %i vertices\n", alpha_.nItems());
	double total, rCost;
// 	Alpha newAlpha = (*alpha_.first());
	Alpha newAlpha;
	newAlpha.copy(alpha_.first());
	for (int n=0; n<newAlpha.nAlpha(); ++n) newAlpha.setAlpha(n, 0.0);
	for (Alpha *alpha = alpha_.first(); alpha != NULL; alpha = alpha->next)
	{
		rCost = 1.0 / alpha->cost(targetSystem, ff);
		newAlpha += (*alpha);//* rCost;
		total += rCost;
	}
	newAlpha /= alpha_.nItems(); //rCost;
// 	stats.copyValues(&newAlpha);
	return newAlpha;
}

// Print statistics for alpha values
void AlphaStore::print()
{
	if (alpha_.nItems() == 0)
	{
		printf("No statistics to print, since no alpha sets were added.\n");
		return;
	}
	
	printf("\nStatistics gathered over %i data for each variable:\n", alpha_.nItems());
	printf("         Minimum        Maximum          S.D.         Average    Parameter\n");
	Parameter **params = alpha_[0]->sourceParameters();
	for (int n=0; n<alpha_[0]->nAlpha(); ++n)
	{
		if (params == NULL) printf("  %i  %13.5e  %13.5e  %13.5e  %13.5e  ???\n", n+1, minima_[n], maxima_[n], sd_[n], mean_[n]);
		else printf("  %i  %13.5e  %13.5e  %13.5e  %13.5e  %s\n", n+1, minima_[n], maxima_[n], sd_[n], mean_[n], params[n]->name());
	}
}

// Copy average values into supplied Alpha
void AlphaStore::copyAverages(Alpha *alpha)
{
	if (alpha->nAlpha() != nAlpha_)
	{
		printf("BUG: Can't copy average values to an uninitialised Alpha.\n");
		return;
	}
	for (int i=0; i<nAlpha_; ++i) alpha->setAlpha(i, mean_[i]);
}
