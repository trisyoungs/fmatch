/*
	*** Genetic Algorithm
	*** genetic.cpp
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

#include "genetic.h"
#include "parameter.h"
#include "constants.h"
#include "ff.h"
#include "system.h"
#include "random.h"
#include "alphastore.h"

// Constructor
Genetic::Genetic(TargetSystem* system, Forcefield* ff)
{
	// Private variables
	nGenes_ = 0;
	forcefield_ = ff;
	targetSystem_ = system;
	betterPointFound_ = FALSE;
}

Genetic::~Genetic()
{
}

// Return (calculating if necessary) cost of specified gene
double Genetic::cost(int n)
{
	// Compute cold cost of gene first, to compare against bestAlpha, then add random T fluctuation
	Alpha *gene = genePool_[n];
	double coldCost = gene->cost(targetSystem_, forcefield_);
	// Check best point in *this* step
	if (coldCost < bestAlpha_.cost(targetSystem_, forcefield_))
	{
		bestAlpha_ = (*gene);
		betterPointFound_ = TRUE;
	}
	// Check best overall point stored
	targetSystem_->updateBestAlpha(gene);
	return coldCost;
}

// Return (calculating if necessary) cost of supplied gene
double Genetic::cost(Alpha *gene)
{
	// Compute cold cost of vertex first, to compare against bestAlpha, then add random T fluctuation
	double coldCost = gene->cost(targetSystem_, forcefield_);
	// Check best point in *this* step
	if (coldCost < bestAlpha_.cost(targetSystem_, forcefield_))
	{
		bestAlpha_ = (*gene);
		betterPointFound_ = TRUE;
	}
	// Check best overall point stored
	targetSystem_->updateBestAlpha(gene);
	return coldCost;
}

// Sort genes according to cost, with lowest first in array
void Genetic::sort()
{
	Alpha *gene, *currentGene = genePool_.first(), *highest = NULL, *nextGene;

	// Selection sort
	int count = 0;
	do
	{
		// Find *highest* value in remaining list
		highest = NULL;
		for (gene = currentGene; gene != NULL; gene = gene->next)
		{
			if (highest == NULL) highest = gene;
			else if (gene->cost(targetSystem_,forcefield_) > highest->cost(targetSystem_,forcefield_)) highest = gene;
		}
		// Move highest to very start of list (which means that the lowest will eventually end up at the start also).
		// Need to find new currentGene before we move the item
		nextGene = (highest == currentGene ? currentGene->next : currentGene);
		// Now move item to start of list
		genePool_.moveToStart(highest);
		// Set new current gene pointer (which might still be the same gene as before...)
		currentGene = nextGene;
	} while (currentGene);
// 	for (gene = genePool_.first(); gene != NULL; gene = gene->next) printf("Sort: Gene %p cost = %f\n", gene, gene->cost(targetSystem_, forcefield_));
}

// Generate initial population
void Genetic::initialise(Alpha &initVertex, int populationSize, double baselengthscale)
{
	nGenes_ = populationSize;
	bestAlpha_ = initVertex;
	nAlpha_ =  initVertex.nAlpha();
	betterPointFound_ = FALSE;
	lengthScale_ = baselengthscale;

	Alpha *gene;
	double lscale, r;
	int n, m;
	for (n=0; n<nGenes_; ++n)
	{
		gene = genePool_.add();
		
		(*gene) = initVertex;
		
		// Randomise all genes except first one
		if (n > 0)
		{
			r = (2.0*Random::number()) - 1.0;
			for (m=0; m<nAlpha_; ++m)
			{
				lscale = lengthScale_ * gene->lengthScale(m);
				gene->multiplyAlpha(m, 1.0+lscale*r);
			}
		}
	}
}

// Return MSD of cost values in current Genetic
double Genetic::costMSD()
{
	// Get average cost of all vertices
	double average = 0.0;
	for (int n=0; n<nGenes_; ++n) average += cost(n);
	average /= nGenes_;
	// Now calculate SD of individual costs with mean value
	double serror = 0.0;
	for (int n=0; n<nGenes_; ++n) serror += (cost(n) - average)*(cost(n) - average);
	return sqrt(serror/nGenes_);
}

// Set specified gene
void Genetic::setGene(int n, Alpha geneInfo)
{
	*(genePool_[n]) = geneInfo;
}

// Return specified gene
Alpha Genetic::gene(int n)
{
	return *(genePool_[n]);
}

// Minimise simplex for supplied configuration and forcefield
Alpha Genetic::evolve(double selectionThreshold, double mutateThreshold)
{
	int n, i, cycle;
	Alpha *gene, *geneA, *geneB;

	// Check nGenes
	if (nGenes_ == 0)
	{
		printf("Error: Genetic Algorithm appears not to have been initialised properly.\n");
		targetSystem_->flagAbort();
		return bestAlpha_;
	}

	// If selection thresold leaves too few genes, we cannot continue
	if ((nGenes_*selectionThreshold) < 3)
	{
		printf("Error: Genetic Algorithm selection threshold is too low, leaving fewer than 3 genes per cycle.\n");
		targetSystem_->flagAbort();
		return bestAlpha_;
	}
	
	// Sort genes by cost
	sort();

	// Perform truncation selection - top selectionThreshold genes will mate to produce  N-2 new genes

	// Delete rejected genes
	for (i=0; i<nGenes_*(1.0-selectionThreshold); ++i) genePool_.removeLast();
	// Determine min, max, and average cost of remaining parents
	double mincost = 1.0e9, maxcost = -1.0e9, avg = 0.0;
	AlphaStore alphaStats;
	for (Alpha *gene = genePool_.first(); gene != NULL; gene = gene->next)
	{
		if (cost(gene) < mincost) mincost = cost(gene);
		if (cost(gene) > maxcost) maxcost = cost(gene);
		avg += cost(gene);
		alphaStats.add(gene, targetSystem_, forcefield_);
	}
	avg /= genePool_.nItems();
	printf("[SOUP]\t\tAfter truncation, %i remaining genes have cost range of %12.6e to %12.6e, with average of %12.6e\n", genePool_.nItems(), mincost, maxcost, avg);
// 	printf("\n[SOUP]\t\tParameter spread over remaining %i genes is :\n", genePool_.nItems());
// 	alphaStats.print(genePool_[0]);
	
	// Move all remaining genes to temporary list
	List<Alpha> oldParents;
	oldParents = genePool_;
	// Clear main pool, leaving best gene (first one) behind, and add on Alpha created from average values
	for (i=1; i<genePool_.nItems(); ++i) genePool_.removeLast();
	gene = genePool_.add();
	(*gene) = (*genePool_[0]);
	alphaStats.copyAverages(gene);
	
	// Generate new genes until pool is full again
	do
	{
		// Randomly select two old parent genes
		geneA = oldParents[int(Random::number()*oldParents.nItems())];
		do
		{
			geneB = oldParents[int(Random::number()*oldParents.nItems())];
		} while (geneA == geneB);

		// Add a new gene, and copy gene[0] to set up structure properly
		gene = genePool_.add();
		(*gene) = (*genePool_[0]);

		// Perform Discrete Recombination
		for (i=0; i<gene->nAlpha(); ++i)
		{
			if (i < 0.5*gene->nAlpha()) gene->setAlpha(i, geneA->alpha(i));
			else gene->setAlpha(i, geneB->alpha(i));
			// Mutate?
			if (Random::number() < mutateThreshold) gene->multiplyAlpha(i, (1.0-lengthScale_) + 2.0 * lengthScale_ * Random::number());
		}
	} while (genePool_.nItems() < nGenes_);

	// Random mutations?
	// TODO
	
	// Sort genes once more (also gets costs of new genes in the process)
	sort();

	return bestAlpha_;
}

// Return whether top two genes are converged within some threshold
bool Genetic::converged(double threshold, double selectionThreshold)
{
	int nKeep = nGenes_ * selectionThreshold;
	for (int n=0; n<nKeep; ++n) if (fabs(cost(genePool_[n])/cost(genePool_[0])-1.0) > threshold) return FALSE;
	return TRUE;
}

// Refresh pool with new genes, using the current best as a basis
void Genetic::refresh()
{
	Alpha *gene;
	double lscale, r;
	int n, m;
	for (n=1; n<nGenes_; ++n)
	{
		gene = genePool_[n];
		
		(*gene) = (*genePool_[0]);
		
		r = (2.0*Random::number()) - 1.0;
		for (m=0; m<nAlpha_; ++m)
		{
			lscale = lengthScale_ * gene->lengthScale(m);
			gene->multiplyAlpha(m, 1.0+lscale*r);
		}
	}
}

// Return whether a better point was found by the algorithm
bool Genetic::betterPointFound()
{
	return betterPointFound_;
}
