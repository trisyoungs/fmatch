/*
	*** Genetic Algorithm 
	*** genetic.h
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

#ifndef FM3_GENETIC_H
#define FM3_GENETIC_H

#include "alpha.h"
#include "list.h"
#include "reflist.h"

// Forward Declarations
class Forcefield;
class TargetSystem;

// Genetic Algorithm
class Genetic
{
	public:
	// Constructor / Destructor
	Genetic(TargetSystem *system, Forcefield *ff);
	~Genetic();
	// Genetic move types
	enum GeneticMove { ReflectionMove, ExpansionMove, OuterContractionMove, InnerContractionMove, ShrinkMove, AllMoves, nGeneticMoves };


	/*
	// Basic Data
	*/
	private:
	// Number of genes in pool
	int nGenes_;
	// Static gene pool
	List<Alpha> genePool_;
	// Number of alpha values in each gene
	int nAlpha_;
	// Variation of initial and newborn genes from reference values
	double lengthScale_;
	// Best alpha encountered during algorithm
	Alpha bestAlpha_;
	// Target system
	TargetSystem *targetSystem_;
	// Target forcefield for minimisation
	Forcefield *forcefield_;
	// Flag specifying whether a better point was found by the Genetic (after minimisation)
	bool betterPointFound_;
	
	private:
	// Return (calculating if necessary) cost of specified gene
	double cost(int n);
	// Return (calculating if necessary) cost of supplied gene
	double cost(Alpha* gene);
	// Sort genes according to cost, with lowest first in array
	void sort();

	public:
	// Generate initial population
	void initialise(Alpha &initVertex, int populationSize, double baselengthscale);
	// Return MSD of cost values in current Genetic
	double costMSD();
	// Set specified gene
	void setGene(int n, Alpha geneInfo);
	// Return specified gene
	Alpha gene(int n);
	// Evolve current gene population by one generation
	Alpha evolve(double selectionThreshold, double mutateThreshold);
	// Return whether top two genes are converged within some threshold
	bool converged(double threshold, double deathRate);
	// Refresh pool with new genes, using the first as a basis
	void refresh();
	// Return whether a better point was found by the algorithm
	bool betterPointFound();
};

#endif
