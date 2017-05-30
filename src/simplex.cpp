/*
	*** Simplex Definition and Methods
	*** simplex.cpp
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

#include "simplex.h"
#include "parameter.h"
#include "constants.h"
#include "ff.h"
#include "system.h"
#include "random.h"
#include "alphastore.h"
#include <string.h>

// Constructor
Simplex::Simplex(TargetSystem* system, Forcefield* ff)
{
	// Private variables
	nAlpha_ = 0;
	nVertices_ = 0;
	forcefield_ = ff;
	targetSystem_ = system;
	vertices_ = NULL;
	betterPointsFound_ = 0;

	// Set move parameters
	rho_ = 1.0;
	chi_ = 2.0;
	gamma_ = 0.5;
	sigma_ = 0.5;
}

Simplex::~Simplex()
{
	if (vertices_ != NULL) delete[] vertices_;
}

// Return (calculating if necessary) cost of specified vertex
double Simplex::cost(int n)
{
	// Compute cold cost of vertex first, to compare against bestAlpha, then add random T fluctuation
	double coldCost = vertices_[n].cost(targetSystem_, forcefield_);
	// Check best point in *this* step
	if (coldCost < bestAlpha_.cost(targetSystem_, forcefield_))
	{
		bestAlpha_ = vertices_[n];
		++betterPointsFound_;
	}
	// Check best overall point stored
	targetSystem_->updateBestAlpha(&vertices_[n]);
	return coldCost;
}

// Return (calculating if necessary) cost of supplied vertex
double Simplex::cost(Alpha &vertex)
{
	// Compute cold cost of vertex first, to compare against bestAlpha, then add random T fluctuation
	double coldCost = vertex.cost(targetSystem_, forcefield_);
	// Check best point in *this* step
	if (coldCost < bestAlpha_.cost(targetSystem_, forcefield_))
	{
		bestAlpha_ = vertex;
		++betterPointsFound_;
	}
	// Check best overall point stored
	targetSystem_->updateBestAlpha(&vertex);
	return coldCost;
}

// Return whether to accept  move based on supplied cost value, vertex, and temperature
bool Simplex::accept(double trialCost, int vertexId, double temperature)
{
	// First, get cost of comparison vertex
	double comparisonCost = vertices_[vertexId].cost(targetSystem_, forcefield_);
	// If supplied cost is lower, then accept it regardless
	if (trialCost < comparisonCost) return TRUE;
	else
	{
		// Accept with some probability...
		double deltaCost = trialCost - comparisonCost;
		if (Random::number() < exp(-deltaCost/temperature))
		{
// 			printf("Accepted bolzmann\n");
			return TRUE;
		}
	}
	return FALSE;
}

// Find extreme cost values in current Simplex at specified vertex temperature
void Simplex::findExtremes(double temperature)
{
	// Set initial values
	vBest_ = 0;
	vWorst_ = 1;
	vNextWorst_ = 2;
	double thisCost;
	int n;
	// Find best, worst, and nextWorst vertices as normal...
	for (n=0; n<nVertices_; ++n)
	{
		thisCost = cost(n);
		if (cost(n) < cost(vBest_)) vBest_ = n;
		else if (cost(n) > cost(vWorst_))
		{
			vNextWorst_ = vWorst_;
			vWorst_ = n;
		}
		else if ((cost(n) > cost(vNextWorst_)) && n != vBest_) vNextWorst_ = n;
	}
	
	// ...and then randomly select points in the Simplex which we will trial as new best, worst, and nextworst points
	do
	{
		n = int(Random::number()*nVertices_);
	} while (n == vWorst_);
	if (Random::number() < exp(-fabs(cost(vBest_)-cost(n))/temperature))
	{
		// Swap points if necessary
		if (vNextWorst_ == n) vNextWorst_ = vBest_;
// 		else if (vWorst_ == n) vWorst_ = vBest_;
		vBest_ = n;
	}
	n = int(Random::number()*nVertices_);
	if (Random::number() < exp(-fabs(cost(vWorst_)-cost(n))/temperature))
	{
		// Swap points if necessary
		if (vNextWorst_ == n) vNextWorst_ = vWorst_;
		else if (vBest_ == n) vBest_ = vWorst_;
		vWorst_ = n;
	}
	do
	{
		n = int(Random::number()*nVertices_);
	} while (n == vWorst_);
	if (Random::number() < exp(-fabs(cost(vNextWorst_)-cost(n))/temperature))
	{
		// Swap points if necessary
		if (vBest_ == n) vBest_ = vNextWorst_;
// 		else if (vWorst_ == n) vWorst_ = vNextWorst_;
		vNextWorst_ = n;
	}
}

// Reflect specified vertex through specified centroid
Alpha Simplex::reflect(Alpha &centroid)
{
	// Reflect specified point : x(new) = (1 + rho) * y - rho * x(worst)    where y = sum x(i) / n, i != worst
	Alpha newVertex(centroid);
	newVertex = centroid * (1.0 + rho_) - vertices_[vWorst_] * rho_;
// 		xnew(:) = (1.0 + rho) * xbar(:) - rho * vertex(nalpha+1,:)
	return newVertex;
}

// Expand Simplex about worst point
Alpha Simplex::expand(Alpha& centroid)
{
	// Expand around highest point : x(new) = (1 + rho * chi) * xbar - rho * chi * x(n+1)    where xbar = sum(i=1,n) x(i) / n
	Alpha newVertex(centroid);
	newVertex = centroid * (1.0 + rho_ * chi_) - vertices_[vWorst_] * rho_ * chi_;
//	xnew = (1.0 + rho * chi) * xbar - rho * chi * vertex(nalpha+1,:)
	return newVertex;
}

// Contract Simplex, giving new vertex outside current polytope
Alpha Simplex::contractOutside(Alpha& centroid)
{
	// Contract simplex (new point outside current polytope) : x(new) = (1 + rho * gamma) * xbar - rho * gamma * vertex(nalpha+1)
	Alpha newVertex(centroid);
	newVertex = centroid * (1.0 + rho_ * gamma_) - vertices_[vWorst_] * (rho_ * gamma_);
	return newVertex;
}

// Contract Simplex, giving new vertex inside current polytope
Alpha Simplex::contractInside(Alpha& centroid)
{
	// Contract simplex (new point inside current polytope) : x(new) = (1 - gamma) * xbar + gamma * vertex(nalpha+1)
	Alpha newVertex(centroid);
	newVertex = centroid * (1.0 - gamma_) - vertices_[vWorst_] * gamma_;
	return newVertex;
}

// Shrink Simplex, contracting around best point (and leaving it as-is)
void Simplex::shrink()
{
	// Shrink vertices of current simplex, leaving only the first (x(1)) as-is: x(n) = x(1) + sigma(x(n) - x(1)),   n=2,nalpha+1
	for (int n=0; n<nVertices_; ++n)
	{
		if (n == vBest_) continue;
		vertices_[n] = (vertices_[n] - vertices_[vBest_]) * sigma_ + vertices_[vBest_];
	}
}

// Initialise starting Simplex in standard way, nudging params
void Simplex::initialise(Alpha &initVertex, double baselengthscale)
{
	nAlpha_ = initVertex.nAlpha();
	nVertices_ = nAlpha_+1;
	bestAlpha_ = initVertex;
	betterPointsFound_ = 0;
	baseLengthScale_ = baselengthscale;

	vertices_ = new Alpha[nVertices_];
	vertexStore_ = new AlphaStore[nVertices_];
	vertices_[0] = initVertex;

	double lscale, r;
	for (int n=1; n<nVertices_; ++n)
	{
		vertices_[n] = vertices_[0];
		r = (2.0*Random::number()) - 1.0;
		for (int m=0; m<nAlpha_; ++m)
		{
			lscale = baseLengthScale_ * vertices_[0].lengthScale(m);
			vertices_[n].multiplyAlpha(m, 1.0+lscale*r);
		}
	}
}

// Initialise simplex by adding one vertex at a time...
void Simplex::addVertex(Alpha &vertex)
{
	// If this is the first vertex it will determine the number of vertices we expect from its alpha population
	if (nVertices_ == 0)
	{
		nAlpha_ = vertex.nAlpha();
		vertices_ = new Alpha[nAlpha_+1];
		vertexStore_ = new AlphaStore[nAlpha_+1];
		vertices_[0] = vertex;
		bestAlpha_ = vertex;
		nVertices_ = 1;
	}
	else if (nVertices_ == (nAlpha_+1))
	{
		printf("BUG: Tried to add too many vertices to Simplex\n");
		return;
	}
	else
	{
		vertices_[nVertices_] = vertex;
		++nVertices_;
	}
}
	
// Return MSD of cost values in current Simplex
double Simplex::costMSD()
{
	// Get average cost of all vertices
	double average = 0.0;
	for (int n=0; n<nVertices_; ++n) average += cost(n);
	average /= nVertices_;
	// Now calculate SD of individual costs with mean value
	double serror = 0.0;
	for (int n=0; n<nVertices_; ++n) serror += (cost(n) - average)*(cost(n) - average);
	return sqrt(serror/nVertices_);
}

// Return whether Simplex has converged
bool Simplex::converged()
{
	return (costMSD() < simplexTolerance_);
}

// Print vertex information
void Simplex::printVertexInformation()
{
	for (int n=0; n<nVertices_; ++n) printf("[SMPLX]\t\tVertex %i has cold cost value = %12.6e\n", n, cost(n));
}

// Print Simplex move information
void Simplex::printMoveInformation()
{
	printf("[SMPLX]\tIn %i moves %i reflections, %i expansions, %i contractions (%i inner, %i outer) and %i shrinks were performed.\n", moveCount_[Simplex::AllMoves], moveCount_[Simplex::ReflectionMove], moveCount_[Simplex::ExpansionMove], moveCount_[Simplex::OuterContractionMove]+moveCount_[Simplex::InnerContractionMove], moveCount_[Simplex::InnerContractionMove], moveCount_[Simplex::OuterContractionMove], moveCount_[Simplex::ShrinkMove]);
}

// Minimise simplex for supplied configuration and forcefield
Alpha Simplex::minimise(int maxMoves, double tolerance, double simplexTemperature)
{
	// Simplex Simulated Annealing following original Nelder-Mead simplex algorithm
	// and incoporating Simulated Annealing. In all the following a positive temperature
	// implies a subtraction from the total cost function of the vertex (i.e. reducing costs)
	// while a negative temperature implies an addition to the total cost function (i.e. increasing costs)
	
	
	int n, move;
	double reflectedCost, newCost;
	Alpha centroid, reflectedVertex, newVertex;
	
	// Check nAlpha
	if (nAlpha_ == 0)
	{
		printf("Error: Simplex appears not to have been initialised properly.\n");
		targetSystem_->flagAbort();
		return newVertex;
	}

	nCycles_ = 0;
	nMoves_ = maxMoves;
	simplexTolerance_ = tolerance;
	betterPointsFound_ = 0;

	// Reset move counters
	for (n=0; n<Simplex::nSimplexMoves; ++n) moveCount_[n] = 0;
	
	// Begin Simplex Minimisation
	for (move=1; move<=nMoves_; ++move)
	{
		++moveCount_[Simplex::AllMoves];

		// Find best, worst, and next-worst points of Simplex
		findExtremes(simplexTemperature);
// 		printf("SIMPLEX : "); for (n=0; n<nVertices_; ++n) printf("%c%12.5e ", vBest_ == n ? 'B' : (vWorst_ == n ? 'W' : (vNextWorst_ == n ? 'N' : ' ')), vertices_[n].cost(targetSystem_, forcefield_)); printf("\n");

		// Calculate centroid for use in simplex move routines
		centroid.zero();
		for (n=0; n<nVertices_; ++n) if (n != vWorst_) centroid += vertices_[n];
		centroid /= nAlpha_;
	
		// First move attempt - Calculate reflection vertex
		reflectedVertex = reflect(centroid);
		reflectedCost = cost(reflectedVertex);
		
		// ... If this point is lower in cost than the next-highest vertex, but higher than the best point, accept it and terminate iteration
// 		if (reflectedCost < cost(vNextWorst_,-temperature) && reflectedCost > cost(vBest_,-temperature))
		if (accept(reflectedCost, vNextWorst_, simplexTemperature) && (!accept(reflectedCost, vBest_, simplexTemperature)))
		{
			vertices_[vWorst_] = reflectedVertex;
			++moveCount_[Simplex::ReflectionMove];
			if (converged()) break;
			continue;
		}
	
		// Stage 3a) - Calculate expansion point (if reflected point is lower than the best point already found)
		// Stage 3b) - Calculate contraction point (if reflected point is worse than the next worst (n-1'th) point)
// 		if (reflectedCost < cost(vBest_, -temperature))
		if (accept(reflectedCost, vBest_, simplexTemperature))
		{
			// Calculate expansion point
			newVertex = expand(centroid);
			newCost = cost(newVertex);
			
			// ... If expanded point is better than reflected point, accept expanded point and terminate iteration.
			// ... If reflected point is better (or equal to) expanded point, accept reflected point and terminate iteration
			if (newCost < reflectedCost)
			{
		// 	        write(0,*) "Expansion is better than reflection, so accepting expansion."
				vertices_[vWorst_] = newVertex;
				++moveCount_[Simplex::ExpansionMove];
				if (converged()) break;
				continue;
			}
			else
			{
		// 	        write(0,*) "Expansion is worse than reflection, so accepting reflection."
				vertices_[vWorst_] = reflectedVertex;
				++moveCount_[Simplex::ReflectionMove];
				if (converged()) break;
				continue;
			}
		}
		else
		{
			// Attempt contractions... if either fails, perform a shrink and continue with next iteration
			// ... If reflected point is bettern than worst point, but worse than next-best (n-1) point, do 'outside' contraction
// 			if ((reflectedCost >= cost(vNextWorst_, -temperature)) && (reflectedCost < cost(vWorst_, -temperature)))
			if (!accept(reflectedCost, vNextWorst_, simplexTemperature) && accept(reflectedCost, vWorst_, simplexTemperature))
			{
				newVertex = contractOutside(centroid);
				newCost = cost(newVertex);
				if (newCost <= reflectedCost)
				{
			// 	  	write(0,*) "Contraction(outside) is better than reflection, so accepting contraction."
					vertices_[vWorst_] = newVertex;
					++moveCount_[Simplex::OuterContractionMove];
					if (converged()) break;
					continue;
				}
				else
				{
					// Perform shrink
					shrink();
					++moveCount_[Simplex::ShrinkMove];
					if (converged()) break;
					continue;
				}
			}
			else
			{
				// ...otherwise do an inside contraction
				newVertex = contractInside(centroid);
				newCost = cost(newVertex);
				// !write(0,"('Ci ',e10.2,40f6.2)") cost_c,vert_c
// 				if (newCost < cost(vWorst_, -temperature))
				if (accept(newCost, vWorst_, simplexTemperature))
				{
					vertices_[vWorst_] = newVertex;
					++moveCount_[Simplex::InnerContractionMove];
					if (converged()) break;
					continue;
				}
				else
				{
					// Perform shrink
					shrink();
					++moveCount_[Simplex::ShrinkMove];
					if (converged()) break;
					continue;
				}
			}
		}
	}

// 	findExtremes(0.0);
	return bestAlpha_;
}

// Minimise simplex for supplied configuration and forcefield
Alpha Simplex::minimiseTGAY(int maxCycles, int movesPerCycle, double tolerance, double simplexTemperature, double vertexTemperature)
{
	// Simplex Simulated Annealing following original Nelder-Mead simplex algorithm
	// and incoporating Simulated Annealing. In all the following a positive temperature
	// implies a subtraction from the total cost function of the vertex (i.e. reducing costs)
	// while a negative temperature implies an addition to the total cost function (i.e. increasing costs)

	int a, n, m, cycle, move, nFailed, lastBest, failLimit = 1e9;
	double reflectedCost, newCost;
	double tSimplex, tVertex, tMultiplier, tDelta;
	Alpha centroid, reflectedVertex, newVertex;
	
	// Check nAlpha
	if (nAlpha_ == 0)
	{
		printf("Error: Simplex appears not to have been initialised properly.\n");
		targetSystem_->flagAbort();
		return newVertex;
	}

	nCycles_ = maxCycles;
	nMoves_ = movesPerCycle;
	simplexTolerance_ = tolerance;
	betterPointsFound_ = 0;

	// Reset move counters
	for (n=0; n<Simplex::nSimplexMoves; ++n) moveCount_[n] = 0;
	
	// Make copy of input vertices (only used with hot vertices)
	Alpha coldVertices_[nVertices_];
	for (n=0; n<nVertices_; ++n) coldVertices_[n] = vertices_[n];
	
	// Store starting temperatures
	tDelta = 1.0 / nCycles_;
	nFailed = 0;
	lastBest = -1;
	
	// Begin Simplex Minimisation
	for (cycle=1; cycle<=nCycles_; ++cycle)
	{
		// Determine cycle temperatures
		tMultiplier = targetSystem_->temperatureMultiplier((cycle-1)*tDelta);
		tSimplex = tMultiplier * simplexTemperature;
		tVertex = tMultiplier * vertexTemperature;
		
		printf("[SMPLX]\tSimplex temperature = %12.5e, vertex temperature = %f\n", tSimplex, tVertex);
		findExtremes(0.0);

		// Reinitialise Simplex?
		if (cycle%failLimit == 0)
		{
			printf("Progress is slow - reinitialising Simplex vertices...\n");
			vertices_[0] = bestAlpha_;
// 				for (n=1; n<nVertices_; ++n) vertices_[n] = vertexStore_[n].predict(targetSystem_, forcefield_);
// 				for (n=0; n<nVertices_; ++n)
// 				{
// 					coldVertices_[n] = vertices_[n];
// 					vertexStore_[n].clear();
// 				}
			double lscale, r;
			for (n=1; n<nVertices_; ++n)
			{
				vertices_[n] = vertices_[0];
				r = Random::number() - 0.5;
				for (int m=0; m<nAlpha_; ++m)
				{
// 					lscale = baseLengthScale_ * vertices_[0].lengthScale(m);
					lscale = 1.0;
					vertices_[n].multiplyAlpha(m, 1.0+lscale*r);
				}
				coldVertices_[n] = vertices_[n];
				vertexStore_[n].clear();
			}
			nFailed = 0;
		}
		
		for (move = 0; move < movesPerCycle; ++move)
		{
			++moveCount_[Simplex::AllMoves];

			// Randomise each vertex a little (from its coldVertex copy) according to the supplied vertex temperature
			for (n=0; n<nVertices_; ++n)
			{
				// Can we predict a better point from past experience?
				if (vertexStore_[n].nData() > 0)
				{
					newVertex = vertexStore_[n].predict(targetSystem_, forcefield_);
// 					printf("Predicted vertex has cost of %11.5e\n", cost(newVertex));
					if (cost(newVertex) < cost(coldVertices_[n]))
					{
						vertices_[n] = newVertex;
						vertexStore_[n].add(&vertices_[n], targetSystem_, forcefield_);
						continue;
					}
				}
				
				vertices_[n] = coldVertices_[n];
				for (m=0; m<nAlpha_; ++m) vertices_[n].multiplyAlpha(m, 1.0+(2.0*Random::number()-1.0)*tVertex);
				cost(vertices_[n]);
// 				printf("Vertex %i : ", n); for (m=0; m<vertices_[n].nAlpha(); ++m) printf("%10.5e ", vertices_[n].alpha(m)); printf("\n");

				vertexStore_[n].add(&vertices_[n], targetSystem_, forcefield_);
			}

			// Find best, worst, and next-worst points of Simplex
			findExtremes(tSimplex);

			// Calculate centroid for use in simplex move routines
			centroid.zero();
			for (n=0; n<nVertices_; ++n) if (n != vWorst_) centroid += vertices_[n];
			centroid /= nAlpha_;
	
			// First move attempt - Calculate reflection vertex
			reflectedVertex = reflect(centroid);
			reflectedCost = cost(reflectedVertex);
			
			// ... If this point is lower in cost than the next-highest vertex, but higher than the best point, accept it and terminate iteration
	// 		if (reflectedCost < cost(vNextWorst_,-temperature) && reflectedCost > cost(vBest_,-temperature))
			if (accept(reflectedCost, vNextWorst_, tSimplex) && (!accept(reflectedCost, vBest_, tSimplex)))
			{
				vertices_[vWorst_] = reflectedVertex;
				coldVertices_[vWorst_] = reflectedVertex;
// 				vertexStore_[vWorst_].clear();
				++moveCount_[Simplex::ReflectionMove];
				++nFailed;
				if (converged()) break;
				continue;
			}
		
			// Stage 3a) - Calculate expansion point (if reflected point is lower than the best point already found)
			// Stage 3b) - Calculate contraction point (if reflected point is worse than the next worst (n-1'th) point)
			if (accept(reflectedCost, vBest_, tSimplex))
			{
				// Reset nFailed counter, since we have found a better point
				nFailed = 0;

				// Calculate expansion point
				newVertex = expand(centroid);
				newCost = cost(newVertex);
				
				// ... If expanded point is better than reflected point, accept expanded point and terminate iteration.
				// ... If reflected point is better (or equal to) expanded point, accept reflected point and terminate iteration
				if (newCost < reflectedCost)
				{
			// 	        write(0,*) "Expansion is better than reflection, so accepting expansion."
					vertices_[vWorst_] = newVertex;
					coldVertices_[vWorst_] = newVertex;
// 					vertexStore_[vWorst_].clear();
					++moveCount_[Simplex::ExpansionMove];
					if (converged()) break;
					continue;
				}
				else
				{
			// 	        write(0,*) "Expansion is worse than reflection, so accepting reflection."
					vertices_[vWorst_] = reflectedVertex;
					coldVertices_[vWorst_] = reflectedVertex;
					++moveCount_[Simplex::ReflectionMove];
					if (converged()) break;
					continue;
				}
			}
			else
			{
				++nFailed;
				// Attempt contractions... if either fails, perform a shrink and continue with next iteration
				// ... If reflected point is bettern than worst point, but worse than next-best (n-1) point, do 'outside' contraction
				if (!accept(reflectedCost, vNextWorst_, tSimplex) && accept(reflectedCost, vWorst_, tSimplex))
				{
					newVertex = contractOutside(centroid);
					newCost = cost(newVertex);
					if (newCost <= reflectedCost)
					{
						vertices_[vWorst_] = newVertex;
						coldVertices_[vWorst_] = newVertex;
// 						vertexStore_[vWorst_].clear();
						++moveCount_[Simplex::OuterContractionMove];
						if (converged()) break;
						continue;
					}
					else
					{
						// Perform shrink
						shrink();
						for (n=0; n<nVertices_; ++n)
						{
							coldVertices_[n] = vertices_[n];
// 							vertexStore_[n].clear();
						}
						++moveCount_[Simplex::ShrinkMove];
						if (converged()) break;
						continue;
					}
				}
				else
				{
					// ...otherwise do an inside contraction
					newVertex = contractInside(centroid);
					newCost = cost(newVertex);
					if (accept(newCost, vWorst_, tSimplex))
					{
						vertices_[vWorst_] = newVertex;
						coldVertices_[vWorst_] = newVertex;
// 						vertexStore_[vWorst_].clear();
						++moveCount_[Simplex::InnerContractionMove];
						if (converged()) break;
						continue;
					}
					else
					{
						// Perform shrink
						shrink();
						for (n=0; n<nVertices_; ++n)
						{
							coldVertices_[n] = vertices_[n];
// 							vertexStore_[n].clear();
						}
						++moveCount_[Simplex::ShrinkMove];
						if (converged()) break;
						continue;
					}
				}
			}
		}
		
		// Print brief cycle cummary
		printf("[SMPLX]\tBest (cycle %i): %11.5e : ", cycle, vertices_[vBest_].cost(targetSystem_, forcefield_)); vertices_[vBest_].print();
	}

	return bestAlpha_;
}

// Return whether a better point was found by the Simplex
bool Simplex::betterPointFound()
{
	return (betterPointsFound_ > 0);
}
