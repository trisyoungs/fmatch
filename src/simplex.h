/*
	*** Simplex Definition and Method
	*** simplex.h
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

#ifndef FM3_SIMPLEX_H
#define FM3_SIMPLEX_H

#include "alpha.h"
#include "list.h"

// Forward Declarations
class Forcefield;
class TargetSystem;
class AlphaStore;

// Simplex Definition
class Simplex
{
	public:
	// Constructor / Destructor
	Simplex(TargetSystem *system, Forcefield *ff);
	~Simplex();
	// Simplex move types
	enum SimplexMove { ReflectionMove, ExpansionMove, OuterContractionMove, InnerContractionMove, ShrinkMove, AllMoves, nSimplexMoves };


	/*
	// Basic Data
	*/
	private:
	// Simplex move parameters
	double rho_, chi_, gamma_, sigma_;
	// Number of alpha values per vertex in simplex (i.e. size of parameterList_)
	int nAlpha_;
	// Number of vertices in vertex array
	int nVertices_;
	// Vertex array
	Alpha *vertices_;
	// Vertex store (used by Anneal2 method only)
	AlphaStore *vertexStore_;
	// Best alpha encountered during simplex
	Alpha bestAlpha_;
	// Target system
	TargetSystem *targetSystem_;
	// Target forcefield for minimisation
	Forcefield *forcefield_;
	// Indices of best, worst, and second worst vertices in current simplex
	int vBest_, vWorst_, vNextWorst_;
	// Maximum number of cycles to perform
	int nCycles_;
	// Number of moves per cycle to perform
	int nMoves_;
	// Convergence criterion for Simplex
	double simplexTolerance_;
	// Counter for Simplex moves
	int moveCount_[nSimplexMoves];
	// Integer count of number of better points found by the Simplex (after minimisation)
	int betterPointsFound_;
	// Base length scale of vertices
	double baseLengthScale_;
	
	private:
	// Return (calculating if necessary) cost of specified vertex
	double cost(int vertexId);
	// Return (calculating if necessary) cost of supplied vertex
	double cost(Alpha &vertex);
	// Return whether to accept  move based on supplied cost value, vertex, and temperature
	bool accept(double trialCost, int vertexId, double temperature);
	// Find extreme cost values in current Simplex
	void findExtremes(double temperature);
	// Reflect worst point in Simplex through centroid
	Alpha reflect(Alpha &centroid);
	// Expand simplex about worst point
	Alpha expand(Alpha &centroid);
	// Contract Simplex, giving new vertex outside current polytope
	Alpha contractOutside(Alpha &centroid);
	// Contract Simplex, giving new vertex inside current polytope
	Alpha contractInside(Alpha &centroid);
	// Shrink Simplex, contracting around best point (and leaving it as-is)
	void shrink();

	public:
	// Initialise starting Simplex in standard way, nudging params
	void initialise(Alpha &initVertex, double lscale);
	// Initialise simplex by adding one vertex at a time...
	void addVertex(Alpha &vertex);
	// Return MSD of cost values in current Simplex
	double costMSD();
	// Return whether Simplex has converged
	bool converged();
	// Print vertex information
	void printVertexInformation();
	// Print Simplex move information
	void printMoveInformation();
	// Perform standard simplex minimisation (at temperature specified)
	Alpha minimise(int maxCycles, double tolerance, double simplexTemperature);
	// Perform alternative Simplex annealing minimisation
	Alpha minimiseTGAY(int nCycles, int movesPerCycle, double tolerance, double simplexTemperature, double vertexTemperature);
	// Return whether a better point was found by the Simplex
	bool betterPointFound();
};

#endif
