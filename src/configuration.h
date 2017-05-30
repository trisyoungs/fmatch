/*
	*** Atomic Configuration
	*** configuration.h
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

#ifndef FM3_CONFIGURATION_H
#define FM3_CONFIGURATION_H

#include "esppoint.h"
#include "dnchar.h"
#include "matrix.h"
#include "interaction.h"

// Forward Declarations
class Species;
class Forcefield;
class KVector;
class AtomPair;
class TargetSystem;

// Atomic Configuration
class Configuration
{
	public:
	// Constructor / Destructor
	Configuration();
	~Configuration();
	// List pointers
	Configuration *prev, *next;
	// Contribution types
	enum ContributionType { AngleContribution, BondContribution, TorsionContribution, VdwContribution, ChargeContribution, nContributionTypes };

	/*
	// Basic Data
	*/
	private:
	// Name of this configuration
	Dnchar name_;

	public:
	// Set name of configuration
	void setName(const char *name);
	// Return name of this configuration
	const char *name();


	/*
	// Atomic Data
	*/
	private:
	// Number of atoms in configuration
	int nAtoms_;
	// Array of atomic coordinates
	Vec3<double> *r_;
	// Array of atomic forces
	Vec3<double> *f_;

	public:
	// Set all atom data
	bool setAtomData(int id, Vec3<double> coordinates, Vec3<double> forces, double weight);
	// Set atom position data
	bool setAtomCoordinates(int id, Vec3<double> coordinates);
	// Set atom forces data
	bool setAtomForces(int id, Vec3<double> forces, double weight);
	// Return number of atoms in configuration
	int nAtoms();
	// Return coordinates of atom specified
	Vec3<double> r(int id);
	// Return array of atomic positions
	Vec3<double> *r();


	/*
	// Unit Cell
	*/
	private:
	// Unit cell type
	int cellType_;
	// Unit cell definition
	Matrix cell_;
	// Unit cell inverse
	Matrix inverseCell_;
	// Reciprocal unit cell
	Matrix reciprocalCell_;
	// Reciprocal cell volume
	double reciprocalVolume_;
	// Reciprocal cell cutoff for Ewald sum
	double reciprocalCutoff_;

	public:
	// Setup unit cell
	void setCell(Matrix &cell, int cellType);
	// Return minimum image vector between supplied atoms
	Vec3<double> mimd(int i, int j);
	// Return reciprocal cell of configuration
	Matrix &reciprocalCell();
	// Return reciprocal cell volume
	double reciprocalVolume();
	// Return reciprocal cell cutoff for Ewald sum
	double reciprocalCutoff();


	/*
	// Additional Fit Data
	*/
	private:
	// Energies to fit
	double fitEnergyValues_[nContributionTypes];
	// Whether to fit specific energies
	bool fitEnergyType_[nContributionTypes];
	// Penalty factor weights for energy fiting
	double fitEnergyWeights_[nContributionTypes];
	// Array of ESP points to fit
	List<ESPPoint> espPoints_;
	
	public:
	// Set energy component to fit
	void setFitEnergy(ContributionType type, double value, double weight = 1.0);
	// Add electrostatic potential point to fit
	void addPotentialPoint(Vec3<double> r, double energy, double weight = 1.0);


	/*
	// Calculated Energy and Forces
	*/
	private:
	// System energy terms
	double energies_[nContributionTypes];
	// Arrays of calculated atomic force contributions from various interaction types
	Vec3<double> **fContributions_;
	// Total atomic forces
	Vec3<double> *fCalc_;
	// Array of weights to use when calculating force SOSE
	double *fWeight_;
	// Integer 'logpoints' at which energy/forces were last calculated
	unsigned int logPoint_[nContributionTypes];

	public:
	// Increase specified energy term
	void addEnergy(ContributionType type, double e);
	// Return specified energy term
	double energy(ContributionType type);
	// Return total energy
	double totalEnergy();
	// Add calculated force to specified atom
	void addForce(ContributionType type, Vec3<double> &fvec, int i);
	// Subtract calculated force from specified atom
	void subtractForce(ContributionType type, Vec3<double> &fvec, int i);
	// Calculate total forces from individual contributions
	void sumForceContributions();
	// Calculate and return sum of squares error in forces (and any additional fitted quantities)
	double totalSOSE();
	// Return calculated force array
	Vec3<double> *fCalc();
	// Return specified logpoint
	unsigned int logPoint(ContributionType type);
	// Update specified logpoint
	void updateLogPoint(Forcefield* ff, Configuration::ContributionType type);


	/*
	// Precalculated Quantities
	*/
	private:
	// Minimum image distance matrix
	double **distanceMatrix_;
	// Minimum image vector matrix
	Vec3<double> **vectorMatrix_;
	// List of all interacting atom pairs in configuration
	AtomPair *atomPairs_;
	// Number of interacting pairs in list
	int nAtomPairs_;
	// List of all excluded atom pairs in configuration
	AtomPair *excludedPairs_;
	// Number of interacting pairs in list
	int nExcludedPairs_;
	// Arrays of reciprocal space vectors for Ewald sum
	Vec3<double> **rCos_, **rSin_;
	// Value of kmax used to create cos and sin arrays
	Vec3<int> kMax_;
	// Array of Ewald kvectors to calculate in most recently prepared configuration
	KVector *kVectors_;
	// Cos and sin vectors for each atom in each kvector, and sums of charge types
	double **atomCos_, **atomSin_, **typeCos_, **typeSin_;
	// Map of atoms to their respective charge types
	int *chargeMap_;
	// Number of kvectors in arrays
	int nKVectors_;

	public:
	// Update/create distance matrix
	void updateDistanceMatrix();
	// Return distance for supplied pair
	double distance(int i, int j);
	// Return mim vector for supplied pair
	Vec3<double> vector(int i, int j);
	// (Re)create atompairs arrays
	void createPairArrays(int natompairs, int nexcludedpairs);
	// Return number of interacting atom pairs
	int nAtomPairs();
	// Return interacting pair list
	AtomPair *atomPairs();
	// Return number of excluded atom pairs
	int nExcludedPairs();
	// Return excluded atom pair list
	AtomPair *excludedPairs();
	// (Re)create kvector arrays
	void createKVectorArrays(Vec3< int > kmax, double alpha, double elecConvert);
	// Update k-vector arrays
	void updateKVectorArrays();
	// Create static kvector/position sums for atom types
	void createStaticKVectorSums(Species *firstSp, int kMax, int nChargeTypes);
	// Return rCos array
	Vec3<double> **rCos();
	// Return rSin array
	Vec3<double> **rSin();
	// Return atomCos array
	double **atomCos();
	// Return atomSin array
	double **atomSin();
	// Return typeCos array
	double **typeCos();
	// Return typeSin array
	double **typeSin();
	// Return atom charge map
	int *chargeMap();
	// Return kVector array
	KVector *kVectors();
	// Return number of kVectors in array
	int nKVectors();
	// Return whether configuration has precalculated atom/type arrays for Ewald sum
	bool hasStaticKVectorSums();
	// Free precalculated data arrays
	void deletePrecalculatedArrays();


	/*
	// Embedded Atom Model Data
	*/
	private:
	// Screening function for each atom pair
	double **Sij_;
	// Derivative of screening function for each atom pair
	double **dSij_drij_;
	// Total background electron density for each atom
	double *rho_;
	// Derivative of background electron density for each atom pair
	double **drho_drij_;
	// Atomic electron densities between atom pairs
	double **rhoa_;
	// Derivative of atomic electron densities between atom pairs
	double **drhoa_drij_;
	// Value of embedding function for each atom
	double *Fi_, *dFi_drho_, **dFi_drhoa_;
	// Pair potential between atom pairs
	double **phi_;
	// Derivative of pair potential between atom pairs
	double **dphi_drij_;

	public:
	// Calculate EAM forces on atoms
	bool calculateEAMForces(Interaction *eamPotential);
	// Calculate EAM forces on atoms (using stored distance matrix)
	bool calculateEAMForcesQuick(Interaction *eamPotential);
	// Delete EAM-related arrays
	void deleteEAMArrays();


	/*
	// Methods
	*/
	public:
	// Initialise all necessary arrays
	bool initialiseAtomArrays(int natoms);
	// Delete existing atom arrays
	void deleteAtomArrays();
	// Replicate atom/force data
	bool replicate(Vec3<int> rep, Species *firstsp);
	// Calculate Ewald sum cutoff
	void calculateEwaldCutoff(Vec3<int> kmax);
	// Reset specified calculated force and energy contributions
	void reset(ContributionType type);
	// Print configuration information, including calculated forces/deltas and energies
	void print(TargetSystem *targetSystem);
	// Print energy information only
	void printEnergy();
};

#endif
