/*
	*** Forcefield Definition
	*** ff.h
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

#ifndef FM3_FORCEFIELD_H
#define FM3_FORCEFIELD_H

#include "interaction.h"
#include "list.h"
#include "vector3.h"
#include "configuration.h"

// Forward Declarations
class Species;
class TargetSystem;
class AtomPair;
class Intra;

// Forcefield Definition
class Forcefield
{
	public:
	// Constructor / Destructor
	Forcefield();
	~Forcefield();


	/*
	// Main Forcefield Definitions
	// Forcefield terms covering all interactions defined in all frames
	*/
	private:
	// List of defined atomic charges
	List<Interaction> charges_;
	// List of defined atom-atom interactions
	List<Interaction> vdw_;
	// List of defined bond terms
	List<Interaction> bonds_;
	// List of defined angle terms
	List<Interaction> angles_;
	// List of defined torsion terms
	List<Interaction> torsions_;
	// Matrix of interactions between atom types
	Interaction ***interactionMatrix_;
	// Size of van der Waals interaction matrix
	int interactionMatrixSize_;
	// Pointer to dependent charge parameter
	Parameter *dependentCharge_;

	public:
	// Search for charge interaction
	Interaction *findCharge(const char *a1);
	// Search for van der Waals interaction (literal string comparison)
	Interaction *findVdw(const char *a1, const char *a2);
	// Search for bond interaction (literal string comparison)
	Interaction *findBond(const char *a1, const char *a2);
	// Search for angle interaction (literal string comparison)
	Interaction *findAngle(const char *a1, const char *a2, const char *a3);
	// Search for torsion interaction (literal string comparison)
	Interaction *findTorsion(const char *a1, const char *a2, const char *a3, const char *a4);
	// Search for van der Waals interaction (literal string comparison)
	Interaction *findVdw(Atom *a1, Atom *a2);
	// Search for bond interaction (literal string comparison)
	Interaction *findBond(Atom *a1, Atom *a2);
	// Search for angle interaction (literal string comparison)
	Interaction *findAngle(Atom *a1, Atom *a2, Atom *a3);
	// Search for torsion interaction (literal string comparison)
	Interaction *findTorsion(Atom *a1, Atom *a2, Atom *a3, Atom *a4);
	// Add charge interaction
	Interaction *addCharge();
	// Add vdw interaction
	Interaction *addVdw();
	// Add bond interaction
	Interaction *addBond();
	// Add angle interaction
	Interaction *addAngle();
	// Add torsion interaction
	Interaction *addTorsion();
	// Return list of specified van der Waals interactions
	Interaction *vdw();
	// Create van der Waals interaction matrix
	bool createInteractionMatrix(List<Atom> &uniqueTypes);
	// Return interaction matrix element
	Interaction *interactionMatrix(int i, int j);
	// Print interaction matrix
	void printInteractionMatrix();
	// Check any combinedVdw definitions
	bool checkCombinedVdw();
	// Set dependent charge
	bool setDependentCharge(Parameter *param);
	// Assign charge to supplied atom
	bool assignCharge(Atom *i, int nMolecules);
	// Search for and assign stored bond definition to supplied intra definition
	bool applyBond(Intra *bond, Atom *i, Atom *j);
	// Search for and assign stored angle definition to supplied intra definition
	bool applyAngle(Intra *angle, Atom *i, Atom *j, Atom *k);
	// Search for and assign stored torsion definition to supplied intra definition
	bool applyTorsion(Intra *torsion, Atom *i, Atom *j, Atom *k, Atom *l);
	// Return number of charge interactions defined
	int nCharges();
	// Return number of vdw interactions defined
	int nVdw();
	// Return number of bond interactions defined
	int nBonds();
	// Return number of angle interactions defined
	int nAngles();
	// Return number of torsion interactions defined
	int nTorsions();
	// Return actual value of charge interaction ID specified
	double charge(int id);


	/*
	// Interaction Parameters used in Forcefield
	*/
	private:
	// List of fixed parameters during optimisation
	List<Parameter> fixedParameters_;
	// List of dependent parameters during optimisation
	List<Parameter> dependentParameters_;
	// List of free (fitted) parameters during optimisation
	List<Parameter> freeParameters_;
	// Logpoint totals for contribution types
	unsigned int logPoints_[Configuration::nContributionTypes];

	public:
	// Add new interaction parameter
	Parameter *addParameter(Parameter::ParameterType type, Parameter::ParameterStrength strength);
	// Return number of fixed parameters
	int nFixedParameters();
	// Return number of dependent parameters
	int nDependentParameters();
	// Return number of free parameters
	int nFreeParameters();
	// Return List of free parameters
	List<Parameter> *freeParameters();
	// Print free parameters
	void printFreeParameters();
	// Print fixed parameters
	void printFixedParameters();
	// Print dependent parameters
	void printDependentParameters();
	// Return cost penalty for free parameters going outside their defined limits
	double freeParameterPenalty();
	// Calculate parameter logpoints
	void calculateLogPoints();
	// Return specified logpoint
	unsigned int logPoint(Configuration::ContributionType type);


	/*
	// Methods
	*/
	public:
	// Setup all contained forcefield terms, check for missing parameters etc.
	bool setup();
	// Check charge setup
	bool checkChargeSetup();
	// Calculate dependent parameters
	bool calculateDependencies();
};

#endif
