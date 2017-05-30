/*
	*** Species Definition
	*** species.cpp
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

#include "species.h"
#include "constants.h"

// Constructor
Species::Species()
{
	// Private variables
	nMolecules_ = 0;
	firstAtom_ = 0;
	elecScalingMatrix_ = NULL;
	vdwScalingMatrix_ = NULL;

	// Public variables
	next = NULL;
	prev = NULL;
}

// Destructor
Species::~Species()
{
	if (elecScalingMatrix_ != NULL)
	{
		for (int n=0; n<atoms_.nItems(); ++n) delete[] elecScalingMatrix_[n];
		delete[] elecScalingMatrix_;
	}
	if (vdwScalingMatrix_ != NULL)
	{
		for (int n=0; n<atoms_.nItems(); ++n) delete[] vdwScalingMatrix_[n];
		delete[] vdwScalingMatrix_;
	}
}

/*
// Basic Data
*/

// Set name of species
void Species::setName(const char *s)
{
	name_ = s;
}

// Return name of species
const char *Species::name()
{
	return name_;
}

// Set number of molecules/occurrences of species
void Species::setNMolecules(int n)
{
	nMolecules_ = n;
}

// Return number of molecules/occurrences of species
int Species::nMolecules()
{
	return nMolecules_;
}

// Add atom
void Species::addAtom(const char *name)
{
	Atom *atom = atoms_.add();
	atom->setName(name);
}

// Return first atom in list
Atom *Species::atoms()
{
	return atoms_.first();
}

// Return name of nth atom
Atom *Species::atom(int index)
{
	if (index < 0 || index >= atoms_.nItems()) printf("Error: Atom index %i out of range for species '%s'\n", index, name_.get());
	else return atoms_[index];
	return NULL;
}

// Return number of defined atoms in single occurrence
int Species::nAtoms()
{
	return atoms_.nItems();
}

// Set index of first atom
void Species::setFirstAtom(int n)
{
	firstAtom_ = n;
}

// Return index of first atoms
int Species::firstAtom()
{
	return firstAtom_;
}

// Set mass of atom specified
bool Species::setMass(int n, double d)
{
	if ((n<0) || (n>=atoms_.nItems()))
	{
		printf("Error: Tried to set the mass of an atom whose internal index (%i) is out of range for this species (%s)\n", n, name_.get());
		return FALSE;
	}
	atoms_[n]->setMass(d);
	return TRUE;
}

// Return mass of atom specified
double Species::mass(int n)
{
	return atoms_[n]->mass();
}

/*
// Bound Interactions
*/

// Add bond interaction to species
void Species::addBond(int i, int j)
{
	Intra *ixn = bonds_.add();
	ixn->setAtoms(i,j);
}

// Return first bond interaction in species
Intra *Species::bonds()
{
	return bonds_.first();
}

// Add angle interaction to species
void Species::addAngle(int i, int j, int k)
{
	Intra *ixn = angles_.add();
	ixn->setAtoms(i,j,k);
}

// Return first angle interaction in species
Intra *Species::angles()
{
	return angles_.first();
}

// Add torsion interaction to species
void Species::addTorsion(int i, int j, int k, int l, double escale, double vscale)
{
	Intra *ixn = torsions_.add();
	ixn->setAtoms(i,j,k,l);
	ixn->setScaleFactors(escale, vscale);
}

// Return first torsion interaction in species
Intra *Species::torsions()
{
	return torsions_.first();
}

/*
// Scaling Matrices
*/

// Create scaling matrices
void Species::createScalingMatrices()
{
	if (elecScalingMatrix_ != NULL) printf("BUG: Electrostatic scaling matrix in species '%s' already exists.\n", name_.get());
	if (vdwScalingMatrix_ != NULL) printf("BUG: van der Waals scaling matrix in species '%s' already exists.\n", name_.get());
	int n, m;
	Intra *bound;
	elecScalingMatrix_ = new double*[atoms_.nItems()];
	vdwScalingMatrix_ = new double*[atoms_.nItems()];
	for (n=0; n<atoms_.nItems(); ++n)
	{
		elecScalingMatrix_[n] = new double[atoms_.nItems()];
		vdwScalingMatrix_[n] = new double[atoms_.nItems()];
		for (m=0; m<atoms_.nItems(); ++m) elecScalingMatrix_[n][m] = 1.0;
		for (m=0; m<atoms_.nItems(); ++m) vdwScalingMatrix_[n][m] = 1.0;
	}
	// Loop over bonds, angles, and torsions, poking in values as we go
	for (bound = bonds_.first(); bound != NULL; bound = bound->next)
	{
		elecScalingMatrix_[bound->i()][bound->j()] = 0.0;
		elecScalingMatrix_[bound->j()][bound->i()] = 0.0;
		vdwScalingMatrix_[bound->i()][bound->j()] = 0.0;
		vdwScalingMatrix_[bound->j()][bound->i()] = 0.0;
	}
	for (bound = angles_.first(); bound != NULL; bound = bound->next)
	{
		elecScalingMatrix_[bound->i()][bound->k()] = 0.0;
		elecScalingMatrix_[bound->k()][bound->i()] = 0.0;
		vdwScalingMatrix_[bound->i()][bound->k()] = 0.0;
		vdwScalingMatrix_[bound->k()][bound->i()] = 0.0;
	}
	for (bound = torsions_.first(); bound != NULL; bound = bound->next)
	{
		elecScalingMatrix_[bound->i()][bound->l()] = bound->elecScale();
		elecScalingMatrix_[bound->l()][bound->i()] = bound->elecScale();
		vdwScalingMatrix_[bound->i()][bound->l()] = bound->vdwScale();
		vdwScalingMatrix_[bound->l()][bound->i()] = bound->vdwScale();
	}
}

// Return electrostatics scaling factor betwen specified atom pair
double Species::elecScalingFactor(int i, int j)
{
	if (i < 0 || i >= atoms_.nItems()) printf("BUG: Index i (%i) is out of range for elecScalingMatrix\n", i);
	else if (j < 0 || j >= atoms_.nItems()) printf("BUG: Index j (%i) is out of range for elecScalingMatrix\n", j);
	else return elecScalingMatrix_[i][j];
	return 0.0;
}

// Return van der Waals scaling factor betwen specified atom pair
double Species::vdwScalingFactor(int i, int j)
{
	if (i < 0 || i >= atoms_.nItems()) printf("BUG: Index i (%i) is out of range for vdwScalingMatrix\n", i);
	else if (j < 0 || j >= atoms_.nItems()) printf("BUG: Index j (%i) is out of range for vdwScalingMatrix\n", j);
	else return vdwScalingMatrix_[i][j];
	return 0.0;
}
