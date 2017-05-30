/*
	*** Species Definition
	*** species.h
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

#ifndef FM3_SPECIES_H
#define FM3_SPECIES_H

#include "dnchar.h"
#include "list.h"
#include "intradef.h"
#include "atom.h"

// Atomic/Molecular Species in Target System
class Species
{
	public:
	// Constructor / Destructor
	Species();
	~Species();
	// List pointers
	Species *prev, *next;


	/*
	// Basic Data
	*/
	private:
	// Name of species
	Dnchar name_;
	// Number of molecules (repeat units) in species
	int nMolecules_;
	// List of atoms (names/masses)
	List<Atom> atoms_;
	// Local (i.e. 0-(N-1)) index of first species atom in configuration
	int firstAtom_;
	
	public:
	// Set name of species
	void setName(const char *s);
	// Return name of species
	const char *name();
	// Set number of molecules/occurrences of species
	void setNMolecules(int n);
	// Return number of molecules/occurrences of species
	int nMolecules();
	// Add atom
	void addAtom(const char *name);
	// Return first atom in list
	Atom *atoms();
	// Return nth atom
	Atom *atom(int index);
	// Return number of defined atoms in single occurrence
	int nAtoms();
	// Set index of first atom
	void setFirstAtom(int n);
	// Return index of first atom
	int firstAtom();
	// Set mass of atom specified
	bool setMass(int n, double d);
	// Return mass of atom specified
	double mass(int n);


	/*
	// Bound Interactions
	*/
	private:
	// List of bonds between atoms (given as integer pairs)
	List<Intra> bonds_;
	// List of angles between atoms (given as integer triples)
	List<Intra> angles_;
	// List of torsions between atoms (given as integer quartets)
	List<Intra> torsions_;

	public:
	// Add bond interaction to species
	void addBond(int i, int j);
	// Return first bond interaction in species
	Intra *bonds();
	// Add angle interaction to species
	void addAngle(int i, int j, int k);
	// Return first angle interaction in species
	Intra *angles();
	// Add torsion interaction to species
	void addTorsion(int i, int j, int k, int l, double escale, double vscale);
	// Return first torsion interaction in species
	Intra *torsions();


	/*
	// Scaling Matrices
	*/
	private:
	// Scaling matrices for vdW and electrostatic interactions between atoms in the same molecule
	double **elecScalingMatrix_, **vdwScalingMatrix_;
	
	public:
	// Create scaling matrices
	void createScalingMatrices();
	// Return electrostatic scaling factor betwen specified atom pair
	double elecScalingFactor(int i, int j);
	// Return van der Waals scaling factor betwen specified atom pair
	double vdwScalingFactor(int i, int j);
};

#endif
