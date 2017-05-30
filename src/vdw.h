/*
	*** Van Der Waals Definition / Calculation
	*** vdw.h
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

#ifndef FM3_VDW_H
#define FM3_VDW_H

// Van der Waals Definition
class VdwInteraction
{
	public:
	// Constructor / Destructor
	VdwInteraction();
	~VdwInteraction();
	// List pointers
	VdwInteraction *prev, *next;


	/*
	// Functional Form
	*/
	public:
	// Allowable van der Waals forms
	enum VdwForm { LennardJones126, LennardJonesAB, EmbeddedAtom, Spline, nVanDerWaalsForms };
	static VdwForm vdwForm(const char *s);

	private:
	// Form of this definition
	VdwForm form_;

	public:
	// Set the form of this definition
	VdwForm form_;


	/*
	// Parameter Data
	*/
	private:
	// References to parameters associated to interaction
	Reflist<Parameter> parameters_;

	public:
	// Add parameter to list
	void addParameter(Parameter *param);


	/*
	// Methods
	*/
	public:
	// Check validity of defined parameters
	bool isValid();
	// Calculate energy and forces using coordinates of atoms supplied
	double calculate(XX);
};

#endif
