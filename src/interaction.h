/*
	*** Interaction Definition
	*** bound.h
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

#ifndef FM3_INTERACTION_H
#define FM3_INTERACTION_H

#define MAXPARAMETERS 13

#include "list.h"
#include "reflist.h"
#include "parameter.h"

// Forward Declarations
class Configuration;
class Atom;
class AtomPair;

// Interaction Data Definition
class InteractionData
{
	public:
	// Keyword name of the interaction
	const char *name;
	// Number of parameters in interaction
	int nParameters;
	// Parameter keywords
	const char *parameters[MAXPARAMETERS];
};

// Bound Interaction
class Interaction
{
	public:
	// Constructor / Destructor
	Interaction();
	~Interaction();
	// List pointers
	Interaction *prev, *next;


	/*
	// Atom Data
	*/
	private:
	// Atom type names to which the interaction applies
	Dnchar types_[4];
	
	private:
	// Return whether strings match, considering wildcard characters
	bool matchWild(const char *s1, const char *s2);
	
	public:
	// Set atom data for interaction
	void setAtoms(const char *a1, const char *a2 = NULL, const char *a3 = NULL, const char *a4 = NULL);
	// Check stored atom names against supplied names (literal name matching)
	bool matches(const char *a1, const char *a2 = NULL, const char *a3 = NULL, const char *a4 = NULL);
	// Check stored atoms against supplied atoms (literal name matching)
	bool matches(Atom *i, Atom *j, Atom *k = NULL, Atom *l = NULL);
	// Check stored atom names against supplied names (parameter matching)
	bool matchesWild(const char *a1, const char *a2, const char *a3 = NULL, const char *a4 = NULL);
	// Check stored atoms against supplied atoms (wildcard matching)
	bool matchesWild(Atom *i, Atom *j, Atom *k = NULL, Atom *l = NULL);
	// Return type name specified
	const char *type(int n);
	// Return formatted atom interaction list text
	const char *atomsText();


	/*
	// Form / Parameter Data
	*/
	public:
	// Type of interaction
	enum InteractionType { NoType, ChargeType, VdwType, CombinedVdwType, BondType, AngleType, TorsionType };
	// Allowable interaction forms and parameter lists
	enum InteractionForm { NoForm, ChargeForm, HarmonicForm, MorseForm, CosineForm, DoubleCosineForm, TripleCosineForm, HarmonicCosineForm, LennardJones126Form, EmbeddedAtomForm, nInteractionForms };
	static InteractionForm interactionForm(const char *s);
	// Combination types (for type_ == CombinedVdwType)
	enum CombinationType { Arithmetic, Geometric };
	static CombinationType combinationType(const char *s);

	private:
	// Integer ID of this interaction in its respective list
	int id_;
	// Type of this interaction
	InteractionType type_;
	// Form of this interaction
	InteractionForm form_;
	// References to parameters associated to interaction (in order added in input file)
	Reflist<Parameter,int> parameters_;
	// Array of pointers to parameter variables, in order specified in InteractionData
	Parameter **values_;
	// Pointers to interactions to combine (only if type_ == CombinedVdwType)
	Interaction *sourceTypeA_, *sourceTypeB_;
	// Number of times this interaction is used in the system
	int usageCount_;

	public:
	// Set ID of interaction
	void setId(int id);
	// Return ID of interaction
	int id();
	// Set type of interaction
	void setType(InteractionType type);
	// Return type of interaction
	InteractionType type();
	// Add parameter to list, specifying (and setting) functional form as well
	bool addParameter(InteractionForm form, Parameter *param);
	// Add all parameters relevant to the specified functional form to the list
	bool addAllParameters(Interaction::InteractionForm form, List<Parameter> &destinationList, Parameter::ParameterType type, Parameter::ParameterStrength strength);
	// Return functional form of interaction
	InteractionForm form();
	// Return functional form text of interaction
	const char *formText();
	// Check validity of defined parameters and finalise function-specific parameter list
	bool setup();
	// Return string containing formatted form/parameter data
	const char *info();
	// Set combination source interactions (only if type_ == CombinedVdwType)
	void setCombinationSource(Interaction *a, Interaction *b);
	// Perform parameter combination
	bool combine();
	// Return whether all parameters involved in the interaction are variable
	bool allParametersVariable();
	// Return number of parameters associated with interaction
	int nParameters();
	// Return list of parameters associated to interaction
	Parameter **parameters();
	// Return nth parameter in list
	Parameter *parameter(int n);
	// Increase usage count
	void increaseUsageCount(int delta = 1);
	// Return charge count
	int usageCount();


	/*
	// Energy / Force Calculation
	*/
	private:
	// Evaluate value and derivative of interaction at supplied value
	bool evaluate(double x, double &fx, double &df_dx);

	public:
	// Calculate bound energy and forces for provided atom pair
	bool calculate(Configuration* cfg, AtomPair *pair);
	// Calculate energy and forces using coordinates of atoms supplied
	bool calculate(Configuration *cfg, int i, int j, int k = -1, int l = -1, double scale = 1.0);
};

#endif
