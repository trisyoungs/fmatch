/*
	*** Forcefield Definition
	*** ff.cpp
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

#include "ff.h"
#include "atom.h"
#include "atompair.h"
#include "species.h"
#include "configuration.h"
#include "system.h"
#include <string.h>

// Constructor
Forcefield::Forcefield()
{
	// Private variables
	interactionMatrix_ = NULL;
	interactionMatrixSize_ = 0;
	dependentCharge_ = NULL;
}

// Destructor
Forcefield::~Forcefield()
{
	if (interactionMatrix_ != NULL)
	{
		for (int n=0; n<interactionMatrixSize_; ++n) delete[] interactionMatrix_[n];
		delete[] interactionMatrix_;
	}
}

/*
// Main Forcefield Definitions
// Forcefield terms covering all interactions defined in all frames
*/

// Search for charge interaction
Interaction *Forcefield::findCharge(const char *a1)
{
	for (Interaction *i = charges_.first(); i != NULL; i = i->next) if (i->matches(a1)) return i;
	return NULL;
}

// Search for bond interaction
Interaction *Forcefield::findVdw(const char *a1, const char *a2)
{
	for (Interaction *i = vdw_.first(); i != NULL; i = i->next) if (i->matches(a1,a2)) return i;
	return NULL;
}

// Search for bond interaction
Interaction *Forcefield::findBond(const char *a1, const char *a2)
{
	for (Interaction *i = bonds_.first(); i != NULL; i = i->next) if (i->matches(a1,a2)) return i;
	return NULL;
}

// Search for bond interaction
Interaction *Forcefield::findAngle(const char *a1, const char *a2, const char *a3)
{
	for (Interaction *i = angles_.first(); i != NULL; i = i->next) if (i->matches(a1,a2,a3)) return i;
	return NULL;
}

// Search for bond interaction
Interaction *Forcefield::findTorsion(const char *a1, const char *a2, const char *a3, const char *a4)
{
	for (Interaction *i = torsions_.first(); i != NULL; i = i->next) if (i->matches(a1,a2,a3,a4)) return i;
	return NULL;
}

// Search for bond interaction
Interaction *Forcefield::findVdw(Atom *a1, Atom *a2)
{
	return findVdw(a1->name(), a2->name());
}

// Search for bond interaction
Interaction *Forcefield::findBond(Atom *a1, Atom *a2)
{
	return findBond(a1->name(), a2->name());
}

// Search for bond interaction
Interaction *Forcefield::findAngle(Atom *a1, Atom *a2, Atom *a3)
{
	return findAngle(a1->name(), a2->name(), a3->name());
}

// Search for bond interaction
Interaction *Forcefield::findTorsion(Atom *a1, Atom *a2, Atom *a3, Atom *a4)
{
	return findTorsion(a1->name(), a2->name(), a3->name(), a4->name());
}

// Add charge interaction
Interaction *Forcefield::addCharge()
{
	Interaction *ixn = charges_.add();
	ixn->setId(charges_.nItems()-1);
	return ixn;
}

// Add vdw interaction
Interaction *Forcefield::addVdw()
{
	Interaction *ixn = vdw_.add();
	ixn->setId(vdw_.nItems()-1);
	return ixn;
}

// Add bond interaction
Interaction *Forcefield::addBond()
{
	Interaction *ixn = bonds_.add();
	ixn->setId(bonds_.nItems()-1);
	return ixn;
}

// Add angle interaction
Interaction *Forcefield::addAngle()
{
	Interaction *ixn = angles_.add();
	ixn->setId(angles_.nItems()-1);
	return ixn;
}

// Add torsion interaction
Interaction *Forcefield::addTorsion()
{
	Interaction *ixn = torsions_.add();
	ixn->setId(torsions_.nItems()-1);
	return ixn;
}

// Return list of specified van der Waals interactions
Interaction *Forcefield::vdw()
{
	return vdw_.first();
}

// Create van der Waals interaction matrix
bool Forcefield::createInteractionMatrix(List<Atom> &uniqueTypes)
{
	if (interactionMatrix_ != NULL) printf("BUG: Interaction matrix already exists.\n");
	int n, m;
	Interaction *ixn;
	interactionMatrixSize_ = uniqueTypes.nItems();
	interactionMatrix_ = new Interaction**[uniqueTypes.nItems()];
	for (n=0; n<uniqueTypes.nItems(); ++n)
	{
		interactionMatrix_[n] = new Interaction*[uniqueTypes.nItems()];
		for (m=0; m<uniqueTypes.nItems(); ++m) interactionMatrix_[n][m] = NULL;
	}
	// So, loop over atom types defined in the uniqueTypes list, searching for relevant defined interactions as we go...
	// Do diagonal (type1 == type2) elements first
	for (n=0; n<interactionMatrixSize_; ++n)
	{
		ixn = findVdw(uniqueTypes[n]->name(), uniqueTypes[n]->name());
		if (ixn == NULL)
		{
			printf("Error: No van der Waals self-interaction defined for type '%s'.\n", uniqueTypes[n]->name());
			return FALSE;
		}
		interactionMatrix_[n][n] = ixn;
	}
	for (n=0; n<interactionMatrixSize_-1; ++n)
	{
		for (m=n+1; m<interactionMatrixSize_; ++m)
		{
			ixn = findVdw(uniqueTypes[n]->name(), uniqueTypes[m]->name());
			if (ixn != NULL)
			{
				interactionMatrix_[n][m] = ixn;
				interactionMatrix_[m][n] = ixn;
			}
			else
			{
				printf("WARNING: No van der Waals interaction defined between types '%s' and '%s'.\n", uniqueTypes[n]->name(), uniqueTypes[m]->name());
			}
		}
	}
	return TRUE;
}

// Return interaction matrix element
Interaction *Forcefield::interactionMatrix(int i, int j)
{
	if (i < 0 || i >= interactionMatrixSize_) printf("BUG: Index i (%i) for interaction matrix is out of bounds.\n", i);
	else if (j < 0 || j >= interactionMatrixSize_) printf("BUG: Index j (%i) for interaction matrix is out of bounds.\n", j);
	else return interactionMatrix_[i][j];
	return NULL;
}

// Print interaction matrix
void Forcefield::printInteractionMatrix()
{
	// Get column headers from diagonal matrix elements
	printf("\t         ");
	for (int n=0; n<interactionMatrixSize_; ++n) printf("%-8s ", interactionMatrix_[n][n]->type(0));
	printf("\n");
	for (int n=0; n<interactionMatrixSize_; ++n)
	{
		printf("\t%-8s ", interactionMatrix_[n][n]->type(0));
		for (int m=0; m<interactionMatrixSize_; ++m)
		{
			if (interactionMatrix_[n][m] == NULL) printf("--       ");
			else if (interactionMatrix_[n][m]->type() == Interaction::CombinedVdwType) printf("*%-7s ", interactionMatrix_[n][m]->formText());
			else printf("%-8s ", interactionMatrix_[n][m]->formText());
		}
		printf("\n");
	}
}

// Check any combinedVdw definitions
bool Forcefield::checkCombinedVdw()
{
	Interaction *ixn, *a, *b;
	for (Interaction *ixn = vdw_.first(); ixn != NULL; ixn = ixn->next)
	{
		if (ixn->type() != Interaction::CombinedVdwType) continue;

		// Check that both van der Waals self-interactions have been defined
		a = findVdw(ixn->type(0), ixn->type(0));
		if (a == NULL)
		{
			printf("Error: Combined van der Waals between '%s' and '%s' requested, but self-interaction for '%s' has not been defined.\n", ixn->type(0), ixn->type(1), ixn->type(0));
			return FALSE;
		}
		b = findVdw(ixn->type(1), ixn->type(1));
		if (b == NULL)
		{
			printf("Error: Combined van der Waals between '%s' and '%s' requested, but self-interaction for '%s' has not been defined.\n", ixn->type(0), ixn->type(1), ixn->type(1));
			return FALSE;
		}

		// Check that both existing interactions are of the same type
		if (a->form() != b->form())
		{
			printf("Error: Types specified for combination ('%s' and '%s') are of differing functional form ('%s' and '%s').\n", ixn->type(0), ixn->type(1), a->formText(), b->formText());
			return FALSE;
		}
		
		// Set source interaction types and create parameters in definition
		ixn->setCombinationSource(a,b);
		if (!ixn->addAllParameters(a->form(), dependentParameters_, Parameter::DependentParameter, Parameter::Weak)) return FALSE;
		printf("\tCreated and added parameters to interaction for combined van der Waals interaction between '%s' and '%s'\n", ixn->type(0), ixn->type(1));
	}
	return TRUE;
}

// Set dependent charge
bool Forcefield::setDependentCharge(Parameter *param)
{
	if (dependentCharge_ == NULL) dependentCharge_ = param;
	else return FALSE;
	return TRUE;
}

// Assign charge to supplied atom
bool Forcefield::assignCharge(Atom *i, int nMolecules)
{
	for (Interaction *ixn = charges_.first(); ixn != NULL; ixn = ixn->next) if (ixn->matches(i->name()))
	{
		i->setChargeParameter(ixn->parameter(0));
		ixn->increaseUsageCount(nMolecules);
		return TRUE;
	}
	return FALSE;
}

// Search for and assign stored bond definition to supplied intra definition
bool Forcefield::applyBond(Intra *bond, Atom *i, Atom *j)
{
	// Search for (wildcard) matching names in defined interactions
	for (Interaction *ixn = bonds_.first(); ixn != NULL; ixn = ixn->next) if (ixn->matchesWild(i,j))
	{
		bond->setInteraction(ixn);
		ixn->increaseUsageCount();
		return TRUE;
	}
	return FALSE;
}

// Search for and assign stored angle definition to supplied intra definition
bool Forcefield::applyAngle(Intra *angle, Atom *i, Atom *j, Atom *k)
{
	// Search for (wildcard) matching names in defined interactions
	for (Interaction *ixn = angles_.first(); ixn != NULL; ixn = ixn->next) if (ixn->matchesWild(i,j,k))
	{
		angle->setInteraction(ixn);
		ixn->increaseUsageCount();
		return TRUE;
	}
	return FALSE;
}

// Search for and assign stored torsion definition to supplied intra definition
bool Forcefield::applyTorsion(Intra *torsion, Atom *i, Atom *j, Atom *k, Atom *l)
{
	// Search for (wildcard) matching names in defined interactions
	for (Interaction *ixn = torsions_.first(); ixn != NULL; ixn = ixn->next) if (ixn->matchesWild(i,j,k,l))
	{
		torsion->setInteraction(ixn);
		ixn->increaseUsageCount();
		return TRUE;
	}
	return FALSE;
}

// Return number of charge interactions defined
int Forcefield::nCharges()
{
	return charges_.nItems();
}

// Return number of vdw interactions defined
int Forcefield::nVdw()
{
	return vdw_.nItems();
}

// Return number of bond interactions defined
int Forcefield::nBonds()
{
	return bonds_.nItems();
}

// Return number of angle interactions defined
int Forcefield::nAngles()
{
	return angles_.nItems();
}

// Return number of torsion interactions defined
int Forcefield::nTorsions()
{
	return torsions_.nItems();
}

// Return actual value of charge interaction ID specified
double Forcefield::charge(int id)
{
	if ((id < 0) || (id >= charges_.nItems())) printf("Internal Error: ID is out of range for retrieval of charge.\n");
	else return charges_[id]->parameter(0)->value();
	return 0.0;
}

/*
// Interaction Parameters used in Forcefield
*/

// Add new interaction parameter
Parameter *Forcefield::addParameter(Parameter::ParameterType type, Parameter::ParameterStrength strength)
{
	Parameter *parameter;
	if (type == Parameter::FixedParameter) parameter = fixedParameters_.add();
	else if (type == Parameter::FreeParameter) parameter = freeParameters_.add();
	else parameter = dependentParameters_.add();
	parameter->setType(type, strength);
	return parameter;
}

// Return number of fixed parameters
int Forcefield::nFixedParameters()
{
	return fixedParameters_.nItems();
}

// Return number of dependent parameters
int Forcefield::nDependentParameters()
{
	return dependentParameters_.nItems();
}

// Return number of free parameters
int Forcefield::nFreeParameters()
{
	return freeParameters_.nItems();
}

// Return List of free parameters
List<Parameter> *Forcefield::freeParameters()
{
	return &freeParameters_;
}

// Print free parameters
void Forcefield::printFreeParameters()
{
	Parameter *param;
	int n = 0;
	for (param = freeParameters_.first(); param != NULL; param = param->next)
	{
// 		if (charges_.contains((Charge*)param)) printf("  Parameter %4i\t    charge %12.5e (on atom %s)\n", ++n, param->value(), param->name());
		printf("  Parameter %4i\t%10s %12.5e (in %s interaction between atoms %s)\n", ++n, param->name(), param->value(), param->parentForm(), param->parentAtoms());
		if (param->hasLimits()) printf("\t\t--> limits are { %s, %s }\n", param->lowerLimitAsText(), param->upperLimitAsText());
	}
}

// Print fixed parameters
void Forcefield::printFixedParameters()
{
	Parameter *param;
	int n = 0;
	for (param = fixedParameters_.first(); param != NULL; param = param->next)
	{
// 		if (charges_.contains((Charge*)param)) printf("  Parameter %4i\t    charge %12.5e (on atom %s)\n", ++n, param->value(), param->name());
		printf("  Parameter %4i\t%10s %12.5e (in %s interaction between atoms %s)\n", ++n, param->name(), param->value(), param->parentForm(), param->parentAtoms());
		if (param->hasLimits()) printf("\t\t--> limits are { %s, %s }\n", param->lowerLimitAsText(), param->upperLimitAsText());
	}
}

// Print dependent parameters
void Forcefield::printDependentParameters()
{
	Parameter *param;
	int n = 0;
	for (param = dependentParameters_.first(); param != NULL; param = param->next)
	{
// 		if (charges_.contains((Charge*)param)) printf("  Parameter %4i\t    charge %12.5e (on atom %s)\n", ++n, param->value(), param->name());
		printf("  Parameter %4i\t%10s %12.5e (in %s interaction between atoms %s)\n", ++n, param->name(), param->value(), param->parentForm(), param->parentAtoms());
		if (param->hasLimits()) printf("\t\t--> limits are { %s, %s }\n", param->lowerLimitAsText(), param->upperLimitAsText());
	}
}

// Return cost penalty for free parameters going outside their defined limits
double Forcefield::freeParameterPenalty()
{
	double penalty = 0, factor;
	for (Parameter *param = freeParameters_.first(); param != NULL; param = param->next)
	{
		if (param->withinLimits()) continue;
		// Bad parameter!  Get its contribution to the penalty function
		factor = param->limitPenalty();
// 		printf("Parameter %s has value %f and limits {%s,%s}, and incurs penalty of %f\n", param->name(), param->value(), param->lowerLimitAsText(), param->upperLimitAsText(), factor*factor);
		penalty += factor*factor;
	}
	return penalty;
}

// Calculate parameter logpoints
void Forcefield::calculateLogPoints()
{
// 	for (int type=0; type<Configuration::nContributionTypes; ++type) logPoints_[type] = 0;
	// Set 'baseline' logpoint values, used to force calculation of forces even if no parameters of a given type are free
	logPoints_[Configuration::BondContribution] = bonds_.nItems();
	logPoints_[Configuration::AngleContribution] = angles_.nItems();
	logPoints_[Configuration::TorsionContribution] = torsions_.nItems();
	logPoints_[Configuration::ChargeContribution] = charges_.nItems();
	logPoints_[Configuration::VdwContribution] = vdw_.nItems();
	// Loop over parameters, adding individual logpoints to relevant totals
	Parameter *param;
	Interaction *ixn;
	for (param = freeParameters_.first(); param != NULL; param = param->next)
	{
		ixn = param->parent();
		if (ixn->type() == Interaction::BondType) logPoints_[Configuration::BondContribution] += param->changeLog();
		else if (ixn->type() == Interaction::AngleType) logPoints_[Configuration::AngleContribution] += param->changeLog();
		else if (ixn->type() == Interaction::TorsionType) logPoints_[Configuration::TorsionContribution] += param->changeLog();
		else if (ixn->type() == Interaction::ChargeType) logPoints_[Configuration::ChargeContribution] += param->changeLog();
		else logPoints_[Configuration::VdwContribution] += param->changeLog();
	}
}

// Return specified logpoint
unsigned int Forcefield::logPoint(Configuration::ContributionType type)
{
	return logPoints_[type];
}

/*
// Methods
*/

// Setup all contained forcefield terms, check for missing parameters etc.
bool Forcefield::setup()
{
	printf("\nFinalising/checking forcefield parameters:\n");
	Interaction *ixn;
	for (ixn = vdw_.first(); ixn != NULL; ixn = ixn->next) if (!ixn->setup()) return FALSE;
	for (ixn = charges_.first(); ixn != NULL; ixn = ixn->next) if (!ixn->setup()) return FALSE;
	for (ixn = bonds_.first(); ixn != NULL; ixn = ixn->next) if (!ixn->setup()) return FALSE;
	for (ixn = angles_.first(); ixn != NULL; ixn = ixn->next) if (!ixn->setup()) return FALSE;
	for (ixn = torsions_.first(); ixn != NULL; ixn = ixn->next) if (!ixn->setup()) return FALSE;
	return TRUE;
}

// Check charge setup
bool Forcefield::checkChargeSetup()
{
	// Count the number of charges to be fitted to see if we need a dependent charge
	int nFreeCharges = 0;
	for (Interaction *ixn = charges_.first(); ixn != NULL; ixn = ixn->next) if (freeParameters_.contains(ixn->parameters()[0])) ++nFreeCharges;
	if (nFreeCharges == 0) printf("\nNo variable charges in forcefield.\n");
	else
	{
		printf("\nThere are %i variable charges in the forcefield.\n", nFreeCharges);
		if (dependentCharge_ == NULL)
		{
			printf("Error: Variable charges have been specified, but no dependent charge has been set.\n");
			return FALSE;
		}
	}

	return TRUE;
}

// Calculate dependent parameters
bool Forcefield::calculateDependencies()
{
	// Dependent charge
	if (dependentCharge_ != NULL)
	{
		// Calculate total charge in model
		double totalCharge = 0.0;
		dependentCharge_->setValue(0.0);
		int depCount = 0;
		Parameter *param;
		for (Interaction *ixn = charges_.first(); ixn != NULL; ixn = ixn->next)
		{
			param = ixn->parameter(0);
			if (param == dependentCharge_) depCount = ixn->usageCount();
			else totalCharge += ixn->usageCount() * param->value();
		}
		dependentCharge_->setValue(-totalCharge/depCount);
// 		printf("Total charge = %f, dependentCharge = %f per %i atoms\n", totalCharge, dependentCharge_->value(), depCount);
	}
	
	// Combined van der Waals parameters
	for (Interaction *ixn = vdw_.first(); ixn != NULL; ixn = ixn->next)
	{
		if (ixn->type() != Interaction::CombinedVdwType) continue;
		if (!ixn->combine())
		{
			printf("Error: Failed to calculate combined vdW parameters for interaction %s.\n", ixn->info());
			return FALSE;
		}
	}
	return TRUE;
}
