/*
	*** Interaction Data and Calculation
	*** interaction.cpp
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

#include "interaction.h"
#include "configuration.h"
#include "atom.h"
#include "atompair.h"
#include <stdio.h>
#include <string.h>

// Interaction Data 
InteractionData interactions[Interaction::nInteractionForms] = {
	{ "_none_", 1,		{ "_dummy_" } },
	{ "charge", 1,		{ "q" } },
	{ "harmonic", 2,	{ "k", "eq" } },
	{ "morse", 3,		{ "e0", "beta", "eq" } },
	{ "cos", 4,		{ "k", "eq", "n", "s" } },
	{ "cos2", 4,		{ "k", "c0", "c1", "c2" } },
	{ "cos3", 3,		{ "c1", "c2", "c3" } },
	{ "harmoniccos", 2,	{ "k", "eq" } },
	{ "lj", 2,		{ "epsilon", "sigma" } },
	{ "eam", 7,		{ "r", "E", "A", "alpha", "beta", "rho", "Z" } }	   // Baskes et al.
// 	{ "eam", 13,		{ "Ec", "phi0", "r0", "alpha", "beta", "gamma", "delta", "c0", "c1", "c2", "c3", "c4", "c5" } }    // Mei et al.
};
Interaction::InteractionForm Interaction::interactionForm(const char* s)
{
	for (int iform=0; iform < Interaction::nInteractionForms; ++iform) if (strcmp(s, interactions[iform].name) == 0) return (Interaction::InteractionForm) iform;
	return Interaction::nInteractionForms;
}

// Constructor
Interaction::Interaction()
{
	// Private variables
	type_ = Interaction::NoType;
	form_ = Interaction::NoForm;
	values_ = NULL;
	sourceTypeA_ = NULL;
	sourceTypeB_ = NULL;
	usageCount_ = 0;
	id_ = -1;

	// Public variables
	prev = NULL;
	next = NULL;
}

// Destructor
Interaction::~Interaction()
{
	if (values_ != NULL) delete[] values_;
}

/*
// Atom Data
*/

// Return whether strings match, considering wildcard characters
bool Interaction::matchWild(const char *s1, const char *s2)
{
	if (strlen(s1) != strlen(s2)) return FALSE;
	for (int n=0; n<strlen(s1); ++n)
	{
		if (s1[n] == s2[n]) continue;
		if (s1[n] == '*' || s2[n] == '*') continue;
		return FALSE;
	}
	return TRUE;
}

// Set atom data for interaction
void Interaction::setAtoms(const char *a1, const char *a2, const char *a3, const char *a4)
{
	types_[0] = a1;
	types_[1] = a2;
	types_[2] = a3;
	types_[3] = a4;
}

// Check stored atom names against supplied names (parameter matching)
bool Interaction::matches(const char *a1, const char *a2, const char *a3, const char *a4)
{
	bool result = FALSE;
	if (a2 == NULL)
	{
		// Match single atom (charge) interaction
		if (types_[0] == a1) result = TRUE;
	}
	else if (a3 == NULL)
	{
		// Match bond interaction i-j or j-i
		if (types_[0] == a1 && types_[1] == a2) result = TRUE;
		else if (types_[1] == a1 && types_[0] == a2) result = TRUE;
	}
	else if (a4 == NULL)
	{
		// Match angle interction i-j-k or k-j-i
		if (types_[1] != a2) result = FALSE;
		else if (types_[0] == a1 && types_[2] == a3) result = TRUE;
		else if (types_[2] == a1 && types_[0] == a3) result = TRUE;
	}
	else
	{
		// Match torsion interaction i-j-k-l or l-k-j-i
		if (types_[0] == a1 && types_[1] == a2 && types_[2] == a3 && types_[3] == a4) result = TRUE;
		else if (types_[3] == a1 && types_[2] == a2 && types_[1] == a3 && types_[0] == a4) result = TRUE;
	}
	return result;
}

// Check stored atoms against supplied names (parameter matching)
bool Interaction::matches(Atom* i, Atom* j, Atom* k, Atom* l)
{
	if (k == NULL) return matches(i->name(), j->name());
	else if (l == NULL) return matches(i->name(), j->name(), k->name());
	else return matches(i->name(), j->name(), k->name(), l->name());
}

// Check stored atoms against supplied atoms (literal name matching)
bool Interaction::matchesWild(const char *a1, const char *a2, const char *a3, const char *a4)
{
	bool result = FALSE;
	if (a3 == NULL)
	{
		// Match bond interaction i-j or j-i
		if (matchWild(types_[0],a1) && matchWild(types_[1],a2)) result = TRUE;
		else if (matchWild(types_[1],a1) && matchWild(types_[0],a2)) result = TRUE;
	}
	else if (a4 == NULL)
	{
		// Match angle interction i-j-k or k-j-i
		if (!matchWild(types_[1],a2)) result = FALSE;
		else if (matchWild(types_[0],a1) && matchWild(types_[2],a3)) result = TRUE;
		else if (matchWild(types_[2],a1) && matchWild(types_[0],a3)) result = TRUE;
	}
	else
	{
		// Match torsion interaction i-j-k-l or l-k-j-i
		if (matchWild(types_[0],a1) && matchWild(types_[1],a2) && matchWild(types_[2],a3) && matchWild(types_[3],a4)) result = TRUE;
		else if (matchWild(types_[3],a1) && matchWild(types_[2],a2) && matchWild(types_[1],a3) && matchWild(types_[0],a4)) result = TRUE;
	}
	return result;
}

// Check stored atoms against supplied atoms (wildcard matching)
bool Interaction::matchesWild(Atom *i, Atom *j, Atom *k, Atom *l)
{
	if (k == NULL) return matchesWild(i->name(), j->name());
	else if (l == NULL) return matchesWild(i->name(), j->name(), k->name());
	else return matchesWild(i->name(), j->name(), k->name(), l->name());
}

// Return type name specified
const char *Interaction::type(int n)
{
	return types_[n];
}

// Return formatted atom interaction list text
const char *Interaction::atomsText()
{
	static Dnchar text;
	text = types_[0];
	for (int n=1; n<4; ++n) if (!types_[n].isEmpty()) text.strcatf("-%s",types_[n].get());
	return text.get();
}

/*
// Form / Parameter Data
*/

// Set ID of interaction
void Interaction::setId(int id)
{
	id_ = id;
}

// Return ID of interaction
int Interaction::id()
{
	return id_;
}

// Set type of interaction
void Interaction::setType(InteractionType type)
{
	type_ = type;
}

// Return type of interaction
Interaction::InteractionType Interaction::type()
{
	return type_;
}

// Add parameter to list
bool Interaction::addParameter(InteractionForm form, Parameter *param)
{
	// Check type of interaction
	if (type_ == Interaction::NoType)
	{
		printf("Error: Type of interaction must be set before a parameter can be added.\n");
		return FALSE;
	}
	
	// Check functional form: First, if one hasn't been set (form_ == 0). Second, if the same form is being set (form_ == form)
	if (form_ == Interaction::NoForm)
	{
		
		if ((form >= 0) && (form <= Interaction::nInteractionForms)) form_ = form;
		else 
		{
			if (type_ == Interaction::VdwType) printf("Error: Functional form ID %i is out of range for van der Waals interaction.\n", form);
			else if (type_ == Interaction::BondType) printf("Error: Functional form ID %i is out of range for bond interaction.\n", form);
			else if (type_ == Interaction::AngleType) printf("Error: Functional form ID %i is out of range for angle interaction.\n", form);
			else if (type_ == Interaction::TorsionType) printf("Error: Functional form ID %i is out of range for torsion interaction.\n", form);
			return FALSE;
		}
	}
	else if (form_ != form)
	{
		if (type_ == Interaction::VdwType) printf("Error: Functional form of van der Waals interaction %s-%s has changed between parameter definitions (was '%s', now '%s').\n", types_[0].get(), types_[1].get(), interactions[form_].name, interactions[form].name);
		else if (type_ == Interaction::BondType) printf("Error: Functional form of bond interaction %s-%s has changed between parameter definitions (was '%s', now '%s').\n", types_[0].get(), types_[1].get(), interactions[form_].name, interactions[form].name); 
		else if (type_ == Interaction::AngleType) printf("Error: Functional form of angle interaction %s-%s-%s has changed between parameter definitions (was '%s', now '%s').\n", types_[0].get(), types_[1].get(), types_[2].get(), interactions[form_].name, interactions[form].name); 
		else if (type_ == Interaction::TorsionType) printf("Error: Functional form of torsion interaction %s-%s-%s-%s has changed between parameter definitions (was '%s', now '%s').\n", types_[0].get(), types_[1].get(), types_[2].get(), types_[3].get(), interactions[form_].name, interactions[form].name);
		return FALSE;
	}

	// The interaction form is valid, so search the form's list of parameters to see if this one exists
	// XXX TODO Spline fitting will need to bypass by some means...
	int n;
	for (n=0; n<interactions[form_].nParameters; ++n) if (strcmp(param->name(), interactions[form_].parameters[n]) == 0) break;
	if (n == interactions[form_].nParameters)
	{
		printf("Error: Parameter '%s' is not present in an interaction of functional form '%s'.\n", param->name(), interactions[form_].name);
		printf("Valid parameter names are:\n\t");
		for (n=0; n<interactions[form_].nParameters; ++n) printf("%s  ", interactions[form_].parameters[n]);
		printf("\n");
		return FALSE;
	}
	// Has this parameter already been defined?
	if (parameters_.containsData(n))
	{
		if (type_ == Interaction::VdwType) printf("Error: van der Waals interaction %s-%s already has parameter '%s' set.\n", types_[0].get(), types_[1].get(), interactions[form_].parameters[n]);
		else if (type_ == Interaction::BondType) printf("Error: bond interaction %s-%s already has parameter '%s' set.\n", types_[0].get(), types_[1].get(), interactions[form_].parameters[n]); 
		else if (type_ == Interaction::AngleType) printf("Error: angle interaction %s-%s-%s already has parameter '%s' set.\n", types_[0].get(), types_[1].get(), types_[2].get(), interactions[form_].parameters[n]); 
		else if (type_ == Interaction::TorsionType) printf("Error: torsion interaction %s-%s-%s-%s already has parameter '%s' set.\n", types_[0].get(), types_[1].get(), types_[2].get(), types_[3].get(), interactions[form_].parameters[n]);
		return FALSE;
	}
	// Add parameter to our reference list, storing parameter index as well
	parameters_.add(param, n);
	param->setParent(this);
	printf("Added parameter '%s' to interaction.\n", param->name());
	return TRUE;
}

// Add all parameters relevant to the specified functional form to the list
bool Interaction::addAllParameters(InteractionForm form, List<Parameter> &destinationList, Parameter::ParameterType type, Parameter::ParameterStrength strength)
{
	if (form_ != Interaction::NoForm)
	{
		printf("Error: Refusing to add all relevant parameters to an interaction which has already been (partially) set.\n");
		return FALSE;
	}
	for (int n=0; n<interactions[form].nParameters; ++n)
	{
		Parameter *par = destinationList.add();
		par->set(interactions[form].parameters[n], 0.0);
		par->setType(type, strength);
		if (!addParameter(form, par))
		{
			printf("Error adding parameter to combined interaction.\n");
			return FALSE;
		}
	}
	return TRUE;
}

// Return functional form of interaction
Interaction::InteractionForm Interaction::form()
{
	return form_;
}

// Return functional form text of interaction
const char *Interaction::formText()
{
	return interactions[form_].name;
}

// Check validity of defined parameters and finalise function-specific parameter list
bool Interaction::setup()
{
	// Check type and form
	if (type_ == Interaction::NoType)
	{
		printf("Error: Type of interaction has not been set.\n");
		return FALSE;
	}
	if (form_ == Interaction::NoForm)
	{
		if (type_ == Interaction::VdwType) printf("Error: Functional form of van der Waals interaction %s-%s has not been set.\n", types_[0].get(), types_[1].get());
		else if (type_ == Interaction::BondType) printf("Error: Functional form of bond interaction %s-%s has not been set.\n", types_[0].get(), types_[1].get());
		else if (type_ == Interaction::AngleType) printf("Error: Functional form of angle interaction %s-%s-%s has not been set.\n", types_[0].get(), types_[1].get(), types_[2].get());
		else if (type_ == Interaction::TorsionType) printf("Error: Functional form of torsion interaction %s-%s-%s-%s has not been set.\n", types_[0].get(), types_[1].get(), types_[2].get(), types_[3].get());
		return FALSE;
	}
	
	// Run over required parameters and build up list of pointers to their values
	int n;
	values_ = new Parameter*[interactions[form_].nParameters];
	for (n=0; n<interactions[form_].nParameters; ++n) values_[n] = NULL;
	n = 0;
	for (Refitem<Parameter,int> *ri = parameters_.first(); ri != NULL; ri = ri->next)
	{
		// Double-check range of data value and corresponding parameter name
		if (ri->data < 0 || ri->data >= interactions[form_].nParameters)
		{
			printf("Error: Stored parameter index is out of range for interaction.\n");
			return FALSE;
		}
		if (strcmp(interactions[form_].parameters[ri->data],ri->item->name()) != 0)
		{
			printf("Error: Stored parameter name does not match corresponding entry in interaction data array.\n");
			return FALSE;
		}
		if (values_[ri->data] != NULL)
		{
			printf("Error: Parameter index has already been stored.\n");
			return FALSE;
		}
		// All is okay, so store pointer and increase counter
		values_[ri->data] = ri->item;
		++n;
	}

	// Check final stored pointer count
	if (n != interactions[form_].nParameters)
	{
		printf("Error: Not all parameters for interaction have been provided (interaction type is '%s', atoms are %s-%s-%s-%s).\n", interactions[form_].name, types_[0].get(), types_[1].get(), types_[2].get(), types_[3].get());
		return FALSE;
	}
	return TRUE;
}

// Return string containing formatted form/parameter data
const char *Interaction::info()
{
	static Dnchar text;
	text.sprintf("%s,  ", interactions[form_].name);
	for (int n=0; n<interactions[form_].nParameters; ++n) text.strcatf("%s=%f  ", values_[n]->name(), values_[n]->value());
	return text.get();
}

// Set combination source interactions (only if type_ == CombinedVdwType)
void Interaction::setCombinationSource(Interaction *a, Interaction *b)
{
	sourceTypeA_ = a;
	sourceTypeB_ = b;
}

// Perform parameter combination
bool Interaction::combine()
{
	if (sourceTypeA_ == NULL || sourceTypeB_ == NULL)
	{
		printf("BUG: One or both combination sources are unset in Interaction.\n");
		return FALSE;
	}
	if (type_ != Interaction::CombinedVdwType)
	{
		printf("BUG: This interaction is not of CombinedVdwType.\n");
		return FALSE;
	}

	// Combine based on functional form
	switch (form_)
	{
		case (Interaction::LennardJones126Form):
			// Parameter order: Epsilon, Sigma
			values_[0]->setValue(sqrt(sourceTypeA_->values_[0]->value() * sourceTypeB_->values_[0]->value()));
			values_[1]->setValue(sqrt(sourceTypeA_->values_[1]->value() * sourceTypeB_->values_[1]->value()));
			break;
		default:
			printf("Error: No code to combine functional form of type '%s'. Complain to programmer.\n", interactions[form_].name);
			return FALSE;
	}
	return TRUE;
}

// Return whether all parameters involved in the interaction are variable
bool Interaction::allParametersVariable()
{
	for (int n=0; n<interactions[form_].nParameters; ++n) if (values_[n]->type() != Parameter::FreeParameter) return FALSE;
	return TRUE;
}

// Return number of parameters associated with interaction
int Interaction::nParameters()
{
	if (form_ == Interaction::NoForm) printf("BUG: Form of interaction type has not been set, so nParameters cannot be returned.\n");
	else return interactions[form_].nParameters;
	return 0;
}

// Return list of parameters associated to interaction
Parameter **Interaction::parameters()
{
	return values_;
}

// Return nth parameter in list
Parameter *Interaction::parameter(int n)
{
	if ((n < 0) || (n >= interactions[form_].nParameters))
	{
		printf("Internal Error: Parameter index %i is out of range for interaction.\n", n);
		return NULL;
	}
	return values_[n];
}

// Increase usage count
void Interaction::increaseUsageCount(int delta)
{
	usageCount_ += delta;
}

// Return charge count
int Interaction::usageCount()
{
	return usageCount_;
}

/*
// Energy / Force Calculation
*/

// Calculate energy and forces for provided atom pair
bool Interaction::calculate(Configuration* cfg, AtomPair *pair)
{
	// Call relevant subroutine to calculate energy
	bool result = FALSE;
	if (type_ == Interaction::BondType || type_ == Interaction::VdwType || type_ == Interaction::CombinedVdwType)
	{
		double u, du_dx;
		Vec3<double> fi;
		
		// Distance between atoms (in Angstroms) is the variable x
		result = evaluate(pair->rij(), u, du_dx);
		u *= pair->vdwScale();
		du_dx *= pair->vdwScale();
		
		// Calculate forces and sum, along with energy
		fi = (pair->vij() / pair->rij()) * -du_dx;
		if (type_ == Interaction::BondType)
		{
			cfg->subtractForce(Configuration::BondContribution, fi, pair->i());
			cfg->addForce(Configuration::BondContribution, fi, pair->j());
			cfg->addEnergy(Configuration::BondContribution, u);
		}
		else
		{
			cfg->subtractForce(Configuration::VdwContribution, fi, pair->i());
			cfg->addForce(Configuration::VdwContribution, fi, pair->j());
			cfg->addEnergy(Configuration::VdwContribution, u);
		}
	}
	else
	{
		printf("BUG: Interaction type cannot be calculated from an AtomPair.\n");
		return FALSE;
	}
	return TRUE;
}

// Calculate energy and forces between specified atoms
bool Interaction::calculate(Configuration* cfg, int i, int j, int k, int l, double scale)
{
	// Call relevant subroutine to calculate energy
	bool result = FALSE;
	if (type_ == Interaction::VdwType || type_ == Interaction::CombinedVdwType)
	{
		double rij, u, du_dx;
		Vec3<double> vij, fi;
		
		// Distance between atoms (in Angstroms) is the variable x
		vij = cfg->mimd(i,j);
		rij = vij.magnitude();
		result = evaluate(rij, u, du_dx);
		u *= scale;
		du_dx *= scale;
		
		// Calculate forces and sum, along with energy
		fi = (vij / rij) * -du_dx;
		cfg->subtractForce(Configuration::VdwContribution, fi, i);
		cfg->addForce(Configuration::VdwContribution, fi, j);
		cfg->addEnergy(Configuration::VdwContribution, u);
	}
	else if (type_ == Interaction::BondType)
	{
		double rij, u, du_dx;
		Vec3<double> vij, fi;
		
		// Distance between atoms (in Angstroms) is the variable x
		vij = cfg->mimd(i,j);
		rij = vij.magnitude();
		result = evaluate(rij, u, du_dx);
		u *= scale;
		du_dx *= scale;
		
		// Calculate forces and sum, along with energy
		fi = (vij / rij) * -du_dx;
		cfg->subtractForce(Configuration::BondContribution, fi, i);
		cfg->addForce(Configuration::BondContribution, fi, j);
		cfg->addEnergy(Configuration::BondContribution, u);
	}
	else if (type_ == Interaction::AngleType)
	{
		double rij, rkj, u, du_dx, dp, angle;
		Vec3<double> vji, vjk, fi, fk;
		
		// Angle between atoms (in radians) is the variable x
		// Calculate minimum image vectors w.r.t. atom j
		vji = cfg->mimd(j,i);
		vjk = cfg->mimd(j,k);
		// Normalise vectors, calculate dot product and angle.
		rij = vji.magAndNormalise();
		rkj = vjk.magAndNormalise();
		dp = vji.dp(vjk);
		if (dp < -1.0) dp = -1.0;
		else if (dp > 1.0) dp = 1.0;
		angle = acos(dp);

		result = evaluate(angle, u, du_dx);
		u *= scale;
		du_dx *= scale;

		// Calculate atomic forces and sum, along with energy
		du_dx *= -1.0 / sin(angle);
		fi = vjk - vji * dp;
		fi *= -du_dx / rij;
		fk = vji - vjk * dp;
		fk *= -du_dx / rkj;
		cfg->addForce(Configuration::AngleContribution, fi, i);
		cfg->subtractForce(Configuration::AngleContribution, fi, j);
		cfg->subtractForce(Configuration::AngleContribution, fk, j);
		cfg->addForce(Configuration::AngleContribution, fk, k);
		cfg->addEnergy(Configuration::AngleContribution, u);
	}
	else if (type_ == Interaction::TorsionType)
	{
		Vec3<double> vji, vjk, vkl, xpj, xpk, dcos_dxpj, dcos_dxpk, fi, fj, fk, fl;
		double rji, rjk, rkl, mag_xpj, mag_xpk, dp, angle, u, du_dx;
		static Matrix dxpj_dij, dxpj_dkj, dxpk_dkj, dxpk_dlk;
		
		// Calculate vectors between atoms
		vji = cfg->mimd(j,i);
		vjk = cfg->mimd(j,k);
		vkl = cfg->mimd(k,l);
		rji = vji.magnitude();
		rjk = vjk.magnitude();
		rkl = vkl.magnitude();
		
// 		printf("i-j-k-l %i-%i-%i-%i MAGs %f %f %f\n",i,j,k,l, mag_ij,mag_kj, mag_lk);
		// Calculate cross products and torsion angle formed (in radians)
		xpj = vji * vjk;
		xpk = vkl * vjk;
		mag_xpj = xpj.magAndNormalise();
		mag_xpk = xpk.magAndNormalise();
		dp = xpj.dp(xpk);
		if (dp < -1.0) dp = -1.0;
		else if (dp > 1.0) dp = 1.0;
		angle = acos(dp);
		
		result = evaluate(angle, u, du_dx);
		u *= scale;
		du_dx *= scale;

		// Calculate atomic forces and sum, along with energy
		du_dx *= -1.0 / sin(angle);
		// Construct matrices for force computation
		dxpj_dij.cyclicPermute(vjk);
		dxpj_dkj.cyclicPermute(-vji);
		dxpk_dkj.cyclicPermute(-vkl);
		dxpk_dlk.cyclicPermute(vjk);
		// Construct derivatives of cos(phi) w.r.t. perpendicular axes
		dcos_dxpj = (xpk - xpj * dp) / mag_xpj;
		dcos_dxpk = (xpj - xpk * dp) / mag_xpk;
		fi.x = du_dx * dcos_dxpj.dp(dxpj_dij.columnAsVec3(0));
		fi.y = du_dx * dcos_dxpj.dp(dxpj_dij.columnAsVec3(1));
		fi.z = du_dx * dcos_dxpj.dp(dxpj_dij.columnAsVec3(2));

		fj.x = du_dx * ( dcos_dxpj.dp( -dxpj_dij.columnAsVec3(0) - dxpj_dkj.columnAsVec3(0) ) - dcos_dxpk.dp(dxpk_dkj.columnAsVec3(0)) );
		fj.y = du_dx * ( dcos_dxpj.dp( -dxpj_dij.columnAsVec3(1) - dxpj_dkj.columnAsVec3(1) ) - dcos_dxpk.dp(dxpk_dkj.columnAsVec3(1)) );
		fj.z = du_dx * ( dcos_dxpj.dp( -dxpj_dij.columnAsVec3(2) - dxpj_dkj.columnAsVec3(2) ) - dcos_dxpk.dp(dxpk_dkj.columnAsVec3(2)) );

		fk.x = du_dx * ( dcos_dxpk.dp( dxpk_dkj.columnAsVec3(0) - dxpk_dlk.columnAsVec3(0) ) + dcos_dxpj.dp(dxpj_dkj.columnAsVec3(0)) );
		fk.y = du_dx * ( dcos_dxpk.dp( dxpk_dkj.columnAsVec3(1) - dxpk_dlk.columnAsVec3(1) ) + dcos_dxpj.dp(dxpj_dkj.columnAsVec3(1)) );
		fk.z = du_dx * ( dcos_dxpk.dp( dxpk_dkj.columnAsVec3(2) - dxpk_dlk.columnAsVec3(2) ) + dcos_dxpj.dp(dxpj_dkj.columnAsVec3(2)) );

		fl.x = du_dx * dcos_dxpk.dp(dxpk_dlk.columnAsVec3(0));
		fl.y = du_dx * dcos_dxpk.dp(dxpk_dlk.columnAsVec3(1));
		fl.z = du_dx * dcos_dxpk.dp(dxpk_dlk.columnAsVec3(2));

		cfg->addForce(Configuration::TorsionContribution, fi, i);
		cfg->addForce(Configuration::TorsionContribution, fj, j);
		cfg->addForce(Configuration::TorsionContribution, fk, k);
		cfg->addForce(Configuration::TorsionContribution, fl, l);
		
		cfg->addEnergy(Configuration::TorsionContribution, u);
// 			printf("i-j-k-l %i-%i-%i-%i DP %16.13f %16.13f %16.13f %16.13f\n",i,j,k,l, mag_xpj, mag_xpk, dp, phi);
// 		if (phi < 0.0) { if (phi > -1e-8) phi = -1e-8; }
// 		else if (phi < 1e-8) phi = 1e-8;
		// Derivative w.r.t. change in torsion angle
// 		dphi_dcosphi = -1.0 / sin(phi);
		// Pathological case where phi = 0...
// 		if (phi < 1E-8) dphi_dcosphi = 0.0;
// 		printf("Error: Torsion interactions are not yet implemented.\n");
// 		result = FALSE;
	}
	return result;
}

// Evaluate value and derivative of interaction at supplied value
bool Interaction::evaluate(double x, double &fx, double &df_dx)
{
	double forcek, eq, e0, expo, beta, c0, c1, c2, c3, delta, n, s, epsilon, sigmax, sigmax6, sigmax12;
	
	// Determine functional form
	switch (form_)
	{
		case (Interaction::HarmonicForm):
			forcek = fabs(values_[0]->value());
			eq = values_[1]->value();
			// F(x) = 0.5 * forcek * (r - eq)**2
			x -= eq;
			fx = 0.5 * forcek * x * x;
			// F'(x) = forcek * (x - eq)
			df_dx = forcek * x;
			break;
		case (Interaction::MorseForm):
			e0 = values_[0]->value();
			beta = values_[1]->value();
			eq = values_[2]->value();
			// F(x) = E0 * (1 - exp( -B(x - r0) ) )**2
			x -= eq;
			expo = exp( -beta * x );
			fx = e0 * ( (1.0-expo)*(1.0-expo) );
			// F'(x) = 2 * B * E0 * (1 - exp( -B(x - r0) ) ) * exp( -B*(x - r0) )
			df_dx = 2.0 * beta * e0 * (1.0 - expo) * expo;
			break;
		case (Interaction::CosineForm):
			forcek = values_[0]->value();
			eq = values_[1]->value() / DEGRAD;
			n = values_[2]->value();
			s = values_[3]->value();
			// F(x) = forcek * (1 + s * cos(n*x - eq))
			fx = forcek * (1.0 + s * cos(n * x - eq));
			// F'(x) = -forcek * n * s * sin(n*x - eq)
			df_dx = -forcek * n * s * sin(n * x - eq);
			break;
		case (Interaction::DoubleCosineForm):
			forcek = values_[0]->value();
			c0 = values_[1]->value();
			c1 = values_[2]->value();
			c2 = values_[3]->value();
			// F(x) = forcek * (C0 + C1 * cos(x) + C2 * cos(2*x))
			fx = forcek * (c0 + c1 * cos(x) + c2 * cos(2.0 * x));
			// F'(x) = -forcek * (c1 * sin(x) + 2 * c2 * sin(2*x))
			df_dx = -forcek * (c1 * sin(x) + 2.0 * c2 * sin(2.0 * x));
			break;
		case (Interaction::TripleCosineForm):
			c1 = values_[0]->value();
			c2 = values_[1]->value();
			c3 = values_[2]->value();
			// F(x) = 0.5 * ( c1*(1+cos(x)) + c2*(1-cos(2*x)) + c3*(1+cos(3*x)) )
			fx = 0.5 * (c1 * (1.0 + cos(x)) + c2 * (1.0 - cos(2.0*x)) + c3 * (1.0 + cos(3.0*x)));
			// F'(x) = 0.5 * ( -c1*sin(x) + 2 * c2*sin(2*x) - 3 * c3*(sin(3*x)) )
			df_dx = 0.5 * ( -c1*sin(x) + 2.0*c2*sin(2.0*x) - 3.0*c3*sin(3.0*x));
			break;
		case (Interaction::HarmonicCosineForm):
			forcek = values_[0]->value();
			eq = values_[1]->value() / DEGRAD;
			delta = cos(x) - cos(eq);
			// F(x) = 0.5 * forcek * (cos(x) - cos(eq)))**2
			fx = 0.5 * forcek * delta * delta;
			// F'(x) = forcek * (cos(x) - cos(eq))) * -sin(x)
			df_dx = -forcek * delta * sin(x);
			break;
		case (Interaction::LennardJones126Form):
			epsilon = values_[0]->value();
			sigmax = values_[1]->value() / x;
			sigmax6 = sigmax*sigmax*sigmax;
			sigmax6 *= sigmax6;
			sigmax12 = sigmax6*sigmax6;
			// F(x) = 4 * epsilon * ( (sigma/r)**12 - (sigma/r)**6 )
			fx = 4.0 * epsilon * (sigmax12 - sigmax6);
			// F'(x) = 48 * epsilon * ( sigma**12/r**13 - 0.5*(sigma**6/r**7) )
			df_dx = 48.0 * epsilon * ( -sigmax12/x + 0.5*(sigmax6/x));
			break;
		default:
			printf("No equation coded for '%s' interaction.\n", interactions[form_].name);
			fx = 0.0;
			df_dx = 0.0;
			return FALSE;
	}
	return TRUE;
}
