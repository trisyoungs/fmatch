/*
	*** System Definition - Input/Output Routines
	*** io.cpp
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

#include "system.h"
#include "lineparser.h"
#include <string.h>

// Input Keywords
const char *InputKeywords[TargetSystem::nInputKeywords] = { "angle", "atom", "bond", "charge", "combine", "configurations", "cubic", "cutoff", "econvert", "end", "ewald", "equalise", "ev", "file", "fit", "fixed", "forcemultiplier", "frames", "masses", "nmols", "nopairs", "orthorhombic", "parallelepiped", "penalty", "relative", "replicate", "retain", "schedule", "seed", "species", "strategy", "system", "threshold", "torsion", "vdw" };
TargetSystem::InputKeyword TargetSystem::inputKeyword(const char *s)
{
	for (int ik=0; ik < TargetSystem::nInputKeywords; ++ik) if (strcmp(s, InputKeywords[ik]) == 0) return (TargetSystem::InputKeyword) ik;
	return TargetSystem::nInputKeywords;
}


// FMatch3 Configurations File Keywords
const char *FMatch3Keywords[TargetSystem::nFMatch3Keywords] = { "configuration", "atoms", "atomsandforces", "cell", "energy", "energyperatom", "esp", "forces" };
TargetSystem::FMatch3Keyword TargetSystem::fmatch3Keyword(const char *s)
{
	for (int ik=0; ik < TargetSystem::nFMatch3Keywords; ++ik) if (strcmp(s, FMatch3Keywords[ik]) == 0) return (TargetSystem::FMatch3Keyword) ik;
	return TargetSystem::nFMatch3Keywords;
}

// Get characters before first occurrence of designated character
const char *beforeChar(const char *s, char delim)
{
        static Dnchar result(1024);
        result.clear();
        for (int i = 0; s[i] != '\0'; i++)
        {
                if (s[i] == delim) break;
                result += s[i];
        }
        return result;
}

// Get characters after first occurrence of designated character
const char *afterChar(const char *s, char delim)
{
        static Dnchar result(1024);
        result.clear();
        bool found = FALSE;
        for (int i = 0; s[i] != '\0'; i++)
        {
                if (found) result += s[i];
                if (s[i] == delim) found = TRUE;
        }
        return result;
}

// Add parameters to interaction
bool TargetSystem::addParameters(LineParser &parser, Interaction *ixn, Parameter::ParameterType type, Interaction::InteractionForm form, Parameter::ParameterStrength strength, int startArg)
{
	Dnchar lowerLimit, upperLimit;
	int n = startArg;
	do 
	{
		Parameter *par = forcefield_.addParameter(type, strength);
		// Get parameter name and value
		par->set(parser.argc(n), parser.argd(n+1));
		n += 2;
		if (!ixn->addParameter(form, par))
		{
			printf("Error adding parameter to interaction.\n");
			return FALSE;
		}

		// Were limits given?
		if (parser.hasArg(n) && (strchr(parser.argc(n),',') != NULL))
		{
			// If this is a dependent parameter, then conplain
			if (type == Parameter::DependentParameter)
			{
				printf("Error: Limits can not be specified for dependent parameters.\n");
				return FALSE;
			}

			// Split arg around ','
			if (strchr(parser.argc(n),',') == NULL)
			{
				printf("Error:: Parameter range must be in the format \"[min,max]\" ('min' or 'max' may be omitted to indicate no limit, but the comma must still be present).\n");
				return FALSE;
			}
			
			lowerLimit = beforeChar(parser.argc(n), ',');
			upperLimit = afterChar(parser.argc(n), ',');
			if (!lowerLimit.isEmpty()) par->setLowerLimit(lowerLimit.asDouble());
			if (!upperLimit.isEmpty()) par->setUpperLimit(upperLimit.asDouble());
			// Increase arg counter
			++n;
		}
		
	} while (n<parser.nArgs());
	return TRUE;
}

// Open specification from supplied file
bool TargetSystem::open(const char *filename)
{
	LineParser parser(filename);
	
	if (!parser.isFileGood())
	{
		printf("Error opening input file '%s'.\n", filename);
		return FALSE;
	}
	
	// Loop over lines in file
	TargetSystem::InputKeyword kwd;
	StrategyStep::StepFunction sf;
	StrategyStep *step;
	LoopStep *loopStep;
	Species *sp;
	Configuration *cfg;
	Interaction::InteractionForm iform;
	Interaction *ixn;
	Parameter *par;
	int n, m, cellType;
	bool done, dependent;
	Parameter::ParameterType parameterType;
	Matrix cell;
	double temp, forceMultiplier;
	Vec3<int> rep(0,0,0);
	bool replicate = FALSE;
	
	// Create strategy stack
	Reflist< List<StrategyStep>, bool> strategyStack;

	do
	{
		parser.getArgsDelim();
		// First argument is always keyword
		kwd = inputKeyword(parser.argc(0));
		switch (kwd)
		{
			// Unrecognised keyword
			case (TargetSystem::nInputKeywords):
				printf("Error: Unrecognised keyword '%s' in input file.\n", parser.argc(0));
				return FALSE;
			// Species definition block
			case (TargetSystem::SpeciesKeyword):
				if (parser.nArgs() != 2)
				{
					printf("Error: 'species' keyword expects exactly one argument in thie context:\n\t[%s]\n", parser.line());
					return FALSE;
				}
				// Second argument on line is species name, so create new species under this title
				sp = species_.add();
				if (parser.hasArg(1)) sp->setName(parser.argc(1));
				if (sp->prev != NULL) sp->setFirstAtom(sp->prev->firstAtom() + sp->prev->nAtoms()*sp->prev->nMolecules());
				done = FALSE;
				// Now cycle over keywords, setting data until we find an 'end' command
				do
				{
					parser.getArgsDelim();
					kwd = inputKeyword(parser.argc(0));
					switch (kwd)
					{
						// Add angle to species
						case (TargetSystem::AngleKeyword):
							if (parser.nArgs() != 4)
							{
								printf("Error: 'angle' keyword expects exactly three arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							sp->addAngle(parser.argi(1)-1, parser.argi(2)-1, parser.argi(3)-1);
							break;
						// Add atom(s) to species
						case (TargetSystem::AtomKeyword):
							if (parser.nArgs() < 2)
							{
								printf("Error: 'atom' keyword expects at least one argument in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							for (n = 1; n<parser.nArgs(); ++n) sp->addAtom(parser.argc(n));
							break;
						// Add bond to species
						case (TargetSystem::BondKeyword):
							if (parser.nArgs() != 3)
							{
								printf("Error: 'bond' keyword expects exactly two arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							sp->addBond(parser.argi(1)-1, parser.argi(2)-1);
							break;
						// Add torsion to species
						case (TargetSystem::TorsionKeyword):
							if (parser.nArgs() != 7)
							{
								printf("Error: 'torsion' keyword expects exactly six arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							sp->addTorsion(parser.argi(1)-1, parser.argi(2)-1, parser.argi(3)-1, parser.argi(4)-1, parser.argd(5), parser.argd(6));
							break;
						// Set mass(es) for atom(s) in species
						case (TargetSystem::MassesKeyword):
							if (parser.nArgs() < 2)
							{
								printf("Error: 'mass' keyword expects at least one argument in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							for (n = 1; n<parser.nArgs(); ++n) if (!sp->setMass(n-1, parser.argd(n))) return FALSE;
							break;
						// Number of molecules
						case (TargetSystem::NMolsKeyword):
							sp->setNMolecules(parser.argi(1));
							break;
						// Unrecognised keyword
						case (TargetSystem::nInputKeywords):
							printf("Error: Unrecognised keyword '%s' found in 'species' block.\n", parser.argc(0));
							return FALSE;
						// End
						case (TargetSystem::EndKeyword):
							done = TRUE;
							break;
						// Misplaced keyword
						default:
							printf("Error: Keyword '%s' used out of context.\n", parser.argc(0));
							return FALSE;
					}
					if (parser.eofOrBlank())
					{
						if (!done)
						{
							printf("Error: Unterminated '%s' block.\n", parameterType == Parameter::FixedParameter ? "fix" : "fit");
							return FALSE;
						}
						done = TRUE;
					}
				} while (!done);
				break;
			// Forcefield definitions
			case (TargetSystem::FitKeyword):
			case (TargetSystem::FixedKeyword):
				// Loop over items in block
				parameterType = (kwd == TargetSystem::FixedKeyword ? Parameter::FixedParameter : Parameter::FreeParameter);
				done = FALSE;
				// Now cycle over keywords, setting data until we find an 'end' command
				do
				{
					parser.getArgsDelim();
					kwd = inputKeyword(parser.argc(0));
					switch (kwd)
					{
						// Add angle function parameter(s)
						case (TargetSystem::AngleKeyword):
							if (parser.nArgs() < 7)
							{
								printf("Error: 'angle' keyword expects at least six arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							ixn = forcefield_.findAngle(parser.argc(1), parser.argc(2), parser.argc(3));
							if (ixn == NULL)
							{
								ixn = forcefield_.addAngle();
								ixn->setAtoms(parser.argc(1), parser.argc(2), parser.argc(3));
								ixn->setType(Interaction::AngleType);
							}
							iform = Interaction::interactionForm(parser.argc(4));
							if (iform == Interaction::nInteractionForms)
							{
								printf("Error: '%s' is not a valid interaction functional form.\n", parser.argc(4));
								return FALSE;
							}
							if (!addParameters(parser, ixn, parameterType, iform, Parameter::Strong, 5))
							{
								printf("Error adding parameters to interaction.\n");
								return FALSE;
							}
							break;
						// Add bond function parameter(s)
						case (TargetSystem::BondKeyword):
							if (parser.nArgs() < 6)
							{
								printf("Error: 'bond' keyword expects at least five arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							ixn = forcefield_.findBond(parser.argc(1), parser.argc(2));
							if (ixn == NULL)
							{
								ixn = forcefield_.addBond();
								ixn->setAtoms(parser.argc(1), parser.argc(2));
								ixn->setType(Interaction::BondType);
							}
							iform = Interaction::interactionForm(parser.argc(3));
							if (iform == Interaction::nInteractionForms)
							{
								printf("Error: '%s' is not a valid interaction functional form.\n", parser.argc(3));
								return FALSE;
							}
							if (!addParameters(parser, ixn, parameterType, iform, Parameter::Strong, 4))
							{
								printf("Error adding parameters to interaction.\n");
								return FALSE;
							}
							break;
						// Add torsion function parameter(s)
						case (TargetSystem::TorsionKeyword):
							if (parser.nArgs() < 8)
							{
								printf("Error: 'torsion' keyword expects at least seven arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							ixn = forcefield_.findTorsion(parser.argc(1), parser.argc(2), parser.argc(3), parser.argc(4));
							if (ixn == NULL)
							{
								ixn = forcefield_.addTorsion();
								ixn->setAtoms(parser.argc(1), parser.argc(2), parser.argc(3), parser.argc(4));
								ixn->setType(Interaction::TorsionType);
							}
							iform = Interaction::interactionForm(parser.argc(5));
							if (iform == Interaction::nInteractionForms)
							{
								printf("Error: '%s' is not a valid interaction functional form.\n", parser.argc(5));
								return FALSE;
							}
							if (!addParameters(parser, ixn, parameterType, iform, Parameter::Moderate, 6))
							{
								printf("Error adding parameters to interaction.\n");
								return FALSE;
							}
							break;
						// Van der waals definition
						case (TargetSystem::VdwKeyword):
							if (parser.nArgs() < 6)
							{
								printf("Error: 'vdw' keyword expects at least five arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							ixn = forcefield_.findVdw(parser.argc(1), parser.argc(2));
							if (ixn == NULL)
							{
								ixn = forcefield_.addVdw();
								ixn->setAtoms(parser.argc(1), parser.argc(2));
								ixn->setType(Interaction::VdwType);
							}
							iform = Interaction::interactionForm(parser.argc(3));
							if (iform == Interaction::nInteractionForms)
							{
								printf("Error: '%s' is not a valid interaction functional form.\n", parser.argc(3));
								return FALSE;
							}
							// For use of EAM potential, only a single atom type is allowed
							if (iform == Interaction::EmbeddedAtomForm)
							{
								if (forcefield_.nVdw() > 1)
								{
									printf("Error: When using EAM potential, only a single atom type is allowed.\n");
									return FALSE;
								}
								else useEam_ = TRUE;
							}
							if (!addParameters(parser, ixn, parameterType, iform, Parameter::Weak, 4))
							{
								printf("Error adding parameters to interaction.\n");
								return FALSE;
							}
							break;
						// Combined Van der waals definition
						case (TargetSystem::CombineKeyword):
							if (parser.nArgs() != 3)
							{
								printf("Error: 'combine' keyword expects exactly two arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							if (strcmp(parser.argc(1), parser.argc(2)) == 0)
							{
								printf("Error: 'combine' cannot specify a self-interaction.\n");
								return FALSE;
							}
							ixn = forcefield_.findVdw(parser.argc(1), parser.argc(2));
							if (ixn == NULL)
							{
								ixn = forcefield_.addVdw();
								ixn->setAtoms(parser.argc(1), parser.argc(2));
								ixn->setType(Interaction::CombinedVdwType);
							}
							break;
						// Charge specification
						case (TargetSystem::ChargeKeyword):
							if (parser.nArgs() < 3)
							{
								printf("Error: 'charge' keyword expects at least two arguments in this context:\n\t[%s]\n", parser.line());
								return FALSE;
							}
							ixn = forcefield_.findCharge(parser.argc(1));
							if (ixn != NULL)
							{
								printf("Error: Charge for atom '%s' has already been defined.\n", parser.argc(1));
								return FALSE;
							}

							// Create new charge
							ixn = forcefield_.addCharge();
							ixn->setAtoms(parser.argc(1));
							ixn->setType(Interaction::ChargeType);
							
							// Is it the dependent charge
							dependent = strcmp(parser.argc(2), "dependent") == 0;

							par = forcefield_.addParameter(dependent ? Parameter::DependentParameter : parameterType, Parameter::Moderate);
							par->set("q", dependent ? 0.0 : parser.argd(2));
							ixn->addParameter(Interaction::ChargeForm, par);
							if (dependent && (!forcefield_.setDependentCharge(par))) return FALSE;
							
							// Were limits given?
							if (parser.hasArg(3) && (strchr(parser.argc(3),',') != NULL))
							{
								// If this is a dependent parameter, then conplain
								if (dependent)
								{
									printf("Error: Limits can not be specified for dependent parameters.\n");
									return FALSE;
								}

								// Split arg around ','
								if (strchr(parser.argc(3),',') == NULL)
								{
									printf("Error:: Parameter range must be in the format \"[min,max]\" ('min' or 'max' may be omitted to indicate no limit, but the comma must still be present).\n");
									return FALSE;
								}
								
								Dnchar lowerLimit = beforeChar(parser.argc(3), ',');
								Dnchar upperLimit = afterChar(parser.argc(3), ',');
								if (!lowerLimit.isEmpty()) par->setLowerLimit(lowerLimit.asDouble());
								if (!upperLimit.isEmpty()) par->setUpperLimit(upperLimit.asDouble());
								// Increase arg counter
								++n;
							}
							else if (parser.hasArg(3))
							{
								printf("Error: Spurious extra data on charge line.\n");
								return FALSE;
							}
							break;
						// Unrecognised keyword
						case (TargetSystem::nInputKeywords):
							printf("Error: Unrecognised keyword '%s' found in '%s' block.\n", parser.argc(0), parameterType == Parameter::FixedParameter ? "fix" : "fit");
							return FALSE;
						// End
						case (TargetSystem::EndKeyword):
							done = TRUE;
							break;
						// Misplaced keyword
						default:
							printf("Error: Keyword '%s' used out of context.\n", parser.argc(0));
							return FALSE;
					}
					if (parser.eofOrBlank())
					{
						if (!done)
						{
							printf("Error: Unterminated '%s' block.\n", parameterType == Parameter::FixedParameter ? "fix" : "fit");
							return FALSE;
						}
						done = TRUE;
					}
				} while (!done);
				break;
			// Configurations
			case (TargetSystem::ConfigurationsKeyword):
				// Loop over items in block
				done = FALSE;
				forceMultiplier = 1.0;
				cellType = 0;
				do
				{
					parser.getArgsDelim();
					kwd = inputKeyword(parser.argc(0));
					switch (kwd)
					{
						// Cubic cell specification
						case (TargetSystem::CubicKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'cubic' keyword expects exactly one argument.\n");
								return FALSE;
							}
							cell.setIdentity();
							cell[0] = parser.argd(1);
							cell[5] = cell[0];
							cell[10] = cell[0];
							cellType = 1;
							break;
						// Orthorhombic cell specification
						case (TargetSystem::OrthorhombicKeyword):
							if (parser.nArgs() != 4)
							{
								printf("Error: 'orthorhombic' keyword expects exactly three arguments.\n");
								return FALSE;
							}
							cell.setIdentity();
							cell[0] = parser.argd(1);
							cell[5] = parser.argd(2);
							cell[10] = parser.argd(3);
							cellType = 2;
							break;
						// Parallelepiped cell specification
						case (TargetSystem::ParallelepipedKeyword):
							if (parser.nArgs() != 7)
							{
								printf("Error: 'parallelepiped' keyword expects exactly six arguments.\n");
								return FALSE;
							}
							// Work in unit vectors. Assume that A lays along x-axis
							cell.setColumn(0,1.0,0.0,0.0,0.0);
							// Assume that B lays in the xy plane. Since A={1,0,0}, cos(gamma) equals 'x' of the B vector.
							temp = cos(parser.argd(6)/DEGRAD);
							cell.setColumn(1,temp,sqrt(1.0 - temp*temp),0.0,0.0);
							// The C vector can now be determined in parts.
							// It's x-component is equal to cos(beta) since {1,0,0}{x,y,z} = {1}{x} = cos(beta)
							cell.setColumn(2,cos(parser.argd(5)/DEGRAD),0.0,0.0,0.0);
							// The y-component can be determined by completing the dot product between the B and C vectors
							cell[9] = ( cos(parser.argd(4)/DEGRAD) - cell[4]*cell[8] ) / cell[5];
							// The z-component is simply the remainder of the unit vector...
							cell[10] = sqrt(1.0 - cell[8]*cell[8] - cell[9]*cell[9]);
							// Lastly, adjust these unit vectors to give the proper cell lengths
							cell.columnMultiply(0,parser.argd(1));
							cell.columnMultiply(1,parser.argd(2));
							cell.columnMultiply(2,parser.argd(3));
							cell.setColumn(3, 0.0, 0.0, 0.0, 1.0);
							cellType = 3;
							break;
						// Force multiplier keyword
						case (TargetSystem::ForceMultiplierKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'forcemultiplier' keyword expects exactly one argument.\n");
								return FALSE;
							}
							forceMultiplier = parser.argd(1);
							break;
						// File containing configurations
						case (TargetSystem::FileKeyword):
							// First argument is filename
							if (parser.nArgs() != 2)
							{
								printf("Error: 'file' keyword expects exactly one argument.\n");
								return FALSE;
							}
							// Read in file...
							if (strstr(parser.argc(1), ".xyzf"))
							{
								printf("Opening XYZ configurations file '%s'...\n", parser.argc(1));
								if (!readConfigurationsXYZ(parser.argc(1), cell, forceMultiplier, cellType, FALSE))
								{
									printf("Error reading configurations file '%s'.\n", parser.argc(1));
									return FALSE;
								}
							}
							else if (strstr(parser.argc(1), ".xyzvf"))
							{
								printf("Opening XYZ configurations file '%s'...\n", parser.argc(1));
								if (!readConfigurationsXYZ(parser.argc(1), cell, forceMultiplier, cellType, TRUE))
								{
									printf("Error reading configurations file '%s'.\n", parser.argc(1));
									return FALSE;
								}
							}
							else if (strstr(parser.argc(1), ".fm"))
							{
								printf("Opening FMatch3 configurations file '%s'...\n", parser.argc(1));
								if (!readConfigurationsFM3(parser.argc(1), cell, forceMultiplier, cellType))
								{
									printf("Error reading configurations file '%s'.\n", parser.argc(1));
									return FALSE;
								}
							}
							else
							{
								printf("Unrecognised extension on configurations file.\n");
								return FALSE;
							}
							break;
						// Unrecognised keyword
						case (TargetSystem::nInputKeywords):
							printf("Error: Unrecognised keyword '%s' found in 'configurations' block.\n", parser.argc(0));
							return FALSE;
						// End
						case (TargetSystem::EndKeyword):
							done = TRUE;
							break;
						// Misplaced keyword
						default:
							printf("Error: Keyword '%s' used out of context.\n", parser.argc(0));
							return FALSE;
					}
					if (parser.eofOrBlank())
					{
						if (!done)
						{
							printf("Error: Unterminated '%s' block.\n", parameterType == Parameter::FixedParameter ? "fix" : "fit");
							return FALSE;
						}
						done = TRUE;
					}
				} while (!done);
				break;
			// Define job strategy
			case (TargetSystem::StrategyKeyword):
				strategyStack.clear();
				strategyStack.add(&strategy_, FALSE);
				// Loop over items in block
				done = FALSE;
				do
				{
					parser.getArgsDelim();
					sf = StrategyStep::stepFunction(parser.argc(0));
					if ((sf == StrategyStep::nStepFunctions) && (strcmp(parser.argc(0),"end") == 0)) break;
					switch (sf)
					{
						// Clear current stored alpha list
						case (StrategyStep::Clear):
							step = new ClearStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(1);
							break;
						// End current loop
						case (StrategyStep::EndLoop):
							// Check the bool parameter of the stack tail, which flags a loop strategy node
							if (strategyStack.last()->data) strategyStack.removeLast();
							else
							{
								printf("Error: Tried to 'endloop' without a current loop.\n");
								return FALSE;
							}
							break;
						// Perform Genetic Algorithm Search
						case (StrategyStep::GeneticMinimise):
							step = new GeneticStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argi(2));
							else step->addParameter(100);
							if (parser.hasArg(3)) step->addParameter(parser.argd(3));
							else step->addParameter(0.1);
							if (parser.hasArg(4)) step->addParameter(parser.argd(4));
							else step->addParameter(0.2);
							if (parser.hasArg(5)) step->addParameter(parser.argd(5));
							else step->addParameter(0.02);
							break;
						// Perform grid search on parameters
						case (StrategyStep::GridMinimise):
							step = new GridMinimiseStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(5);
							if (parser.hasArg(2)) step->addParameter(parser.argd(2));
							else step->addParameter(0.02);
							break;
						// Create loop
						case (StrategyStep::Loop):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'loop' step requires number of iterations to be provided.\n");
								return FALSE;
							}
							loopStep = new LoopStep();
							strategyStack.last()->item->own(loopStep);
							strategyStack.add(loopStep->strategy(), TRUE);
							loopStep->addParameter(parser.argi(1));
							break;
						// Predict
						case (StrategyStep::Predict):
							step = new PredictStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(1);
							if (parser.hasArg(2)) step->addParameter(parser.argi(2));
							else step->addParameter(5);
							if (parser.hasArg(3)) step->addParameter(parser.argi(3));
							else step->addParameter(5);
							break;
						// Print parameters
						case (StrategyStep::Print):
							step = new PrintStep();
							strategyStack.last()->item->own(step);
							if (parser.nArgs() == 2)
							{
								int info = 0;
								if (strcmp(parser.argc(1), "current") == 0) info = 1;
								else if (strcmp(parser.argc(1), "best") == 0) info = 2;
								if (info == 0)
								{
									printf("Error: Unrecognised option to 'print'.\n");
									return FALSE;
								}
								step->addParameter(info);
							}
							break;
						// Randomise input forcefield
						case (StrategyStep::Randomise):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'randomise' step requires degree of randomisation to be provided.\n");
								return FALSE;
							}
							step = new RandomisationStep();
							strategyStack.last()->item->own(step);
							step->addParameter(parser.argd(1));
							break;
						// Perform parameter rattle
						case (StrategyStep::Rattle):
							step = new RattleStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argd(2));
							else step->addParameter(0.05);
							break;
						// Steepest Descent
						case (StrategyStep::SD):
							step = new SteepestDescentStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argd(2));
							else step->addParameter(1.0);
							break;
						// Select parameter set to optimise
						case (StrategyStep::Select):
							if (parser.nArgs() < 2)
							{
								printf("Error: 'select' step requires at least one parameter strength to be provided.\n");
								return FALSE;
							}
							step = new SelectStep();
							strategyStack.last()->item->own(step);
							for (n=1; n<parser.nArgs(); ++n)
							{
								// Is the argument a number?
								if (parser.argi(n) != 0) step->addParameter(-parser.argi(n));
								else 
								{
									Parameter::ParameterStrength ps = Parameter::parameterStrength(parser.argc(n));
									if (ps == Parameter::nParameterStrengths)
									{
										printf("Error: Unrecognised parameter strength '%s'.\n", parser.argc(n));
										return FALSE;
									}
									step->addParameter(ps);
								}
							}
							break;
						// Simplex Anneal
						case (StrategyStep::SimplexAnneal):
							step = new SimplexAnnealStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argi(2));
							else step->addParameter(100);
							if (parser.hasArg(3)) step->addParameter(parser.argd(3));
							else step->addParameter(1.0E-3);
							if (parser.hasArg(4)) step->addParameter(parser.argd(4));
							else step->addParameter(0.1);
							break;
						// Simplex Annealing 2
						case (StrategyStep::SimplexAnneal2):
							step = new SimplexAnneal2Step();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argi(2));
							else step->addParameter(10);
							if (parser.hasArg(3)) step->addParameter(parser.argd(3));
							else step->addParameter(10.0);
							if (parser.hasArg(4)) step->addParameter(parser.argd(4));
							else step->addParameter(0.1);
							if (parser.hasArg(5)) step->addParameter(parser.argd(5));
							else step->addParameter(0.1);
							break;
						// Simplex Minimise
						case (StrategyStep::SimplexMinimise):
							step = new SimplexMinimiseStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argi(2));
							else step->addParameter(100);
							if (parser.hasArg(3)) step->addParameter(parser.argd(3));
							else step->addParameter(1.0E-3);
							if (parser.hasArg(4)) step->addParameter(parser.argd(4));
							else step->addParameter(0.1);
							break;
						// Store current best alpha
						case (StrategyStep::Store):
							step = new StoreStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(1);
							break;
						// Super Simplex
						case (StrategyStep::SuperSimplex):
							step = new SuperSimplexStep();
							strategyStack.last()->item->own(step);
							if (parser.hasArg(1)) step->addParameter(parser.argi(1));
							else step->addParameter(100);
							if (parser.hasArg(2)) step->addParameter(parser.argi(2));
							else step->addParameter(0.1);
							break;
						case (StrategyStep::HiddenFeature):
							step = new HiddenFeatureStep();
							strategyStack.last()->item->own(step);
							break;
						// Unrecognised keyword
						default:
							printf("Error: Step '%s' not recognised.\n", parser.argc(0));
							return FALSE;
					}
					if (parser.eofOrBlank())
					{
						if (!done)
						{
							printf("Error: Unterminated 'strategy' block.\n");
							return FALSE;
						}
						done = TRUE;
					}
				} while (!done);
				break;
			// Setup system parameters
			case (TargetSystem::SystemKeyword):
				// Loop over items in block
				done = FALSE;
				do
				{
					parser.getArgsDelim();
					kwd = inputKeyword(parser.argc(0));
					switch (kwd)
					{
						// Flag for combination of missing vdw parameters
						case (TargetSystem::CombineKeyword):
							if (parser.nArgs() != 1)
							{
								printf("Error: 'combine' keyword expects zero arguments in this context.\n");
								return FALSE;
							}
							combineMissingVdw_ = TRUE;
							printf("Missing van der Waals interactions will be generated according to combination rules.\n");
							break;
						// Set interatomic cutoff distance
						case (TargetSystem::CutoffKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'cutoff' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							cutoff_ = parser.argd(1);
							break;
						// Set electrostatic conversion factor
						case (TargetSystem::EConvertKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'econvert' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							elecConvert_ = 1389.35444426359172669289 * parser.argd(1); 
							break;
						// Turn on equalisation of parameter length scales
						case (TargetSystem::EqualiseKeyword):
							if (parser.nArgs() != 1)
							{
								printf("Error: 'equalise' keyword expects no arguments in this context.\n");
								return FALSE;
							}
							equalise_ = TRUE;
							break;
						// Set electronvolts as main unit of energy
						case (TargetSystem::ElectronVoltsKeyword):
							if (parser.nArgs() != 1)
							{
								printf("Error: 'ev' keyword expects no arguments in this context.\n");
								return FALSE;
							}
							energyUnit_ = TargetSystem::ElectronVolt;
							break;
						// Set Ewald sum parameters (and turn on Ewald sum)
						case (TargetSystem::EwaldKeyword):
							if (parser.nArgs() != 5)
							{
								printf("Error: 'ewald' keyword expects exactly four arguments in this context.\n");
								return FALSE;
							}
							ewaldAlpha_ = parser.argd(1);
							ewaldKMax_.set(parser.argi(2), parser.argi(3), parser.argi(4));
							useEwaldSum_ = TRUE;
							break;
						// Set frame range
						case (TargetSystem::FramesKeyword):
							if (parser.nArgs() < 3)
							{
								printf("Error: 'frames' keyword expects at least two arguments in this context.\n");
								return FALSE;
							}
							firstConfiguration_ = parser.argi(1)-1;
							lastConfiguration_ = parser.argi(2)-1;
							if (parser.hasArg(3)) configurationSkip_ = parser.argi(3);
							if (parser.hasArg(4)) configurationStackSize_ = parser.argi(4);
							break;
						// Set that no atom pair list should be used
						case (TargetSystem::NoPairsKeyword):
							if (parser.nArgs() != 1)
							{
								printf("Error: 'nopairs' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							useAtomPairs_ = FALSE;
							break;
						// Set parameter penalty factor
						case (TargetSystem::PenaltyKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'penalty' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							penaltyFactor_ = parser.argd(1);
							break;
						// Set unit cell replication values
						case (TargetSystem::ReplicateKeyword):
							if (parser.nArgs() != 4)
							{
								printf("Error: 'replicate' keyword expects three arguments in this context.\n");
								return FALSE;
							}
							rep.set(parser.argi(1), parser.argi(2), parser.argi(3));
							replicate = TRUE;
							break;
						// Set parameter retention flag
						case (TargetSystem::RetainKeyword):
							if (parser.nArgs() != 1)
							{
								printf("Error: 'retain' keyword expects no arguments in this context.\n");
								return FALSE;
							}
							retain_ = TRUE;
							break;
						// Set temperature schedule
						case (TargetSystem::ScheduleKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'schedule' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							if (TargetSystem::temperatureSchedule(parser.argc(1)) == TargetSystem::nTemperatureSchedules)
							{
								printf("Error: '%s' is not the name of a valid temperature schedule.\n", parser.argc(1));
								return FALSE;
							}
							schedule_ = TargetSystem::temperatureSchedule(parser.argc(1));
							break;
						// Set random seed
						case (TargetSystem::SeedKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'seed' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							seed_ = parser.argi(1);
							break;
						// Set cost threshold
						case (TargetSystem::ThresholdKeyword):
							if (parser.nArgs() != 2)
							{
								printf("Error: 'threshold' keyword expects exactly one argument in this context.\n");
								return FALSE;
							}
							threshold_ = parser.argd(1);
							break;
						// Unrecognised keyword
						case (TargetSystem::nInputKeywords):
							printf("Error: Unrecognised keyword '%s' found in 'configurations' block.\n", parser.argc(0));
							return FALSE;
						// End
						case (TargetSystem::EndKeyword):
							done = TRUE;
							break;
						// Misplaced keyword
						default:
							printf("Error: Keyword '%s' used out of context.\n", parser.argc(0));
							return FALSE;
					}
					if (parser.eofOrBlank())
					{
						if (!done)
						{
							printf("Error: Unterminated 'system' block.\n");
							return FALSE;
						}
						done = TRUE;
					}
				} while (!done);
				break;
			default:
				printf("Error: Keyword '%s' is not valid at this point in the file.\n", parser.argc(0));
				return FALSE;
		}
		
	} while (!parser.eofOrBlank());
	parser.closeFile();
	
	// Setup system - perform all checks and construct relevant forcefield terms and lists
	printf("\nChecking and setting up....\n\n");
	Atom *i, *j;
	Intra *b;

	// Step 1 - check number of atoms inferred by species matches number of atoms in configurations, and replicate if requested
	nAtoms_ = 0;
	for (sp = species_.first(); sp != NULL; sp = sp->next) nAtoms_ += sp->nAtoms()*sp->nMolecules();
	printf("Total number of atoms inferred by species definitions = %i\n", nAtoms_);
	if (replicate && (cellType == 0))
	{
		printf("Error: Can't replicate a system without a proper cell definition.\n");
		return FALSE;
	}
	for (Configuration *cfg = configurations_.first(); cfg != NULL; cfg = cfg->next)
	{
		++n;
		if (cfg->nAtoms() != nAtoms_)
		{
			printf("Error: Configuration %i contains %i atoms, which differs from the above value.\n", n, cfg->nAtoms());
			return FALSE;
		}

		// Replicate atomic coordinates if requested
		if (replicate && (!cfg->replicate(rep, species_.first()))) return FALSE;
	}
	if (replicate)
	{
		// Increase number of molecules for each defined species
		int factor = rep.x*rep.y*rep.z;
		for (sp = species_.first(); sp != NULL; sp = sp->next) sp->setNMolecules(sp->nMolecules()*factor);
		
		// Calculate new nAtoms_
		nAtoms_ *= factor;
	}

	// Step 2 - construct scaling matrices for each species
	for (sp = species_.first(); sp != NULL; sp = sp->next)
	{
		printf("Next species is %s\n", sp->name());
		sp->createScalingMatrices();
		printf("Created intramolecular scaling matrices for species '%s'\n", sp->name());
	}

	// Step 3 - Add missing combined vdW parameters (if requested), check forms of combinedVdw interactions, and set up parameters
	if (combineMissingVdw_)
	{
		// Loop over aom pairs in all species, checking for relevant vdW interaction
		for (sp = species_.first(); sp != NULL; sp = sp->next)
		{
			for (i = sp->atoms(); i != NULL; i = i->next)
				for (j=i->next; j != NULL; j = j->next)
				{
					if (forcefield_.findVdw(i->name(), j->name())) continue;
					ixn = forcefield_.addVdw();
					ixn->setAtoms(i->name(), j->name());
					ixn->setType(Interaction::CombinedVdwType);
					printf("Added combined van der Waals interaction between '%s' and '%s' in species '%s'\n", i->name(), j->name(), sp->name());
				}
		}
	}
	if (!forcefield_.checkCombinedVdw()) return FALSE;

	// Step 4 - Build unique atom name (type) list, assign type IDs back to atoms, and construct van der Waals interaction matrix
	for (sp = species_.first(); sp != NULL; sp = sp->next)
	{
		// Loop over atoms in species
		for (i = sp->atoms(); i != NULL; i = i->next)
		{
			for (j = uniqueAtoms_.first(); j != NULL; j = j->next) if (strcmp(i->name(), j->name()) == 0) break;
			// If no match found, add a new entry. Then, assign type id of new (or matched) atom to 'i'
			if (j == NULL)
			{
				j = uniqueAtoms_.add();
				j->setName(i->name());
				j->setTypeId(uniqueAtoms_.nItems()-1);
			}
			// Assign type id
			i->setTypeId(j->typeId());
		}
	}
	printf("\nConstructed unique type list:\n");
	for (j = uniqueAtoms_.first(); j != NULL; j = j->next) printf("\t%i\t%s\n", j->typeId(), j->name());
	printf("\nBuilding van der Waals interaction matrix...\n");
	if (!forcefield_.createInteractionMatrix(uniqueAtoms_)) return FALSE;
	forcefield_.printInteractionMatrix();
	printf("Creating type map...\n");
	typeMap_ = new int[nAtoms_];
	n = 0;
	for (sp = species_.first(); sp != NULL; sp = sp->next)
	{
		for (int m=0; m<sp->nMolecules(); ++m)
		{
			for (i = sp->atoms(); i != NULL; i = i->next) typeMap_[n++] = i->typeId();
		}
	}

	// Step 5 - check that all necessary parameters have been supplied to forcefield interactions
	if (!forcefield_.setup()) return FALSE;

	// Step 6 - check that charges have been supplied for all atoms (if useEwaldSum_)
	if (useEwaldSum_)
	{
		for (sp = species_.first(); sp != NULL; sp = sp->next)
		{
			// Loop over atoms in species
			for (i = sp->atoms(); i != NULL; i = i->next)
			{
				// Assign charge parameter corresponding to name
				if (!forcefield_.assignCharge(i, sp->nMolecules()))
				{
					printf("Error: Charge for atom '%s' in species '%s' has not been provided in the forcefield.\n", i->name(), sp->name());
					return FALSE;
				}
			}
		}
	}
	if (!forcefield_.checkChargeSetup()) return FALSE;

	// Step 7 - link interactions to intramolecular terms defined in species (and check atom ids)
	for (sp = species_.first(); sp != NULL; sp = sp->next)
	{
		Interaction *ixn;
		// Bonds
		for (Intra *bond = sp->bonds(); bond != NULL; bond = bond->next)
		{
			// Check atom ids are within range
			for (n=0; n<2; ++n) if (bond->atom(n) < 0 || bond->atom(n) >= sp->nAtoms())
			{
				printf("Error: Atom index specified in bond for species '%s' is out of range (%i vs %i atoms in single unit)\n", sp->name(), bond->atom(n)+1, sp->nAtoms());
				return FALSE;
			}
			// Search for matching interaction definition
			if (!forcefield_.applyBond(bond, sp->atom(bond->i()), sp->atom(bond->j())))
			{
				printf("Error: Failed to find bond definition for %s-%s specified in species '%s'.\n", sp->atom(bond->i())->name(), sp->atom(bond->j())->name(), sp->name());
				return FALSE;
			}
		}
		// Angles
		for (Intra *angle = sp->angles(); angle != NULL; angle = angle->next)
		{
			// Check atom ids are within range
			for (n=0; n<3; ++n) if (angle->atom(n) < 0 || angle->atom(n) >= sp->nAtoms())
			{
				printf("Error: Atom index specified in angle for species '%s' is out of range (%i vs %i atoms in single unit)\n", sp->name(), angle->atom(n)+1, sp->nAtoms());
				return FALSE;
			}
			// Search for matching interaction definition
			if (!forcefield_.applyAngle(angle, sp->atom(angle->i()), sp->atom(angle->j()), sp->atom(angle->k())))
			{
				printf("Error: Failed to find angle definition for %s-%s-%s specified in species '%s'.\n", sp->atom(angle->i())->name(), sp->atom(angle->j())->name(), sp->atom(angle->k())->name(), sp->name());
				return FALSE;
			}
		}
		// Torsions
		for (Intra *torsion = sp->torsions(); torsion != NULL; torsion = torsion->next)
		{
			// Check atom ids are within range
			for (n=0; n<4; ++n) if (torsion->atom(n) < 0 || torsion->atom(n) >= sp->nAtoms())
			{
				printf("Error: Atom index specified in torsion for species '%s' is out of range (%i vs %i atoms in single unit)\n", sp->name(), torsion->atom(n)+1, sp->nAtoms());
				return FALSE;
			}
			// Search for matching interaction definition
			if (!forcefield_.applyTorsion(torsion, sp->atom(torsion->i()), sp->atom(torsion->j()), sp->atom(torsion->k()), sp->atom(torsion->l())))
			{
				printf("Error: Failed to find torsion definition for %s-%s-%s-%s specified in species '%s'.\n", sp->atom(torsion->i())->name(), sp->atom(torsion->j())->name(), sp->atom(torsion->k())->name(), sp->atom(torsion->l())->name(), sp->name());
				return FALSE;
			}
		}
		printf("Found all necessary intramolecular terms for species '%s'.\n", sp->name());
	}

	// Step 8 - calculate reciprocal space cutoffs for all configurations, if Ewald sum is in use
	if (useEwaldSum_)
	{
		for (Configuration *cfg = configurations_.first(); cfg != NULL; cfg = cfg->next) cfg->calculateEwaldCutoff(ewaldKMax_);
	}

	// Step 9 - Create global scaling matrices for vdW and electrostatics
	printf("Creating global scaling matrices...\n");
	vdwScaleMatrix_ = new double*[nAtoms_];
	elecScaleMatrix_ = new double*[nAtoms_];
	for (n=0; n<nAtoms_; ++n)
	{
		vdwScaleMatrix_[n] = new double[nAtoms_];
		elecScaleMatrix_[n] = new double[nAtoms_];
		// Set intial elements here
		for (m=0; m<nAtoms_; ++m)
		{
			vdwScaleMatrix_[n][m] = 1.0;
			elecScaleMatrix_[n][m] = 1.0;
		}
		// Set diagonals to zero
		vdwScaleMatrix_[n][n] = 0.0;
		elecScaleMatrix_[n][n] = 0.0;
	}
	// Fill in intramolecular parts
	n = 0;
	for (sp = species_.first(); sp != NULL; sp = sp->next)
	{
		for (int m=0; m<sp->nMolecules(); ++m)
		{
			for (int ii = 0; ii < sp->nAtoms(); ++ii)
			{
				for (int jj = 0; jj < sp->nAtoms(); ++jj)
				{
					vdwScaleMatrix_[n+ii][n+jj] = sp->vdwScalingFactor(ii,jj);
					elecScaleMatrix_[n+ii][n+jj] = sp->elecScalingFactor(ii,jj);
				}
			}
			n += sp->nAtoms();
		}
	}

	// Print summary of species and parameter information before we exit
	n = 1;
	for (sp = species_.first(); sp != NULL; sp = sp->next)
	{
		printf("\nSpecies %i : %s\n", n, sp->name());
		m = 0;
		for (i = sp->atoms(); i != NULL; i = i->next) printf("\tAtom    %-3i  %-5s  TypeId = %i\n", ++m, i->name(), i->typeId());
		for (b = sp->bonds(); b != NULL; b = b->next) printf("\tBond    %-3i %-3i %-5s %-5s %s\n", b->i()+1, b->j()+1, sp->atom(b->i())->name(), sp->atom(b->j())->name(), b->interaction()->info());
		for (b = sp->angles(); b != NULL; b = b->next) printf("\tAngle   %-3i %-3i %-3i %-5s %-5s %-5s %s\n", b->i()+1, b->j()+1, b->k()+1, sp->atom(b->i())->name(), sp->atom(b->j())->name(), sp->atom(b->k())->name(), b->interaction()->info());
		for (b = sp->torsions(); b != NULL; b = b->next) printf("\tTorsion %-3i %-3i %-3i %-3i %-5s %-5s %-5s %-5s %s\n", b->i()+1, b->j()+1, b->k()+1, b->l()+1, sp->atom(b->i())->name(), sp->atom(b->j())->name(), sp->atom(b->k())->name(), sp->atom(b->l())->name(), b->interaction()->info());
	}

	// Print parameter information
	forcefield_.calculateDependencies();
	printf("\nParameter Info:\n");
	printf("\t%i free parameters over all forcefield terms:\n", forcefield_.nFreeParameters());
	forcefield_.printFreeParameters();
	printf("\t%i dependent parameters over all forcefield terms:\n", forcefield_.nDependentParameters());
	forcefield_.printDependentParameters();
	printf("\t%i fixed parameters over all forcefield terms:\n", forcefield_.nFixedParameters());
	forcefield_.printFixedParameters();

	// Check configuration range
	if ((firstConfiguration_ < 0) || (firstConfiguration_ >= configurations_.nItems()))
	{
		printf("Error: Starting configuration index %i is out of range (allowed values are 1 to %i).\n", firstConfiguration_+1, configurations_.nItems());
		return FALSE;
	}
	if (lastConfiguration_ == -1) lastConfiguration_ = configurations_.nItems()-1;
	if ((lastConfiguration_ >= configurations_.nItems()) || (lastConfiguration_ < firstConfiguration_))
	{
		printf("Error: Final configuration index %i is out of range (allowed values are %i to %i).\n", lastConfiguration_+1, firstConfiguration_+1, configurations_.nItems());
		return FALSE;
	}
	if (configurationSkip_ < 1)
	{
		printf("Error: Skip between configuration frames is less than 1 (%i).\n", configurationSkip_);
		return FALSE;
	}
	else if ((lastConfiguration_ != firstConfiguration_) && ((lastConfiguration_-firstConfiguration_)/configurationSkip_ == 0))
	{
		printf("\nWARNING: Configuration skip value is set such that only one configuration in the specified range will be considered.\n");
	}
	
	// Write summary of system control parameters
	printf("\nConfiguration range is %i to %i, with a skip value of %i\n", firstConfiguration_+1, lastConfiguration_+1, configurationSkip_);
	printf("Configurations will be considered (stacked) %i at a time\n", configurationStackSize_);
	printf("Interatomic cutoff distance (Angstroms)          = %12.6f\n", cutoff_);
	printf("Electrostatic conversion factor ()               = %f\n", elecConvert_);
	printf("Equalisation of parameter length scales          = %s\n", equalise_ ? "On" : "Off" );
	if (useEwaldSum_)
	{
		printf("Calculation method for electrostatics            = Ewald\n");
		printf("    Ewald alpha                                  = %12.5e\n", ewaldAlpha_);
		printf("    Ewald kMax                                   = %i,%i,%i\n", ewaldKMax_.x, ewaldKMax_.y, ewaldKMax_.z); 
	}
	else printf("Calculation method for electrostatics            = None\n");
	printf("Penalty multiplier for parameters outside limits = %12.5e\n", penaltyFactor_); 
	printf("Parameter retention between configuration stacks = %s\n", retain_ ? "On" : "Off" );
	printf("Temperature schedule for annealing runs          = %s\n", TargetSystem::temperatureSchedule(schedule_));
	printf("Random seed                                      = %u\n", seed_); 
	printf("Threshold cost at which to stop optimisation     = %11.5e\n", threshold_);
	printf("Assumed unit of energy                           = %s\n", energyUnit_ == TargetSystem::KiloJouleMol ? "kJ mol-1" : "eV");

	// Store original free forcefield parameters so we can revert back to them if necessary
	originalAlpha_.set(forcefield_.freeParameters());
	
	setup_ = TRUE;

	return TRUE;
}
