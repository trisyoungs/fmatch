/*
	*** System Definition
	*** system.cpp
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
#include "simplex.h"
#include "atompair.h"
#include "alphastore.h"
#include "md.h"
#include <string.h>
#include <time.h>
#include <ctime>

// Temperature schedule types
const char *TemperatureScheduleKeywords[TargetSystem::nTemperatureSchedules] = { "linear", "hilly" };
TargetSystem::TemperatureSchedule TargetSystem::temperatureSchedule(const char *s)
{
	for (int ts=0; ts < TargetSystem::nTemperatureSchedules; ++ts) if (strcmp(s, TemperatureScheduleKeywords[ts]) == 0) return (TargetSystem::TemperatureSchedule) ts;
	return TargetSystem::nTemperatureSchedules;
}
const char *TargetSystem::temperatureSchedule(TargetSystem::TemperatureSchedule ts)
{
	return TemperatureScheduleKeywords[ts];
}

// Constructor
TargetSystem::TargetSystem()
{
	// Private variables
	typeMap_ = NULL;
	setup_ = FALSE;
	seed_ = (unsigned)time( NULL );
	cutoff_ = 15.0;
	useEwaldSum_ = FALSE;
	ewaldAlpha_ = 0.0;
	elecConvert_ = 1389.35444426359172669289;
	combineMissingVdw_ = FALSE;
	firstConfiguration_ = 0;
	lastConfiguration_ = -1;
	configurationSkip_ = 1;
	configurationStackSize_ = 1;
	penaltyFactor_ = 1000.0;
	equalise_ = FALSE;
	retain_ = FALSE;
	useAtomPairs_ = TRUE;
	useStoredForces_ = TRUE;
	useStaticKVectors_ = TRUE;
	useEam_ = FALSE;
	schedule_ = TargetSystem::LinearSchedule;
	threshold_ = 100.0;
	energyUnit_ = TargetSystem::KiloJouleMol;
	md.setTargetSystem(this);

	abort_ = FALSE;
}

// Destructor
TargetSystem::~TargetSystem()
{
	if (typeMap_ != NULL) delete[] typeMap_;
}

/*
// Species Definition
*/

 // Return first item in species list
Species *TargetSystem::species()
{
	return species_.first();
}

// Return atom vdW type name for specified atom
const char *TargetSystem::atomTypeName(int n)
{
	return uniqueAtoms_[typeMap_[n]]->name();
}

/*
// Configurations
*/

// Prepare specified configuration for energy/force calculation
bool TargetSystem::addConfigurationToStack(Configuration *cfg)
{
	int a1, a2, m1, m2, i, j, n, pass, m2finish;
	int nAtomPairs, nExcludedPairs;
	AtomPair *atomPairs = NULL, *excludedPairs = NULL;
	Species *sp1, *sp2;
	double escale, vscale, rij;
	Vec3<double> vij;
	
	printf("\nPreparing configuration '%s' for energy/force calculation.\n", cfg->name());
	configurationStack_.add(cfg);

	// Create atom pair list
	// First pass - determine size of list
	// Second pass - generate list
	if (useAtomPairs_) for (pass = 0; pass < 2; ++pass)
	{
		nAtomPairs = 0;
		nExcludedPairs = 0;
		for (sp1 = species_.first(); sp1 != NULL; sp1 = sp1->next)
		{
			a1 = sp1->firstAtom();
			// Loop over species 1 molecules
			for (m1 = 0; m1 < sp1->nMolecules(); ++m1)
			{
				// Loop over species 1 atoms
				for (i = 0; i<sp1->nAtoms(); ++i)
				{
					// Loop over species 2
					for (sp2 = sp1; sp2 != NULL; sp2 = sp2->next)
					{
						a2 = sp2->firstAtom() + sp2->nAtoms() * (sp1 == sp2 ? m1 : 0);
						// Loop over species 2 molecules
						for (m2 = (sp1 == sp2 ? m1 : 0); m2 < sp2->nMolecules(); ++m2)
						{
							// Loop over species 2 atoms
							for (j = 0; j<sp2->nAtoms(); ++j)
							{
								// Calculate minimum image vector and distance
								vij = cfg->mimd(i+a1,j+a2);
								rij = vij.magnitude();
								// Check distance
								if (rij > cutoff_) continue;

								// Check intramolecular scale matrix for sp1 == sp2 && m1 == m2
								if (sp1 == sp2 && m1 == m2)
								{
									// Check for same atom
									if (i >= j) continue;
									escale = sp1->elecScalingFactor(i,j);
									vscale = sp1->vdwScalingFactor(i,j);
									// If scale factor is not 1.0, add to excluded list and continue
									if ((escale < 0.99 || vscale < 0.99))
									{
										if (pass == 1) excludedPairs[nExcludedPairs].set(a1+i,a2+j,sp1->atom(i), sp1->atom(j), vij, rij, escale, vscale);
										++nExcludedPairs;
									}
									if ((escale > 0.01 || vscale > 0.01))
									{
										if (pass == 1) atomPairs[nAtomPairs].set(a1+i,a2+j,sp1->atom(i), sp1->atom(j), vij, rij, escale, vscale);
										++nAtomPairs;
									}
									
									continue;
								}
								else
								{
									escale = 1.0;
									vscale = 1.0;
								}
								
								// Check van der Waals interaction definition
								if (forcefield_.interactionMatrix(sp1->atom(i)->typeId(),sp2->atom(j)->typeId()) == NULL) continue;

								// All is OK - Increment / add to list
								if (pass == 1) atomPairs[nAtomPairs].set(a1+i,a2+j,sp1->atom(i),sp2->atom(j), vij, rij, escale, vscale);
								++nAtomPairs;
							}
							a2 += sp2->nAtoms();
						}
					}
				}
				a1 += sp1->nAtoms();
			}
		}
		// End of first pass? If so, create array
		if (pass == 0)
		{
			cfg->createPairArrays(nAtomPairs, nExcludedPairs);
			atomPairs = cfg->atomPairs();
			excludedPairs = cfg->excludedPairs();
		}
	}
	else cfg->updateDistanceMatrix();

	// Create kvector, cos and sin arrays for Ewald sum
	if (useEwaldSum_)
	{
		cfg->createKVectorArrays(ewaldKMax_, ewaldAlpha_, elecConvert_);
		// Can also use static sums here for speed
		if (useStaticKVectors_) cfg->createStaticKVectorSums(species_.first(), ewaldKMax_.max(), forcefield_.nCharges());
	}

	return TRUE;
}
	
// Prepare nth configuration for energy / force calculation
bool TargetSystem::addConfigurationToStack(int cfgId)
{
	if (cfgId < 0 || cfgId >= configurations_.nItems())
	{
		printf("BUG: Configuration index %i is out of range in addConfigurationToStack.\n", cfgId);
		return FALSE;
	}
	return addConfigurationToStack(configurations_[cfgId]);
}

// Clear stack
void TargetSystem::clearStack()
{
	for (Refitem<Configuration,double> *ri = configurationStack_.first(); ri != NULL; ri = ri->next)
	{
		ri->item->deleteEAMArrays();
		ri->item->deletePrecalculatedArrays();
	}
	configurationStack_.clear();
}

// Calculate forces for configurations on current stack
bool TargetSystem::stackForces()
{
	int a1, a2, m1, m2, i, j, k, l, n;
	Species *sp1, *sp2;
	Intra *bound;
	Interaction *ixn;
	AtomPair *pair;
	Configuration *cfg;
	double totalSOSE = 0.0;
	
	// Check for no configurations
	if (configurationStack_.nItems() == 0)
	{
		printf("Error: No configurations in stack when calculating energy and forces.\n");
		return FALSE;
	}

	// Calculate forcefield logpoints
	forcefield_.calculateLogPoints();

	// Calculate energy, forces, and SOSE in forces for all frames
	// Only calculate terms which need updating (based on current forcefield logPoints)
	for (Refitem<Configuration,double> *ri = configurationStack_.first(); ri != NULL; ri = ri->next)
	{
		cfg = ri->item;
		
		/*
		// Van der Waals
		*/
		if (cfg->logPoint(Configuration::VdwContribution) != forcefield_.logPoint(Configuration::VdwContribution))
		{
			cfg->reset(Configuration::VdwContribution);
			if (useStoredForces_) cfg->updateLogPoint(&forcefield_, Configuration::VdwContribution);
			if (useEam_)
			{
				// DANGEROUS - EAM potential should be first (and only) potential in list
				if (useAtomPairs_) cfg->calculateEAMForces(forcefield_.vdw());
				else cfg->calculateEAMForcesQuick(forcefield_.vdw());
			}
			else if (useAtomPairs_)
			{
				pair = cfg->atomPairs();
				// 	printf("NAtom pairs to calculate = %i\n", nAtomPairs_);
				for (n=0; n<cfg->nAtomPairs(); ++n)
				{
					ixn = forcefield_.interactionMatrix(pair->typeI(), pair->typeJ());
					if (ixn == NULL) continue;
					if (!ixn->calculate(cfg, pair)) return FALSE;
					++pair;
				}
			}
			else
			{
				for (i=0; i<nAtoms_-1; ++i)
				{
					for (j=i+1; j<nAtoms_; ++j)
					{
						if (cfg->distance(i,j) > cutoff_) continue;
						ixn = forcefield_.interactionMatrix(typeMap_[i],typeMap_[j]);
						if (ixn == NULL) continue;
						if (!ixn->calculate(cfg, i, j, -1, -1, vdwScaleMatrix_[i][j])) return FALSE;
					}
				}
			}
		}

		/*
		// Electrostatics
		*/
		if (cfg->logPoint(Configuration::ChargeContribution) != forcefield_.logPoint(Configuration::ChargeContribution))
		{
			cfg->reset(Configuration::ChargeContribution);
			if (useStoredForces_) cfg->updateLogPoint(&forcefield_, Configuration::ChargeContribution);

			if (useEwaldSum_)
			{
				ewaldReal(cfg);
				if (cfg->hasStaticKVectorSums()) quickEwaldReciprocal(cfg);
				else ewaldReciprocal(cfg);
				ewaldCorrect(cfg);
			}
// 			else if (useCoulombSum_)
		}

		/*
		// Intramolecular term calculation
		*/
		bool calcBonds = FALSE, calcAngles = FALSE, calcTorsions = FALSE;
		if (cfg->logPoint(Configuration::BondContribution) != forcefield_.logPoint(Configuration::BondContribution))
		{
			cfg->reset(Configuration::BondContribution);
			if (useStoredForces_) cfg->updateLogPoint(&forcefield_, Configuration::BondContribution);
			calcBonds = TRUE;
		}
		if (cfg->logPoint(Configuration::AngleContribution) != forcefield_.logPoint(Configuration::AngleContribution))
		{
			cfg->reset(Configuration::AngleContribution);
			if (useStoredForces_) cfg->updateLogPoint(&forcefield_, Configuration::AngleContribution);
			calcAngles = TRUE;
		}
		if (cfg->logPoint(Configuration::TorsionContribution) != forcefield_.logPoint(Configuration::TorsionContribution))
		{
			cfg->reset(Configuration::TorsionContribution);
			if (useStoredForces_) cfg->updateLogPoint(&forcefield_, Configuration::TorsionContribution);
			calcTorsions = TRUE;
		}
		
		for (sp1 = species_.first(); sp1 != NULL; sp1 = sp1->next)
		{
			// Bonds
			if (calcBonds) for (bound = sp1->bonds(); bound != NULL; bound = bound->next)
			{
				i = bound->i();
				j = bound->j();
				a1 = sp1->firstAtom();
				// Loop over species molecules
				for (m1 = 0; m1 < sp1->nMolecules(); ++m1)
				{
					if (!bound->interaction()->calculate(cfg, a1+i, a1+j)) return FALSE;
					a1 += sp1->nAtoms();
				}
			}
			// Angles
			if (calcAngles) for (bound = sp1->angles(); bound != NULL; bound = bound->next)
			{
				i = bound->i();
				j = bound->j();
				k = bound->k();
				a1 = sp1->firstAtom();
				// Loop over species molecules
				for (m1 = 0; m1 < sp1->nMolecules(); ++m1)
				{
					if (!bound->interaction()->calculate(cfg, a1+i, a1+j, a1+k)) return FALSE;
					a1 += sp1->nAtoms();
				}
			}
			// Torsions
			if (calcTorsions) for (bound = sp1->torsions(); bound != NULL; bound = bound->next)
			{
				i = bound->i();
				j = bound->j();
				k = bound->k();
				l = bound->l();
				a1 = sp1->firstAtom();
				// Loop over species molecules
				for (m1 = 0; m1 < sp1->nMolecules(); ++m1)
				{
					if (!bound->interaction()->calculate(cfg, a1+i, a1+j, a1+k, a1+l)) return FALSE;
					a1 += sp1->nAtoms();
				}
			}
		}
		
		// Totalise force contributions for configuration (so fCalc_ are accurate)
		cfg->sumForceContributions();
	}
	return TRUE;
}
	
// Calculate SOSE for current configuration stack
double TargetSystem::stackSOSE()
{
	// Sum SOSE for each configuration in stack
	double totalSOSE = 0.0;
	Configuration *cfg;
	for (Refitem<Configuration,double> *ri = configurationStack_.first(); ri != NULL; ri = ri->next)
	{
		cfg = ri->item;
		
		// Increment SOSE
		totalSOSE += cfg->totalSOSE();
	}
	
	// Add on contributions from parameters going over their defined limits
	double penalty = penaltyFactor_ * forcefield_.freeParameterPenalty();
// 	printf("totalSOSE:  cost = %f, penalty = %f, total = %f\n", totalSOSE, penalty, totalSOSE+penalty);
	totalSOSE += penalty;
	return totalSOSE;
}

/*
// Job Strategy
*/

// Return whether system has been set up
bool TargetSystem::setup()
{
	return setup_;
}

// Return random seed
unsigned int TargetSystem::seed()
{
	return seed_;
}

// Set Ewald sum parameters
void TargetSystem::setEwald(double alpha, int kx, int ky, int kz)
{
	useEwaldSum_ = TRUE;
	ewaldAlpha_ = alpha;
	ewaldKMax_.set(kx,ky,kz);
}

// Return whether Ewald sum is in use
bool TargetSystem::useEwaldSum()
{
	return useEwaldSum_;
}

// Return electrostatic conversion factor
double TargetSystem::elecConvert()
{
	return elecConvert_;
}

// Return whether equalisation of parameter length scales is enabled
bool TargetSystem::equalise()
{
	return equalise_;
}

// Return whether to retain best set of parameters into next configuration cycle
bool TargetSystem::retain()
{
	return retain_;
}

// Return temperature multiplier based on supplied run coordinate and current schedule
double TargetSystem::temperatureMultiplier(double x)
{
	static double Pi = 3.141593;
	double result = 1.0;
	switch (schedule_)
	{
		case (TargetSystem::LinearSchedule):
			result = 1.0 - x;
			break;
		case (TargetSystem::HillySchedule):
			result = ((cos(x*5*Pi)+1.0)*(1.0-x)*0.5) * exp(-2.0*x);
			break;
		default:
			printf("Internal Error: No temperature schedule defined.\n");
			result = 1.0;
			break;
	}
	return result;
}

// Return threshold value for cost at which to stop minimisation
double TargetSystem::threshold()
{
	return threshold_;
}

// Return energy unit to assume throughout
TargetSystem::EnergyUnit TargetSystem::energyUnit()
{
	return energyUnit_;
}

/*
// Methods
*/

// Print job strategy
void TargetSystem::printStrategy()
{
	printf("\nPer-configuration job strategy consists of %i steps:\n\n", strategy_.nItems());
	int i = 0;
	for (StrategyStep *step = strategy_.first(); step != NULL; step = step->next)
	{
		printf("%-3i\t", ++i);
		step->print();
	}
}

// Return original fittable alpha
Alpha TargetSystem::originalAlpha()
{
	return originalAlpha_;
}

// Update best alpha (if cost of supplied alpha is lower)
bool TargetSystem::updateBestAlpha(Alpha* trialAlpha)
{
	double newCost = trialAlpha->cost(this, &forcefield_);
	// Is this the first update? If so, store passed alpha, else compare with current values
	if (bestAlpha_.nAlpha() == 0) printf("\t\t... %11.5e (init)\n", newCost);
	else if (newCost < bestAlpha_.cost(this, &forcefield_))
	{
		double perAtom = newCost/nAtoms_;
		double real = sqrt(perAtom/3.0);
		printf("\t\t... %11.5e  (%11.5e atom-1, %11.5e %s A-1 xyz-1)\n", newCost, perAtom, real, energyUnit_ == TargetSystem::KiloJouleMol ? "kJ mol-1" : "eV");
	}
	else return FALSE;

	// If we're here then we must need to update the bestAlpha_...
	// Poke trialAlpha cost into forcefield just in case...
	trialAlpha->poke();
	// Retrieve original alpha set...
	bestAlpha_ = originalAlpha_;
	// ...and read in current (now best) forcefield values
	bestAlpha_.peek();
	return TRUE;
}

// Return best alpha parameters found for current configuration set
Alpha TargetSystem::bestAlpha()
{
	return bestAlpha_;
}

// Clear stored alpha list
void TargetSystem::clearAlphaList(int id)
{
	if ((id < 0) || (id > 9))
	{
		printf("Internal Error: AlphaList id (%i) is out of range.\n", id);
		return;
	}
	alphaList_[id].clear();
}

// Add alpha to list
void TargetSystem::storeAlpha(int id, Alpha *alpha)
{
	if ((id < 0) || (id > 9))
	{
		printf("Internal Error: AlphaList id (%i) is out of range.\n", id);
		return;
	}
	Alpha *storedAlpha = alphaList_[id].insert(NULL);
	storedAlpha->copy(alpha);
// 	while (alphaList_[id].nItems() > alphaListLimit_) alphaList_[id].removeLast();
}

// Return stored list of Alpha
Alpha *TargetSystem::alphaList(int id)
{
	if ((id < 0) || (id > 9))
	{
		printf("Internal Error: AlphaList id (%i) is out of range.\n", id);
		return NULL;
	}
	return alphaList_[id].first();
}

// Return number of best alpha in alphalist
int TargetSystem::alphaListSize(int id)
{
	if ((id < 0) || (id > 9))
	{
		printf("Internal Error: AlphaList id (%i) is out of range.\n", id);
		return 0;
	}
	return alphaList_[id].nItems();
}

// Flag to abort ASAP
void TargetSystem::flagAbort()
{
	abort_ = TRUE;
}

// Return whether to abort ASAP
bool TargetSystem::abort()
{
	return abort_;
}

// Run Strategy on Configurations
int TargetSystem::run(TargetSystem::JobType job)
{
	double t;
	int n, loop, cfgNo;
	bool done = FALSE, foundBetter;
	int frameCount;
	Configuration *cfg;
	Alpha currentAlpha;

	// Print strategy (if job == Optimise)
	if (job == TargetSystem::Optimise) printStrategy();

	// If the jobtype is 'DummyRun' then quit here
	if (job == TargetSystem::DummyRun) return 0;

	// If jobType is 'MolecularDynamics' do it now...
	if (job == TargetSystem::MolecularDynamics)
	{
		// Turn off force storage here...
		useAtomPairs_ = FALSE;
		useStoredForces_ = FALSE;
		useStaticKVectors_ = FALSE;
		// Create single-configuration stack
		cfg = configurations_[0];
		if (!addConfigurationToStack(cfg)) return -1;
		// Initialise system
		md.initialise(cfg);
		// Set starting velocities
		md.generateVelocities();
		// Run MD
		md.run();
		return 0;
	}

	// Create statistics object for Alpha values
	AlphaStore statistics;

	// Main 'loop' over configurations
	cfgNo = 0;
	do
	{
		// Have we used up all available configurations?
		if (cfgNo > lastConfiguration_)
		{
			printf("\nNo more configurations available.\n\n");
			break;
		}
		
		// Select next configuration(s)
		printf("\n**************************************\nConfiguration(s) selected:\n");
		for (n=0; n<configurationStackSize_; ++n)
		{
			if (cfgNo > lastConfiguration_)
			{
				// Ran out of configurations...
				printf("Not enough frames left to fulful stack size requirement.\nExiting...\n");
				done = TRUE;
				break;
			}
			cfg = configurations_[cfgNo];
			printf("\t%5i\n%s\n", cfgNo+1, cfg->name());
			
			// Prepare configuration for calculation
			if (!addConfigurationToStack(cfg)) return -1;
			
			cfgNo += configurationSkip_;
		}
		if (done) break;
		printf("**************************************\n");
		
		// What are we supposed to be doing?
		switch (job)
		{
			// Perform optimisation
			case (TargetSystem::Optimise):
				// Retrieve original set of parameters, or retain best fit from last run?
				if (retain_ && (bestAlpha_.nAlpha() != 0))
				{
					currentAlpha = bestAlpha_;
					printf("\nBest parameters from the previous configuration run will be used as a starting point.\n");
				}
				else
				{
					currentAlpha = originalAlpha_;
					printf("\nStarting parameters reset to those in input forcefield.\n");
					// Clear global best alpha
					bestAlpha_.clear();
				}
				currentAlpha.poke();

				// Calculate relative length scales?
				if (equalise_) currentAlpha.calculateLengthScales(this, &forcefield_);
				else currentAlpha.resetLengthScales();

				// Now loop over strategy steps
				for (StrategyStep *step = strategy_.first(); step != NULL; step = step->next)
				{
					// Print out next step information and execute
					step->print();
					currentAlpha = step->execute(currentAlpha, &forcefield_, this);
// 					currentAlpha.print();

					// Abort?
					if (abort_) break;
					
					// Threshold reached?
					if (bestAlpha_.cost(this, &forcefield_) < threshold_)
					{
						printf("Cost threshold has been reached - current bestAlpha have cost of %11.5e (threshold = %11.5e)\n", bestAlpha_.cost(this, &forcefield_), threshold_);
						break;
					}
				}

				// Calculate and print current SOSE for stack
				stackForces();
				printf("\n\tFinal Stack SOSE(f) = %f\n", stackSOSE());
				frameCount = 1;
				for (Refitem<Configuration,double> *ri = configurationStack_.first(); ri != NULL; ri = ri->next, ++frameCount) printf("\t\t%2i   %f\n", frameCount, ri->item->totalSOSE());

				// Abort?
				if (abort_) break;

				// Accumulate optimised forcefield into statistics
				printf("[FINAL]\t"); bestAlpha_.print();
				
				// Append to output file
				static ofstream outputFile("fm3.best");
				if (outputFile.is_open() && outputFile.good())
				{
					Dnchar s(-1," %11.5e : ", bestAlpha_.cost(this, &forcefield_));
					for (int n=0; n<bestAlpha_.nAlpha(); ++n)
					{
						Parameter **params = bestAlpha_.sourceParameters();
						s.strcatf("%8.4e (%s) ", bestAlpha_.alpha(n), params[n] == NULL ? "???" : params[n]->name());
					}
					outputFile << s.get() << "\n";
				}
				else printf("!!!! Unable to write best parameters to output file.\n");
				
				bestAlpha_.poke();
				statistics.add(&bestAlpha_, this, &forcefield_);
				break;
			// Test forces
			case (TargetSystem::TestForces):
				stackForces();
				printf("\tStack SOSE(f) = %f\n", stackSOSE());
				for (Refitem<Configuration,double> *ri = configurationStack_.first(); ri != NULL; ri = ri->next) ri->item->print(this);
				break;
			// Test supplied parameters
			case (TargetSystem::TestParameters):
				currentAlpha = originalAlpha_;
				currentAlpha.peek();
				currentAlpha.print();
				break;
		}

		// Clear configuration stack
		clearStack();
	} while (!done);

	// Aborted?
	if (abort_)
	{
		printf("Optimisation aborted.\n");
		return -1;
	}

	// Print statistics
	statistics.print();
	
	return 0;
}
