/*
	*** System Definition
	*** system.h
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

#ifndef FM3_SYSTEM_H
#define FM3_SYSTEM_H

#include "species.h"
#include "configuration.h"
#include "interaction.h"
#include "strategy.h"
#include "ff.h"
#include "lineparser.h"
#include "md.h"

// Target System for Fitting Procedure
class TargetSystem
{
	public:
	// Constructor / Destructor
	TargetSystem();
	~TargetSystem();
	// Job type
	enum JobType { DummyRun, MolecularDynamics, Optimise, TestForces, TestParameters };
	// Temperature schedule types
	enum TemperatureSchedule { LinearSchedule, HillySchedule, nTemperatureSchedules };
	static TemperatureSchedule temperatureSchedule(const char *s);
	static const char *temperatureSchedule(TemperatureSchedule ts);
	// Energy units
	enum EnergyUnit { KiloJouleMol, ElectronVolt };


	/*
	// Input Keywords
	*/
	public:
	enum InputKeyword { AngleKeyword, AtomKeyword, BondKeyword, ChargeKeyword, CombineKeyword, ConfigurationsKeyword, CubicKeyword, CutoffKeyword, EConvertKeyword, EndKeyword, EwaldKeyword, EqualiseKeyword, ElectronVoltsKeyword, FileKeyword, FitKeyword, FixedKeyword, ForceMultiplierKeyword, FramesKeyword, MassesKeyword, NMolsKeyword, NoPairsKeyword, OrthorhombicKeyword, ParallelepipedKeyword, PenaltyKeyword, RelativeKeyword, ReplicateKeyword, RetainKeyword, ScheduleKeyword, SeedKeyword, SpeciesKeyword, StrategyKeyword, SystemKeyword, ThresholdKeyword, TorsionKeyword, VdwKeyword, nInputKeywords };
	static InputKeyword inputKeyword(const char *s);
	
	
	/*
	// FMatch3 Configurations File Keywords
	*/
	public:
	enum FMatch3Keyword { ConfigurationKeyword, AtomsKeyword, AtomsAndForcesKeyword, CellKeyword, EnergyKeyword, EnergyPerAtomKeyword, ESPKeyword, ForcesKeyword, nFMatch3Keywords };
	static FMatch3Keyword fmatch3Keyword(const char *s);


	/*
	// Species Definition
	// Description of Molecular/Atomic Species in System
	*/
	private:
	// Total number of atoms inferred by species
	int nAtoms_;
	// List of unique atom types defined within species
	List<Atom> uniqueAtoms_;
	// List of species in the system
	List<Species> species_;
	// Map of atom's VDW types
	int *typeMap_;
	// Global vdw and electrostatic scaling matrices
	double **vdwScaleMatrix_, **elecScaleMatrix_;

	public:
	// Return first item in species list
	Species *species();
	// Return atom vdW type name for specified atom
	const char *atomTypeName(int n);


	/*
	// Main Forcefield definition and Ewald routines
	*/
	private:
	// Main forcefield definitions
	Forcefield forcefield_;
	// Error function and complementary error function
	double erfc(double);
        double erf(double);
	// Calculate real-space Ewald sum energy and forces
	void ewaldReal(Configuration *cfg);
	// Calculate reciprocal-space Ewald sum energy and forces
	void ewaldReciprocal(Configuration *cfg);
	// Calculate quick reciprocal-space Ewald sum energy and forces (from stored cfg sums)
	void quickEwaldReciprocal(Configuration *cfg);
	// Calculate corrections to Ewald sum
	void ewaldCorrect(Configuration* cfg);


	/*
	// Configurations
	// List of configurations / frames containing coordinates and target forces
	*/
	private:
	// List of configurations to fit
	List<Configuration> configurations_;
	// Pointer to current fit target
	Reflist<Configuration,double> configurationStack_;
	// Configuration range values for fit, skip value, and stack size
	int firstConfiguration_, lastConfiguration_, configurationSkip_, configurationStackSize_;
	
	private:
	// Read in frames from XYZ file specified
	bool readConfigurationsXYZ(const char* filename, Matrix& cell, double forceMultiplier, int cellType, bool hasVelocities);
	// Read in frames from FMatch3 file specified
	bool readConfigurationsFM3(const char* filename, Matrix& cell, double forceMultiplier, int cellType);
	
	public:
	// Prepare specified configuration for energy/force calculation
	bool addConfigurationToStack(Configuration *cfg);
	// Prepare n'th configuration for energy/force calculation
	bool addConfigurationToStack(int n);
	// Clear stack
	void clearStack();
	// Calculate forces for configurations on current stack
	bool stackForces();
	// Calculate SOSE for current configuration stack
	double stackSOSE();


	/*
	// System Setup
	*/
	private:
	// Whether the system has been succesfully setup
	bool setup_;
	// Random seed
	unsigned int seed_;
	// Real-space cutoff for intermolecular atomic interactions
	double cutoff_;
	// Whether Ewald sum will be calculated
	bool useEwaldSum_;
	// Ewald kmax
	Vec3<int> ewaldKMax_;
	// Ewald alpha parameter
	double ewaldAlpha_;
	// Electrostatic conversion factor
	double elecConvert_;
	// Whether to generate missing van der Waals terms through combination rules
	bool combineMissingVdw_;
	// Penalty factor for parameters which go outside their defined limits
	double penaltyFactor_;
	// Whether equalisation of parameter length scales is enabled
	bool equalise_;
	// Whether to retain best set of parameters into next configuration cycle
	bool retain_;
	// Whether to use precalculated atom pair array or standard distance matrix
	bool useAtomPairs_;
	// Whether to use stored forces (if valid) to speed up calculation
	bool useStoredForces_;
	// Whether to use static kvector arrays to speed up Ewald sum
	bool useStaticKVectors_;
	// Whether EAM potential is used (only supports single atom type in this version)
	bool useEam_;
	// Temperature schedule to use in annealing runs
	TemperatureSchedule schedule_;
	// Threshold value for cost at which to stop minimisation
	double threshold_;
	// Energy unit to assume throughout
	EnergyUnit energyUnit_;

	private:
	// Add parameters to interaction
	bool addParameters(LineParser& parser, Interaction* ixn, Parameter::ParameterType type, Interaction::InteractionForm form, Parameter::ParameterStrength strength, int startArg);

	public:
	// Open specification from supplied file and perform all necessary checks on parameters
	bool open(const char *filename);
	// Return whether system has been set up correctly
	bool setup();
	// Return random seed
	unsigned int seed();
	// Set Ewald sum parameters
	void setEwald(double alpha, int kx, int ky, int kz);
	// Return whether Ewald sum is in use
	bool useEwaldSum();
	// Return electrostatic conversion factor
	double elecConvert();
	// Return whether equalisation of parameter length scales is enabled
	bool equalise();
	// Return whether to retain best set of parameters into next configuration cycle
	bool retain();
	// Return temperature multiplier based on supplied run coordinate and current schedule
	double temperatureMultiplier(double x);
	// Return threshold value for cost at which to stop minimisation
	double threshold();
	// Return energy unit to assume throughout
	EnergyUnit energyUnit();


	/*
	// Job Strategy
	*/
	private:
	// Minimisation strategy
	List<StrategyStep> strategy_;
	// List of all original fittable alpha values in current optimisation
	Alpha originalAlpha_;
	// Best set of alpha found for current configuration set
	Alpha bestAlpha_;
	// Lists of stored alpha values
	List<Alpha> alphaList_[10];
	// Flag to abort minimisation strategy at earliest available opportunity
	bool abort_;

	public:
	// Print job strategy
	void printStrategy();
	// Return original fittable alpha with original parameters
	Alpha originalAlpha();
	// Update best alpha (if cost of supplied alpha is lower)
	bool updateBestAlpha(Alpha *trialAlpha);
	// Return best alpha parameters found for current configuration set
	Alpha bestAlpha();
	// Clear stored alpha list
	void clearAlphaList(int id);
	// Add alpha to list
	void storeAlpha(int id, Alpha *alpha);
	// Return stored list of Alpha
	Alpha *alphaList(int id);
	// Return number of best alpha in alphalist
	int alphaListSize(int id);
	// Flag to abort ASAP
	void flagAbort();
	// Return whether to abort ASAP
	bool abort();
	// Perform requested task
	int run(TargetSystem::JobType job);


	/*
	// Molecular Dynamics
	*/
	public:
	// Molecular dynamics engine
	MD md;
};

#endif
