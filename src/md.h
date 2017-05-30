/*
	*** Molecular Dynamics Engine
	*** md.h
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

#ifndef FM3_MD_H
#define FM3_MD_H

#include "vector3.h"

// Forward Declarations
class Forcefield;
class TargetSystem;
class Configuration;

// Molecular Dynamics Engine
class MD
{
	public:
	// Constructor / Destructor
	MD();
	~MD();


	/*
	// Basic Data
	*/
	private:
	// Supplied target system pointer
	TargetSystem *targetSystem_;
	// Target configuration
	Configuration *configuration_;
	// Number of atoms in simulation
	int nAtoms_;
	// Previous atomic positions (current positions are stored in Configuration)
	Vec3<double> *rLast_;
	// Atomic velocities and centre of mass velocity
	Vec3<double> *v_, vCom_;
	// Atomic accelerations
	Vec3<double> *a_;
	// Atomic masses
	double *mass_;
	// Requested and instantaneous simulation temperatures
	double t_, tInstant_;
	// Number of steps to run for
	int nSteps_;
	// Timestep
	double deltaT_;
	// Kinetic energy of system
	double eKinetic_;
	
	private:
	// Calculate kinetic energy of system
	void kinetic();
	// Propagate system using Verlet integration
	void verlet();
	// Propagate system using velocity Verlet integration
	void velocityVerletA();
	void velocityVerletB();
	// Print step summary
	void print(int step);

	public:
	// Set targetsystem
	void setTargetSystem(TargetSystem* system);
	// Initialise arrays etc.
	void initialise(Configuration *cfg);
	// Set basic run information
	void set(int nSteps, double temperature, double deltaT);
	// Generate initial velocities
	void generateVelocities();
	// Run simulation
	int run();
};

#endif
