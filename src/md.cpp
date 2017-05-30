/*
	*** Molecular Dynamics Engine
	*** md.cpp
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

#include "md.h"
#include "random.h"
#include "constants.h"
#include "configuration.h"
#include "system.h"

// Constructor
MD::MD()
{
	// Private variables
	targetSystem_ = NULL;
	rLast_ = NULL;
	nAtoms_ = 0;
	v_ = NULL;
	a_ = NULL;
	mass_ = NULL;
	nSteps_ = 100;
	t_ = 300.0;
	deltaT_ = 0.001;
}

// Destructor
MD::~MD()
{
	delete[] a_;
	delete[] rLast_;
	delete[] v_;
	delete[] mass_;
}

// Generate initial velocities
void MD::generateVelocities()
{
	double ke, massSum = 0.0;
	int i;
	// Initialise random atomic velocities according to requested temperature
	printf("Initialising particle velocities for T = %11.5e K\n", t_);
	vCom_.zero();
	for (i=0; i<nAtoms_; ++i)
	{
		v_[i].x = exp(Random::number()-0.5) / sqrt(TWOPI);
		v_[i].y = exp(Random::number()-0.5) / sqrt(TWOPI);
		v_[i].z = exp(Random::number()-0.5) / sqrt(TWOPI);
		vCom_ += v_[i] * mass_[i];
		massSum += mass_[i];
	}

	// Remove velocity shift
	vCom_ /= massSum;
	for (i=0; i<nAtoms_; ++i) v_[i] -= vCom_;

	// Rescale velocities for desired temperature
	kinetic();
	for (i=0; i<nAtoms_; ++i) v_[i] *= sqrt(t_ / tInstant_);
}

// Calculate kinetic energy of system
void MD::kinetic()
{
	// Determine total kinetic energy and instantaneous temperature of particles
	// As a check, also calculate total velocity of system (md_sumv)
	// KE = SUM(i) [ 0.5 * mass(i) * v(t)**2 ]
	// T = KE / ( 3/2 * N * kb)
	int i;
	eKinetic_ = 0.0;
	vCom_.zero();
	for (i = 0; i<nAtoms_; ++i)
	{
		eKinetic_ += 0.5 * mass_[i] * v_[i].dp(v_[i]);
		vCom_ += v_[i] * mass_[i];
	}

	// Finalise total kinetic energy (10J mol-1)
// 	eKinetic_ *= 0.5;

	// Calculate instantaneous temperature
	// J = kg m2 s-2  -->   10 J = g Ang2 ps-2
	// If ke is in units of [g mol-1 Angstroms2 ps-2] then must use kb in units of 10 J mol
 	tInstant_ = eKinetic_ * 2.0 / ( 3.0 * nAtoms_ * KB_10JMOLK);
	
	// Convert kinetic energy from 10 J/mol to desired units
	if (targetSystem_->energyUnit() == TargetSystem::ElectronVolt) eKinetic_ /= 9648.536;		// 10 J/mol -> eV
	else eKinetic_ /= 100.0;		// 10 J/mol -> kJ/mol
}

// Propagate system using Verlet integration
void MD::verlet()
{
	// Propagate system according to Verlet algorithm
	double deltaTSq_ = deltaT_*deltaT_, twoDeltaT_ = 2.0 * deltaT_;
	Vec3<double> rNew, a, *f, *r, force;
	f = configuration_->fCalc();
	r = configuration_->r();

	for (int i=0; i<nAtoms_; ++i)
	{
		// Get force (in correct units)
		force = f[i] * (targetSystem_->energyUnit() == TargetSystem::ElectronVolt ? 9648.536 : 100.0);

		// r(t+dt) = 2*r(t) - r(t-dt) + f(t)/m * dt**2
		a = force / mass_[i];
		rNew = r[i]*2.0 - rLast_[i] + a * deltaTSq_;
		
		// v(d+dt) = (r(t+dt) - r(t-dt)) / 2dt
		v_[i] = (rNew - rLast_[i]) / twoDeltaT_;
		
		// Store new values
		rLast_[i] = r[i];
		r[i] = rNew;
	}
	
// 	! Fold new atomic positions into the unit cell
// 	call md_foldatoms()
// 	! PBC old positions so we can calculate velocity properly on next iteration
// 	call md_pbcatoms(md_rlast, md_r)
// 	end subroutine md_verlet
}

// Propagate system using velocity Verlet integration (part A)
void MD::velocityVerletA()
{
	// Propagate system according to Verlet algorithm
	double deltaTSq = deltaT_*deltaT_;
	Vec3<double> *r;
	r = configuration_->r();

	// A:  r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt**2
	// A:  v(t+dt/2) = v(t) + 0.5*a(t)*dt
	// B:  a(t+dt) = F(t+dt)/m
	// B:  v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
	for (int i=0; i<nAtoms_; ++i)
	{
		// Propagate positions (by whole step)...
		r[i] += v_[i]*deltaT_ + a_[i]*0.5*deltaTSq;
		
		// ...velocities (by half step)...
		v_[i] += a_[i]*0.5*deltaT_;
	}
// 	! Fold new atomic positions into the unit cell
// 	call md_foldatoms()
}

// Propagate system using velocity Verlet integration (part A)
void MD::velocityVerletB()
{
	// Propagate system according to Verlet algorithm
	Vec3<double> force, *f = configuration_->fCalc();

	// A:  r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt**2
	// A:  v(t+dt/2) = v(t) + 0.5*a(t)*dt
	// B:  a(t+dt) = F(t+dt)/m
	// B:  v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
	for (int i=0; i<nAtoms_; ++i)
	{
		// Get force (in correct units of 10J/mol)
		force = f[i] * (targetSystem_->energyUnit() == TargetSystem::ElectronVolt ? 9648.536 : 100.0);

		// Determine new accelerations
		a_[i] = force / mass_[i];

		// ..and finally velocities again (by second half-step)
		v_[i] += a_[i]*0.5*deltaT_;
	}
}

// Print step summary
void MD::print(int step)
{
	if (step == 0)
	{
		printf("                                                                        Energies, eV\n");
		printf("Step      T, K        Total        Kinetic     Potential   P.E. / Atom   P.E. (emb)    P.E. (pair)    COM Velocity\n");
	}
	double ePotential = configuration_->totalEnergy();
	printf("%-6i %11.5e %12.5e %12.5e %12.5e %12.5e %8.5f %8.5f %8.5f\n", step, tInstant_, ePotential+eKinetic_, eKinetic_, ePotential, ePotential/nAtoms_, vCom_.x, vCom_.y, vCom_.z);
}

// Set targetsystem
void MD::setTargetSystem(TargetSystem *system)
{
	targetSystem_ = system;
}

// Initialise simulation, providing velocities at temperature requested
void MD::initialise(Configuration *cfg)
{
	// Initialise arrays
	configuration_ = cfg;
	nAtoms_ = cfg->nAtoms();
	rLast_ = new Vec3<double>[nAtoms_];
	v_ = new Vec3<double>[nAtoms_];
	a_ = new Vec3<double>[nAtoms_];
	mass_ = new double[nAtoms_];

	// Grab atomic masses
	int m, count = 0;
	Atom *i;
	for (Species *sp = targetSystem_->species(); sp != NULL; sp = sp->next)
	{
		for (m = 0; m<sp->nMolecules(); ++m)
		{
			for (i = sp->atoms(); i != NULL; i = i->next)
			{
				mass_[count] = i->mass();
				if (mass_[count] < 0.0) printf("Error: Mass for atom %i in system has not been set in species.\n", count);
				++count;
			}
		}
	}
}

// Set basic run information
void MD::set(int nSteps, double temperature, double deltaT)
{
	// Set variables
	nSteps_ = nSteps;
	t_ = temperature;
	deltaT_ = deltaT;
}

// Run simulation
int MD::run()
{
	Dnchar s;

	if (!targetSystem_->stackForces()) return -1;

	// Calculate previous atomic positions and accelerations (for use by basic Verlet and velocity algorithms)
	Vec3<double> *r = configuration_->r(), *f = configuration_->fCalc(), force;
	for (int i=0; i<nAtoms_; ++i)
	{
		rLast_[i] = r[i] - v_[i]*deltaT_;
		
		force = f[i] * (targetSystem_->energyUnit() == TargetSystem::ElectronVolt ? 9648.536 : 100.0);
		a_[i] = force / mass_[i];
	}

	// Calculate potential energy and forces, and then kinetic energy and temperature
	kinetic();
	print(0);

	ofstream outputFile, trajectoryFile;
	trajectoryFile.open("md.traj.xyzvf");
	outputFile.open("md.out");
	outputFile << "# Step    Total        Temperature  Kinetic      Potential    vdW          Elec         Bond         Angle        Torsion\n";
	
	trajectoryFile << nAtoms_ << "\n" << "Title  0\n";
	for (int i=0; i<nAtoms_; ++i)
	{
		s.sprintf("%-8s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n", targetSystem_->atomTypeName(i), r[i].x, r[i].y, r[i].z, v_[i].x, v_[i].y, v_[i].z, f[i].x, f[i].y, f[i].z);
		trajectoryFile << s.get();
	}
	
	for (int step = 1; step <=nSteps_; ++step)
	{
		velocityVerletA();
		configuration_->updateDistanceMatrix();
		if (targetSystem_->useEwaldSum()) configuration_->updateKVectorArrays();
		targetSystem_->stackForces();
// 		verlet();
		velocityVerletB();
		kinetic();
		print(step);
// 		printf("%f  %f  %f  %f  %f\n", r[1].x, configuration_->energy(Configuration::VdwContribution), f[1].x, f[1].y, f[1].z);
// 		r[1].print();
// 		f[0].print();
// 		f[1].print();
// 		r[1].x += 0.01;

		// Write trajectory xyz file
		if (step%10 == 0)
		{
			trajectoryFile << nAtoms_ << "\n" << "Title " << step << "\n";
			for (int i=0; i<nAtoms_; ++i)
			{
				s.sprintf("%-8s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n", targetSystem_->atomTypeName(i), r[i].x, r[i].y, r[i].z, v_[i].x, v_[i].y, v_[i].z, f[i].x, f[i].y, f[i].z);
				trajectoryFile << s.get();
			}
		}
		
		// Write output file
		if (step%10 == 0)
		{
			double ePotential = configuration_->totalEnergy();
			s.sprintf("%7i  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", step, ePotential+eKinetic_, tInstant_, eKinetic_, ePotential, configuration_->energy(Configuration::VdwContribution), configuration_->energy(Configuration::ChargeContribution), configuration_->energy(Configuration::BondContribution), configuration_->energy(Configuration::AngleContribution), configuration_->energy(Configuration::TorsionContribution));
			outputFile << s;
		}
	}
	outputFile.close();
	trajectoryFile.close();
	return 0;
}
