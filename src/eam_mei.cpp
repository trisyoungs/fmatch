/*
	*** Embedded atom model energy/forces
	*** eam.cpp
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

#include "configuration.h"

// Calculate EAM forces on atoms
bool Configuration::calculateEAMForces(Interaction *eamPotential)
{
	// Create arrays if they don't already exist
	int n, m, i, j, k;
	if (Sij_ == NULL)
	{
		Sij_ = new double*[nAtoms_];
		dSij_drij_ = new double*[nAtoms_];
		rho_ = new double[nAtoms_];
		drho_drij_ = new double*[nAtoms_];
		rhoa_ = new double*[nAtoms_];
		drhoa_drij_ = new double*[nAtoms_];
		Fi_ = new double[nAtoms_];
		dFi_drho_ = new double[nAtoms_];
		phi_ = new double*[nAtoms_];
		dphi_drij_ = new double*[nAtoms_];
		for (n=0; n<nAtoms_; ++n)
		{
			Sij_[n] = new double[nAtoms_];
			dSij_drij_[n] = new double[nAtoms_];
			drho_drij_[n] = new double[nAtoms_];
			rhoa_[n] = new double[nAtoms_];
			drhoa_drij_[n] = new double[nAtoms_];
			phi_[n] = new double[nAtoms_];
			dphi_drij_[n] = new double[nAtoms_];
		}
	}
	
	// Grab parameters from specified potential
	double Ec = eamPotential->parameter(0)->value();
	double phi0 = eamPotential->parameter(1)->value();
	double r0 = eamPotential->parameter(2)->value();
	double alpha = eamPotential->parameter(3)->value();
	double beta = eamPotential->parameter(4)->value();
	double gamma = eamPotential->parameter(5)->value();
	double delta = eamPotential->parameter(6)->value();
	double c[6];
	for (n=0; n<6; ++n) c[n] = eamPotential->parameter(7+n)->value();
	
	// Potential Form 2 - Mei et al.
	double x, omxsq, omxcubed, dphi, rc, rn, drhoa, rhoe = 1.0;
	double term, adivb, rdivr, dpart1, dpart2, rij, rpowk, rdivrpow, rdivrpowm1;
	Vec3<double> vij, fpair, fembed;
	static double roots[4] = { 0.0, 1.0,  1.414213562, 1.732050808 }, s[4] = { 0.0, 12.0, 6.0, 24.0 };
	static double base_rc = 1.95, base_rn = 1.75;
	static bool eam_cutoffs = TRUE;

	// Clear data arrays before we start..
	for (i=0; i<nAtoms_; ++i)
	{
		rho_[i] = 0.0;
		Fi_[i] = 0.0;
		dFi_drho_[i] = 0.0;
		for (j=0; j<nAtoms_; ++j)
		{
			drho_drij_[i][j] = 0.0;
			rhoa_[i][j] = 0.0;
		}
	}

	// Calculate cutoff function values for each atomic pair
	// [Mei et al, PRb, 43, 4653 (1991), eq 18a/18b]
	rc = base_rc * r0;
	rn = base_rn * r0;

	if (eam_cutoffs)
	{
		for (i=0; i<nAtoms_-1; ++i)
		{
			for (j=i+1; j<nAtoms_; ++j)
			{
// 				do i = 1,natoms-1
// 				do j = i+1,natoms
				vij = mimd(i, j);
				rij = vij.magnitude();
				x = (rij-rn)/(rc-rn);
				if (rij <= rn)
				{
					Sij_[i][j] = 1.0;
					dSij_drij_[i][j] = 0.0;
				}
				else if (rij >= rc)
				{
					Sij_[i][j] = 0.0;
					dSij_drij_[i][j] = 0.0;
				}
				else
				{
					omxsq = (1.0-x)*(1.0-x);
					omxcubed = omxsq*(1.0-x);
					Sij_[i][j] = omxcubed * (1.0 + 3.0*x + 6.0*x*x);
					dSij_drij_[i][j] = omxcubed * (3.0 + 12.0*x) - 3.0*omxsq * (1.0 + 3.0*x + 6.0*x*x);
					dSij_drij_[i][j] /= (rc-rn);
				}
				Sij_[j][i] = Sij_[i][j];
				dSij_drij_[j][i] = dSij_drij_[i][j];
// 				printf("XX %i %i %f\n", i+1, j+1, Sij_[i][j]);
			}
		}
	}
	else
	{
		for (i=0; i<nAtoms_; ++i)
			for (j=0; j<nAtoms_; ++j)
			{
				Sij_[i][j] = 1.0;
				dSij_drij_[i][j] = 0.0;
			}
	}

	// Calculate background electron density for each atom, and atomic densities for each atom pair neighbour distance
	//
	//  Eqs 2 / 3 : Background electronic density around individual atom i
	//
	//	rho_i = SUM  f(rij)			(Eq 2)
	//		j<>i
	//			 l   c(l) ( r0 )l
	//	f(rij) = rhoe * SUM  ---- ( -- )	(Eq 3)
	//		 	0,5   12  ( r  )
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vij = mimd(i, j);
			rij = vij.magnitude();
				
			rpowk = 1.0;
			for (k=0; k<6; ++k)
			{
				// SPEED Precalculate (r0/rij) power terms
				rhoa_[i][j] += rhoe * (c[k] / 12.0) * rpowk;
				drhoa_drij_[i][j] -= rhoe * (c[k] / 12.0) * k * rpowk / rij;
				rpowk *= (r0/rij);
			}
			rhoa_[j][i] = rhoa_[i][j];
			drhoa_drij_[j][i] = drhoa_drij_[i][j];

			rho_[i] += rhoa_[i][j] * Sij_[i][j];
			rho_[j] += rhoa_[j][i] * Sij_[j][i];

			drho_drij_[i][j] = drhoa_drij_[i][j]*Sij_[i][j] + rhoa_[i][j]*dSij_drij_[i][j];
			drho_drij_[j][i] = drho_drij_[i][j];
// 	printf("rho %i %i %f %f %f %f\n", i+1, j+1, rho_[i], rho_[j], Sij_[i][j], drho_drij_[i][j]);
		}
	}

	// Determine embedding energy for each atom
	//
	//  Eq 5 : Embedding energy per atom
	//		     [     alpha    ( rho  ) ] ( rho  ) alpha/beta
	//  F_i(rho(i)) = -Ec [ 1 - ----- ln ( ---- ) ] ( ---- )
	//		     [     beta     ( rhoe ) ] ( rhoe )
	//	
	//			      m
	//		+ 0.5 phi0 * SUM s(m) exp( -(sqrt(m)-1) * gamma )
	//			     1,3
	//
	//		  [ 				      delta    ( rho  ) ]
	//		* [ 1 + (sqrt(m)-1) * delta - sqrt(m) ----- ln ( ---- ) ]
	//		  [				      beta     ( rhoe ) ]
	//
	//		  ( rho  ) sqrt(m) * gamma / beta
	//		* ( ---- )
	//		  ( rhoe )
	adivb = alpha / beta;
	for (i=0; i<nAtoms_; ++i)
	{
		if (rho_[i] > 1.0e-10)
		{
			// SPEED Precalc powers, log
			rdivr = rho_[i] / rhoe;
			rdivrpow = pow(rdivr, adivb);
			rdivrpowm1 = pow(rdivr, adivb-1.0);
			Fi_[i] = -Ec * (1.0 - adivb * log(rdivr)) * rdivrpow;
			dFi_drho_[i] = (-Ec*(1.0-adivb*log(rdivr))*(adivb*rdivrpowm1)) + (Ec*adivb*rdivrpow*(rhoe/rho_[i]));
// 			printf("%i %f %f\n", i+1, Fi_[i], dFi_drho_[i]);
			for (m=1; m<4; ++m)
			{
				term = 0.5 * phi0 * s[m] * exp(-(roots[m]-1.0)*gamma);
				Fi_[i] += term * (1.0 + (roots[m]-1.0)*delta - roots[m]*(delta/beta)*log(rdivr)) * pow(rdivr,roots[m]*gamma/beta);
				dpart1 = (1.0+(roots[m]-1.0)*delta-roots[m]*(delta/beta)*log(rdivr)) * (roots[m]*gamma/beta)*pow(rdivr,(roots[m]*gamma/beta)-1.0);
				dpart2 = -pow(rdivr,roots[m]*gamma/beta) * (roots[m]*delta/beta)*rhoe/rho_[i];
				dFi_drho_[i] += term * (dpart1 + dpart2);
			}
		}
		else
		{
			Fi_[i] = 0.0;
			// write(0,"(a,i4,a)") "WARNING - rhobar for atom ",i," is zero...."
			dFi_drho_[i] = 0.0;
		}
	}

	// Calculate pair potential between atoms
	//
	// Eq 4
	//			 [	     ( r      ) ]     [        ( r	) ]
	//   phi_ij(rij) = -phi0 [ 1 + delta ( -- - 1 ) ] exp [ -gamma ( -- - 1 ) ]
	//			 [	     ( r0     ) ]     [        ( r0	) ]
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vij = mimd(i, j);
			rij = vij.magnitude();
			rdivr = rij/r0;

			term = -phi0 * (1.0 + delta*(rdivr-1.0)) * exp(-gamma*(rdivr-1.0));
			phi_[i][j] = term * Sij_[i][j];
			phi_[j][i] = phi_[i][j];

			dphi = -phi0*(1.0+delta*(rdivr-1.0))*exp(-gamma*(rdivr-1.0))*(-gamma/r0)-phi0*exp(-gamma*(rdivr-1.0))*delta/r0;
			dphi_drij_[i][j] = dphi * Sij_[i][j] + term * dSij_drij_[i][j];
			dphi_drij_[j][i] = dphi_drij_[i][j];

		}
	}

	// Calculate forces on atoms
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vij = mimd(i, j);
			rij = vij.magnitude();
			vij /= rij;
			
			// From pair potential
// 			fvec(1:3) = -dphi_drij(i,j) * rij(i,j,1:3) / rij(i,j,4)
// 			eam_forces_phi(i,1:3) = eam_forces_phi(i,1:3) - fvec(1:3)
// 			eam_forces_phi(j,1:3) = eam_forces_phi(j,1:3) + fvec(1:3)
			fpair = vij * -dphi_drij_[i][j];
			fContributions_[Configuration::VdwContribution][i] -= fpair;
			fContributions_[Configuration::VdwContribution][j] += fpair;

			// From embedding function
// 			fvec(1:3) = -(dFi_drho(i) * drho_drij(i,j) + dFi_drho(j) * drho_drij(i,j)) * rij(i,j,1:3) / rij(i,j,4)
// 			eam_forces_F(i,1:3) = eam_forces_F(i,1:3) - fvec(1:3)
// 			eam_forces_F(j,1:3) = eam_forces_F(j,1:3) + fvec(1:3)
			fembed = vij * -(dFi_drho_[i] * drho_drij_[i][j] + dFi_drho_[j] * drho_drij_[i][j]);
			fContributions_[Configuration::VdwContribution][i] -= fembed;
			fContributions_[Configuration::VdwContribution][j] += fembed;
		}
	}

// 	eam_forces_total(:,:) = eam_forces_F(:,:) + eam_forces_phi(:,:)

	// Calculate total potential energy
	double eam_energy_phi = 0.0;
	double eam_energy_F = 0.0;
	
	for (i=0; i<nAtoms_; ++i) 
	{
		energies_[Configuration::VdwContribution] += Fi_[i];
		eam_energy_F += Fi_[i];
	}
	
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			energies_[Configuration::VdwContribution] += phi_[i][j];
			eam_energy_phi += phi_[i][j];
		}
	}
// 	printf("%f %f %f\n", mimd(0, 1).magnitude(), eam_energy_phi, eam_energy_F);
// 	eam_energy_total = eam_energy_phi + eam_energy_F
	return TRUE;
}
