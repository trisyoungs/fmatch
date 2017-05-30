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
		dFi_drhoa_ = new double*[nAtoms_];
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
			dFi_drhoa_[n] = new double[nAtoms_];
		}
	}
	
	double rc, rn, rij, x, omxsq, omxcubed, dphi, astar, Eu, phibar, dEu, dphibar;
	Vec3<double> vij, fpair, fembed;
	static double base_rc = 1.95, base_rn = 1.75;
	
	// Grab parameters from specified potential
	double r0 = eamPotential->parameter(0)->value();
	double E0 = eamPotential->parameter(1)->value();
	double A0 = eamPotential->parameter(2)->value();
	double alpha = eamPotential->parameter(3)->value();
	double beta = eamPotential->parameter(4)->value();
	double rhoscale = eamPotential->parameter(5)->value();
	double Z0 = eamPotential->parameter(6)->value();
	
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
			dFi_drhoa_[i][j] = 0.0;
		}
	}

	// Calculate cutoff function values for each atomic pair
	// [Mei et al, PRb, 43, 4653 (1991), eq 18a/18b]
	rc = base_rc * r0;
	rn = base_rn * r0;

	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
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
		}
	}

	// Calculate background electron density for each atom, and atomic densities for each atom pair neighbour distance
	// Jelinek et al, Cond. Matt. arxiv:cond-mat/0610602v4
	//  Eq 4 : Background electronic density around individual atom i
	//		     rho_i^(0)
	//	 rhobar(i) = --------- * G(Gamma_i)
	//		      rho_i^0
	//
	//	rho_i^(k) = SUM  rho_j^a(k)(rij) * Sij			(Eq 9a, for k = 0)
	//		    j<>i
	//
	//	  rho_i^0 = rho_i0 * Z_i0 * G(Gamma_i^ref)		(Eq 7)
	//
	//				[		(  rij	    ) ]
	// rho_i^a(k)(rij) = rho_i0 * exp[ -beta_i^(k) * ( ----- - 1 ) ]	(Eq 10)
	//				[		( r_i^0	    ) ]
	//
	//	    r_i^0 = nearest-neighbour distance in reference structure (EAM parameter 1, r0)
	//      beta_i^(k) = element-dependent parameters (EAM parameter 5, beta)
	//	   rho_i0 = element-dependent density scaling (should be an EAM parameter 6, rhoscale)
	//	     Z_i0 = First nearest-neighbour coordination number of reference system (EAM parameter 7, Z0)
	//	      Sij = cutoff function between atom, calculated above
	//	Neglecting higher-order terms (i.e. only k = 0 considered), Gamma_i = 0 and G(Gamma_i) = 1 (also for Gamma_i^ref)

	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vij = mimd(i,j);
			rij = vij.magnitude();

			rhoa_[i][j] = rhoscale * exp(-beta*((rij/r0)-1.0));
			rhoa_[j][i] = rhoa_[i][j];

			drhoa_drij_[i][j] = rhoscale * exp( -beta*((rij/r0)-1.0) ) * (-beta/r0);
			drhoa_drij_[j][i] = drhoa_drij_[i][j];

			rho_[i] += rhoa_[i][j] * Sij_[i][j];
			rho_[j] += rhoa_[j][i] * Sij_[j][i];

			drho_drij_[i][j] = drhoa_drij_[i][j]*Sij_[i][j] + rhoa_[i][j]*dSij_drij_[i][j];
			drho_drij_[j][i] = drho_drij_[i][j];
		}
	}

	// If k = 0 only, then Eq 7 reduces to 'rho_i0 * Z_i0'
	for (i=0; i<nAtoms_; ++i)
	{
		rho_[i] /= rhoscale * Z0;
		for (j=0; j<nAtoms_; ++j) drho_drij_[i][j] /= rhoscale * Z0;
	}
// 	rhobar(:) = rhobar(:) / (rhoscale * Z0)
// 	drhobar_drij(:,:) = drhobar_drij(:,:) / (rhoscale * Z0)

	// Determine embedding energy for each atom
	// Jelinek et al, Cond. Matt. arxiv:cond-mat/0610602v4
	//
	//  Eq 3 : Embedding energy per atom
	//
	//  F_i(rhobar(i)) = Ai * Ei * rhobar(i) * ln(rhobar(i))
	//
	//	        Ai = Energy scaling parameter (EAM parameter 3, A0)
	//		Ei = Sublimation energy (EAM parameter 2, E0)
	for (i=0; i<nAtoms_; ++i)
	{
		if (rho_[i] > 1.0e-10)
		{
			Fi_[i] = A0 * E0 * rho_[i] * log(rho_[i]);
			dFi_drho_[i] = (A0*E0)*(1.0 + log(rho_[i]));
		}
		else
		{
			// Note: Derivative at rhobar(i) = 0 should be -INF
		// 	    write(0,"(a,i4,a)") "WARNING - rhobar for atom ",i," is zero...."
			dFi_drho_[i] = -(A0*E0);
		}
	}
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			if (rhoa_[i][j] > 1.0e-10) dFi_drhoa_[i][j] = A0 * E0 * (1.0 + log(rhoa_[i][j])) * drhoa_drij_[i][j];
			else dFi_drhoa_[i][j] = -(A0 * E0);
			dFi_drhoa_[j][i] = dFi_drhoa_[i][j];
		}
	}

	// Calculate pair potential between atoms
	// Jelinek et al, Cond. Matt. arxiv:cond-mat/0610602v4
	//
	// Eq 12 : Pair potentials
	//
	//    phi_ij(rij) = phibar_ij(rij) * Sij
	//
	//		    1  [		      ( Zij		    )	   ( Zij		 ) ]
	// phibar_ij(rij) = --- [ 2 * E_ij^u(rij) - Fi ( --- * rho_j^a0(rij) ) - Fj ( --- * rho_j^a0(rij) ) ]  (Eq 13)
	//		   Zij (		      ( Zi		    )      ( Zj		         ) ]
	//
	//	  E_ij^u = -Eij * (1 + a_ij(rij)) * exp(-a_ij(rij))						  (Eq 14)
	//
	//			     ( rij      )
	//	   a_ij = alpha_ij * ( ---- - 1 )								  (Eq 15)
	//			     ( rij0     )
	//
	//	   E_ij = Sublimation energy (EAM parameter 2 for pure metals, E0)
	//      alpha_ij = Exponential decay factor for universal energy term (EAM parameter 4, alpha)
	//	   Z_ij = "Depends on the structure of the reference system" (EAM parameter 7, Z0)
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			rij = mimd(i,j).magnitude();

			astar = alpha*((rij/r0) - 1.0);

			Eu = -E0 * (1.0 + astar) * exp(-astar);

			x = A0 * E0 * rhoa_[i][j] * log(rhoa_[i][j]);	// * (Zij / Zi)
			phibar = (2.0 * Eu - x - x) / Z0;
			
			phi_[i][j] = phibar * Sij_[i][j];
			phi_[j][i] = phi_[i][j];

			// Derivatives
			dEu = E0 * astar * exp(-astar) * (alpha / r0);
			x = dFi_drhoa_[i][j];	// * (Zij / Zi);
			dphibar = (2.0 * dEu - x - x) / Z0;

			dphi_drij_[i][j] = dphibar * Sij_[i][j] + phibar * dSij_drij_[i][j];
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
			fpair = vij * -dphi_drij_[i][j];
			fContributions_[Configuration::VdwContribution][i] -= fpair;
			fContributions_[Configuration::VdwContribution][j] += fpair;

			// From embedding function
// 			fvec(1:3) = -(dFi_drhobar(i) * drhobar_drij(i,j) + dFi_drhobar(j) * drhobar_drij(i,j)) * rij(i,j,1:3) / rij(i,j,4)
			fembed = vij * -(dFi_drho_[i] * drho_drij_[i][j] + dFi_drho_[j] * drho_drij_[i][j]);
			fContributions_[Configuration::VdwContribution][i] -= fembed;
			fContributions_[Configuration::VdwContribution][j] += fembed;
		}
	}

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

// Calculate EAM forces on atoms (using stored distance matrix)
bool Configuration::calculateEAMForcesQuick(Interaction *eamPotential)
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
		dFi_drhoa_ = new double*[nAtoms_];
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
			dFi_drhoa_[n] = new double[nAtoms_];
		}
	}
	
	double rc, rn, rij, x, omxsq, omxcubed, dphi, astar, Eu, phibar, dEu, dphibar;
	Vec3<double> vij, fpair, fembed;
	static double base_rc = 1.95, base_rn = 1.75;
	
	// Grab parameters from specified potential
	double r0 = eamPotential->parameter(0)->value();
	double E0 = eamPotential->parameter(1)->value();
	double A0 = eamPotential->parameter(2)->value();
	double alpha = eamPotential->parameter(3)->value();
	double beta = eamPotential->parameter(4)->value();
	double rhoscale = eamPotential->parameter(5)->value();
	double Z0 = eamPotential->parameter(6)->value();
	
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
			dFi_drhoa_[i][j] = 0.0;
		}
	}

	// Calculate cutoff function values for each atomic pair
	// [Mei et al, PRb, 43, 4653 (1991), eq 18a/18b]
	rc = base_rc * r0;
	rn = base_rn * r0;

	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			rij = distanceMatrix_[i][j];

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
		}
	}

	// Calculate background electron density for each atom, and atomic densities for each atom pair neighbour distance
	// Jelinek et al, Cond. Matt. arxiv:cond-mat/0610602v4
	//  Eq 4 : Background electronic density around individual atom i
	//		     rho_i^(0)
	//	 rhobar(i) = --------- * G(Gamma_i)
	//		      rho_i^0
	//
	//	rho_i^(k) = SUM  rho_j^a(k)(rij) * Sij			(Eq 9a, for k = 0)
	//		    j<>i
	//
	//	  rho_i^0 = rho_i0 * Z_i0 * G(Gamma_i^ref)		(Eq 7)
	//
	//				[		(  rij	    ) ]
	// rho_i^a(k)(rij) = rho_i0 * exp[ -beta_i^(k) * ( ----- - 1 ) ]	(Eq 10)
	//				[		( r_i^0	    ) ]
	//
	//	    r_i^0 = nearest-neighbour distance in reference structure (EAM parameter 1, r0)
	//      beta_i^(k) = element-dependent parameters (EAM parameter 5, beta)
	//	   rho_i0 = element-dependent density scaling (should be an EAM parameter 6, rhoscale)
	//	     Z_i0 = First nearest-neighbour coordination number of reference system (EAM parameter 7, Z0)
	//	      Sij = cutoff function between atom, calculated above
	//	Neglecting higher-order terms (i.e. only k = 0 considered), Gamma_i = 0 and G(Gamma_i) = 1 (also for Gamma_i^ref)

	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vij = vectorMatrix_[i][j];
			rij = distanceMatrix_[i][j];

			rhoa_[i][j] = rhoscale * exp(-beta*((rij/r0)-1.0));
			rhoa_[j][i] = rhoa_[i][j];

			drhoa_drij_[i][j] = rhoscale * exp( -beta*((rij/r0)-1.0) ) * (-beta/r0);
			drhoa_drij_[j][i] = drhoa_drij_[i][j];

			rho_[i] += rhoa_[i][j] * Sij_[i][j];
			rho_[j] += rhoa_[j][i] * Sij_[j][i];

			drho_drij_[i][j] = drhoa_drij_[i][j]*Sij_[i][j] + rhoa_[i][j]*dSij_drij_[i][j];
			drho_drij_[j][i] = drho_drij_[i][j];
		}
	}

	// If k = 0 only, then Eq 7 reduces to 'rho_i0 * Z_i0'
	for (i=0; i<nAtoms_; ++i)
	{
		rho_[i] /= rhoscale * Z0;
		for (j=0; j<nAtoms_; ++j) drho_drij_[i][j] /= rhoscale * Z0;
	}
// 	rhobar(:) = rhobar(:) / (rhoscale * Z0)
// 	drhobar_drij(:,:) = drhobar_drij(:,:) / (rhoscale * Z0)

	// Determine embedding energy for each atom
	// Jelinek et al, Cond. Matt. arxiv:cond-mat/0610602v4
	//
	//  Eq 3 : Embedding energy per atom
	//
	//  F_i(rhobar(i)) = Ai * Ei * rhobar(i) * ln(rhobar(i))
	//
	//	        Ai = Energy scaling parameter (EAM parameter 3, A0)
	//		Ei = Sublimation energy (EAM parameter 2, E0)
	for (i=0; i<nAtoms_; ++i)
	{
		if (rho_[i] > 1.0e-10)
		{
			Fi_[i] = A0 * E0 * rho_[i] * log(rho_[i]);
			dFi_drho_[i] = (A0*E0)*(1.0 + log(rho_[i]));
		}
		else
		{
			// Note: Derivative at rhobar(i) = 0 should be -INF
		// 	    write(0,"(a,i4,a)") "WARNING - rhobar for atom ",i," is zero...."
			dFi_drho_[i] = -(A0*E0);
		}
	}
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			if (rhoa_[i][j] > 1.0e-10) dFi_drhoa_[i][j] = A0 * E0 * (1.0 + log(rhoa_[i][j])) * drhoa_drij_[i][j];
			else dFi_drhoa_[i][j] = -(A0 * E0);
			dFi_drhoa_[j][i] = dFi_drhoa_[i][j];
		}
	}

	// Calculate pair potential between atoms
	// Jelinek et al, Cond. Matt. arxiv:cond-mat/0610602v4
	//
	// Eq 12 : Pair potentials
	//
	//    phi_ij(rij) = phibar_ij(rij) * Sij
	//
	//		    1  [		      ( Zij		    )	   ( Zij		 ) ]
	// phibar_ij(rij) = --- [ 2 * E_ij^u(rij) - Fi ( --- * rho_j^a0(rij) ) - Fj ( --- * rho_j^a0(rij) ) ]  (Eq 13)
	//		   Zij (		      ( Zi		    )      ( Zj		         ) ]
	//
	//	  E_ij^u = -Eij * (1 + a_ij(rij)) * exp(-a_ij(rij))						  (Eq 14)
	//
	//			     ( rij      )
	//	   a_ij = alpha_ij * ( ---- - 1 )								  (Eq 15)
	//			     ( rij0     )
	//
	//	   E_ij = Sublimation energy (EAM parameter 2 for pure metals, E0)
	//      alpha_ij = Exponential decay factor for universal energy term (EAM parameter 4, alpha)
	//	   Z_ij = "Depends on the structure of the reference system" (EAM parameter 7, Z0)
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			rij = distanceMatrix_[i][j];

			astar = alpha*((rij/r0) - 1.0);

			Eu = -E0 * (1.0 + astar) * exp(-astar);

			x = A0 * E0 * rhoa_[i][j] * log(rhoa_[i][j]);	// * (Zij / Zi)
			phibar = (2.0 * Eu - x - x) / Z0;
			
			phi_[i][j] = phibar * Sij_[i][j];
			phi_[j][i] = phi_[i][j];

			// Derivatives
			dEu = E0 * astar * exp(-astar) * (alpha / r0);
			x = dFi_drhoa_[i][j];	// * (Zij / Zi);
			dphibar = (2.0 * dEu - x - x) / Z0;

			dphi_drij_[i][j] = dphibar * Sij_[i][j] + phibar * dSij_drij_[i][j];
			dphi_drij_[j][i] = dphi_drij_[i][j];

		}
	}

	// Calculate forces on atoms
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vij = vectorMatrix_[i][j];
			rij = distanceMatrix_[i][j];
			vij /= rij;
			
			// From pair potential
// 			fvec(1:3) = -dphi_drij(i,j) * rij(i,j,1:3) / rij(i,j,4)
			fpair = vij * -dphi_drij_[i][j];
			fContributions_[Configuration::VdwContribution][i] -= fpair;
			fContributions_[Configuration::VdwContribution][j] += fpair;

			// From embedding function
// 			fvec(1:3) = -(dFi_drhobar(i) * drhobar_drij(i,j) + dFi_drhobar(j) * drhobar_drij(i,j)) * rij(i,j,1:3) / rij(i,j,4)
			fembed = vij * -(dFi_drho_[i] * drho_drij_[i][j] + dFi_drho_[j] * drho_drij_[i][j]);
			fContributions_[Configuration::VdwContribution][i] -= fembed;
			fContributions_[Configuration::VdwContribution][j] += fembed;
		}
	}

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

// Delete EAM-related arrays
void Configuration::deleteEAMArrays()
{
	if (Sij_ != NULL) for (int n=0; n<nAtoms_; ++n)
	{
		if (Sij_[n] != NULL) delete[] Sij_[n];
		if (dSij_drij_[n] != NULL) delete[] dSij_drij_[n];
		if (drho_drij_[n] != NULL) delete[] drho_drij_[n];
		if (rhoa_[n] != NULL) delete[] rhoa_[n];
		if (drhoa_drij_[n] != NULL) delete[] drhoa_drij_[n];
		if (phi_[n] != NULL) delete[] phi_[n];
		if (dphi_drij_[n] != NULL) delete[] dphi_drij_[n];
		if (dFi_drhoa_[n] != NULL) delete[] dFi_drhoa_[n];
	}
	if (Sij_ != NULL) delete[] Sij_;
	if (dSij_drij_ != NULL) delete[] dSij_drij_; 
	if (rho_ != NULL) delete[] rho_;
	if (drho_drij_ != NULL) delete[] drho_drij_; 
	if (rhoa_ != NULL) delete[] rhoa_;
	if (drhoa_drij_ != NULL) delete[] drhoa_drij_;
	if (Fi_ != NULL) delete[] Fi_;
	if (dFi_drho_ != NULL) delete[] dFi_drho_;
	if (dFi_drhoa_ != NULL) delete[] dFi_drhoa_;
	if (phi_ != NULL) delete[] phi_;
	if (dphi_drij_ != NULL) delete[] dphi_drij_;
}
