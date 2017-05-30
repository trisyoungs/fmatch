/*
	*** Ewald sum energy / force calculation
	*** ewald.cpp
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

// Ref: "Long-Range Interactions in Many-Particle Simulation", P. Gibbon and G. Sutmann
//	In "Quantum Simulations of Complex Many-Body Systems; From Theory to Algorithms"
//	NIC Series, Vol. 10, ISBN 3-00-009057, pp. 467-506, 2002.

#include <math.h>
#include "system.h"
#include "vector3.h"
#include "atompair.h"
#include "configuration.h"
#include "species.h"
#include "kvector.h"

#define SQRTPI 1.772453850905515
#define TWOPI 6.283185307179586

// Error Function
double TargetSystem::erfc(double x)
{
	// Approximation to the complementary error function.
	// Ref: Abramowitz and Stegun, Handbook of Mathematical Functions,
	//      National Bureau of Standards, Formula 7.1.26
	static double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
	double t, tp, result;
	t = 1.0 / ( 1.0 + p * x );
	tp = t * ( a1 + t * ( a2 + t * ( a3 + t * ( a4 + t * a5 ) ) ) );
	result = tp * exp(-(x*x));
	return result;
}

// Complementary error function
double TargetSystem::erf(double x)
{
	return (1.0 - erfc(x));
}

// Ewald Sum Real-space contributions

// Energy:
//				       N   N		    erfc(alpha * |rij + n|) 
// 		   E(real) = 0.5 * E'  E   E  q(i) * q(j) * -----------------------
// (p 474)			   n  i=1 j=1			   |rij + n|
// 'n' is box vector - here we only consider the minimimum image coordinates of the atoms in the central box (n=0)
// Factor of 1/2 is not required in the summation since the sums go from i=0,N-1 and j=i,N

// Forces:
//			    N-1  N  q(i) * q(j)				2*alpha*rij			       -->
//		F(real) = E' E   E  ----------- * ( erfc(alpha * rij) + ----------- * exp(-(alpha*rij)**2) ) * rij
//			  n i=1 j>i   rij**3				   sqrtpi


// Calculate Ewald sum energy and forces
void TargetSystem::ewaldReal(Configuration *cfg)
{
	// Calculate real-space contributions to the Ewald sum.
	int i,j,aoff,m1,con, n;
	Vec3<double> fi, vij;
	double rij, energy, alpharij, factor, qiqj, qqrij3, qi;

	energy = 0.0;

	// Real-space contribution to energy and forces
	if (useAtomPairs_)
	{
		// Use the stored pairs list...
		AtomPair *pair = cfg->atomPairs();
		for (n=0; n<cfg->nAtomPairs(); ++n)
		{
	// 		printf("AtomPairEwald = %i\n", n);
			// Grab some numbers
			rij = pair->rij();
			qiqj = pair->chargeI() * pair->chargeJ() * pair->elecScale();
			alpharij = ewaldAlpha_ * rij;
			
			// Energy
			energy += qiqj * erfc(alpharij) / rij;

			// Forces
			factor = erfc(alpharij) + 2.0*alpharij/SQRTPI * exp(-(alpharij*alpharij));
			qqrij3 = factor * qiqj / (rij * rij * rij);
			qqrij3 *= elecConvert_;
			fi = pair->vij() * qqrij3;
			cfg->subtractForce(Configuration::ChargeContribution, fi, pair->i());
			cfg->addForce(Configuration::ChargeContribution, fi, pair->j());
			
			++pair;
		}
	}
	else
	{
		for (i=0; i<nAtoms_-1; ++i)
		{
			qi = forcefield_.charge(typeMap_[i]);
			for (j=i+1; j<nAtoms_; ++j)
			{
				rij = cfg->distance(i,j);
				if (rij > cutoff_) continue;
				vij = cfg->vector(i,j);
				
				qiqj = qi * forcefield_.charge(typeMap_[j]) * elecScaleMatrix_[i][j];
				alpharij = ewaldAlpha_ * rij;
				
				// Energy
				energy += qiqj * erfc(alpharij) / rij;

				// Forces
				factor = erfc(alpharij) + 2.0*alpharij/SQRTPI * exp(-(alpharij*alpharij));
				qqrij3 = factor * qiqj / (rij * rij * rij);
				qqrij3 *= elecConvert_;
				fi = vij * qqrij3;
				cfg->subtractForce(Configuration::ChargeContribution, fi, i);
				cfg->addForce(Configuration::ChargeContribution, fi, j);
			}
		}
	}

	// Sum energy
// 	printf("Real = %f\n", energy * elecConvert_);
	cfg->addEnergy(Configuration::ChargeContribution, energy * elecConvert_);
}

// Ewald Reciprocal Space Calculations:
// Energy:
//			    2*pi    N   N					    -ksq	 1
//		E(recip) =  ---- E' E   E  q(i) * q(j) * exp(ik.(rj - ri)) * exp( --------- ) * ---
//			    L**3 k i=1 j=1					  4*alphasq	ksq
//
// Forces:
//				 N
//		  F(recip) = E   E q(j) 
//			    k/=0 j

void TargetSystem::ewaldReciprocal(Configuration *cfg)
{
	// Calculate the reciprocal contribution of all atoms to the Ewald sum.
	// Note that all exponential factors are precalculated and contained within the KVector list
	int kx, ky, kz, i, n, m, k;
	double xycos, xysin, xyzcos, xyzsin, energy, sumcos, sumsin;
	static double *atomcos = NULL, *atomsin = NULL;
	static int lastNAtoms = 0;
	Vec3<double> **rCos, **rSin, fi, rvec;

	// Get pointers to cos and sin terms for atoms
	int kmax = ewaldKMax_.max();
	// Create/resize temporary arrays
	if (lastNAtoms < cfg->nAtoms() || atomcos == NULL)
	{
		if (atomcos == NULL) delete[] atomcos;
		if (atomsin == NULL) delete[] atomsin;
		lastNAtoms = cfg->nAtoms();
		atomcos = new double[lastNAtoms];
		atomsin = new double[lastNAtoms];
	}
	//printf("Cutoffsq = %f  (%f)\n",cutoffsq,sqrt(cutoffsq));

	energy = 0.0;
	KVector *kvec = cfg->kVectors();
	rCos = cfg->rCos();
	rSin = cfg->rSin();
	for (k=0; k<cfg->nKVectors(); ++k)
	{
		// Now sum contributions over atoms
		sumcos = 0.0;
		sumsin = 0.0;
		i = 0;
		kx = kvec->x();
		ky = kvec->y();
		kz = kvec->z();
		for (Species *sp = species_.first(); sp != NULL; sp = sp->next)
		{
			for (m=0; m<sp->nMolecules(); ++m)
			{
				for (n=0; n<sp->nAtoms(); ++n)
				{
					// Calculate k-vector (x*y)
					xycos = rCos[abs(kx)][i].x * rCos[abs(ky)][i].y - rSin[kmax+kx][i].x * rSin[kmax+ky][i].y;
					xysin = rCos[abs(kx)][i].x * rSin[kmax+ky][i].y + rSin[kmax+kx][i].x * rCos[abs(ky)][i].y;
					// Calculate k-vector (xy*z);
					xyzcos = xycos * rCos[abs(kz)][i].z - xysin * rSin[kmax+kz][i].z;
					xyzsin = xycos * rSin[kmax+kz][i].z + xysin * rCos[abs(kz)][i].z;
					atomcos[i] = sp->atom(n)->charge() * xyzcos;
					atomsin[i] = sp->atom(n)->charge() * xyzsin;
					sumcos += atomcos[i];
					sumsin += atomsin[i];
					++i;
				}
			}
		}
// 		printf("k/sumcos/simsin = %i %f %f\n", k, sumcos, sumsin);
		// Calculate energy contributions
		energy += (sumcos*sumcos + sumsin*sumsin) * kvec->factor();
		
		// Calculate forces
		rvec = kvec->vector() * 2.0 * kvec->factor();
		for (i=0; i<nAtoms_; ++i)
		{
// 			fi = kvec->vector() * (2.0 * exp1 * (atomsin[i]*sumcos - atomcos[i]*sumsin)) * kvec->factor();
			fi = rvec * (atomsin[i]*sumcos - atomcos[i]*sumsin);
			cfg->addForce(Configuration::ChargeContribution, fi, i);
		}
		++kvec;
	}
	
	// Sum energy
// 	printf("Recip = %f\n", energy);
	cfg->addEnergy(Configuration::ChargeContribution, energy);
}

void TargetSystem::quickEwaldReciprocal(Configuration *cfg)
{
	// Calculate the reciprocal contribution of all atoms to the Ewald sum, using precalculated sums in configuration
	int kx, ky, kz, i, n, m, k;
	double energy;
	double sumcos, sumsin, **atomCos, **atomSin, **typeCos, **typeSin;
	int *chargeMap;
	Vec3<double> fi, rvec;

	// Get pointers to cos and sin terms for atoms and types
	energy = 0.0;
	KVector *kvec = cfg->kVectors();
	atomCos = cfg->atomCos();
	atomSin = cfg->atomSin();
	typeCos = cfg->typeCos();
	typeSin = cfg->typeSin();
	chargeMap = cfg->chargeMap();

	// Loop over kvectors
	for (k=0; k<cfg->nKVectors(); ++k)
	{
		// Calculate total sums from sums partitioned into charge types
		sumcos = 0.0;
		sumsin = 0.0;
		for (n=0; n<forcefield_.nCharges(); ++n)
		{
			sumcos += typeCos[k][n] * forcefield_.charge(n);
			sumsin += typeSin[k][n] * forcefield_.charge(n);
		}
// 		printf("k/sumcos/simsin = %i %f %f\n", k, sumcos, sumsin);
		// Calculate energy contributions from the interactions of patterns
		energy += (sumcos*sumcos + sumsin*sumsin) * kvec->factor();

		// Calculate forces
		rvec = kvec->vector() * 2.0 * kvec->factor();
		for (i=0; i<nAtoms_; ++i)
		{
// 			fi = kvec->vector() * (2.0 * exp1 * (atomsin[i]*sumcos - atomcos[i]*sumsin)) * kvec->factor();
			fi = rvec * (atomSin[k][i]*sumcos - atomCos[k][i]*sumsin) * forcefield_.charge(chargeMap[i]);
			cfg->addForce(Configuration::ChargeContribution, fi, i);
		}
		++kvec;
	}
	
	// Sum energy
// 	printf("Recip = %f\n", energy);
	cfg->addEnergy(Configuration::ChargeContribution, energy);
}

// Ewald Corrections
//			    alpha   N
//		E(self) = - ------  E  q(i) * q(i) 
//			    sqrtpi i=1
//
//				    * *		      erf(alpha(rij))
//		E(molecular) = - E  E E q(i) * q(j) * ---------------
//				 m  i j			    rij
// Sums over i=* and j=* indicate excluded interactions, i.e. bond i-j, angle i-x-j and torsion i-x-x-j.

void TargetSystem::ewaldCorrect(Configuration* cfg)
{
	int i, j, n, m;
	double molcorrect, energy, qi, qiqj, rij, chargesum, alpharij, factor, qqrij3;
	Vec3<double> fi, vij;
	
	// Self-interaction correction
	chargesum = 0.0;
	for (Species *sp = species_.first(); sp != NULL; sp = sp->next)
	{
		for (m=0; m<sp->nMolecules(); ++m)
		{
			for (n=0; n<sp->nAtoms(); ++n) chargesum += sp->atom(n)->charge() * sp->atom(n)->charge();
		}
	}
// 	printf("Ewald Self Energy = %f\n", -(ewaldAlpha_/SQRTPI) * chargesum * elecConvert_);
	cfg->addEnergy(Configuration::ChargeContribution, -(ewaldAlpha_/SQRTPI) * chargesum * elecConvert_);
// 	energy = (alpha/SQRTPI) * chargesum * prefs.elecConvert();
// 	estore->add(EnergyStore::EwaldSelfEnergy,energy,id_);

	// Correct the reciprocal Ewald sum for intramolecular interactions
	energy = 0.0;
	if (useAtomPairs_)
	{
		AtomPair *pair = cfg->excludedPairs();
		for (n=0; n<cfg->nExcludedPairs(); ++n)
		{
			// Energy correction
			rij = pair->rij();
			qiqj = pair->chargeI() * pair->chargeJ() * (1.0 - pair->elecScale());
			energy += qiqj * erf(ewaldAlpha_*rij)/rij;
			
			// Force correction
			alpharij = ewaldAlpha_ * rij;
			factor = erf(alpharij) - 2.0*alpharij/SQRTPI * exp(-(alpharij*alpharij));
			factor *= qiqj * elecConvert_ / (rij * rij * rij);
			fi = pair->vij() * factor;
			cfg->addForce(Configuration::ChargeContribution, fi, pair->i());
			cfg->subtractForce(Configuration::ChargeContribution, fi, pair->j());
			++pair;
		}
	}
	else
	{
		for (i=0; i<nAtoms_-1; ++i)
		{
			qi = forcefield_.charge(typeMap_[i]);
			for (j=i+1; j<nAtoms_; ++j)
			{
				if (elecScaleMatrix_[i][j] > 0.999) continue;
				rij = cfg->distance(i,j);
				if (rij > cutoff_) continue;
				vij = cfg->vector(i,j);
				
				qiqj = qi * forcefield_.charge(typeMap_[j]) * (1.0 - elecScaleMatrix_[i][j]);

				// Energy
				energy += qiqj * erf(ewaldAlpha_*rij)/rij;

				// Forces
				alpharij = ewaldAlpha_ * rij;
				factor = erf(alpharij) - 2.0*alpharij/SQRTPI * exp(-(alpharij*alpharij));
				factor *= qiqj * elecConvert_ / (rij * rij * rij);
				fi = vij * factor;
				cfg->addForce(Configuration::ChargeContribution, fi, i);
				cfg->subtractForce(Configuration::ChargeContribution, fi, j);
			}
		}
	}

// 	printf("Intermol correction energy = %f\n", -energy * elecConvert_);
	cfg->addEnergy(Configuration::ChargeContribution, -energy * elecConvert_);
// 	estore->add(EnergyStore::EwaldMolecularEnergy,energy,id_);
}
