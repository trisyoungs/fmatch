/*
	*** Atomic Configuration
	*** configuration.cpp
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
#include "species.h"
#include "ff.h"
#include "atompair.h"
#include "kvector.h"
#include "system.h"
#include <math.h>

// Constructor
Configuration::Configuration()
{
	// Private variables
	r_ = NULL;
	f_ = NULL;
	fWeight_ = NULL;
	fContributions_ = NULL;
	fCalc_ = NULL;
	cellType_ = 0;
	for (int n=0; n<Configuration::nContributionTypes; ++n)
	{
		energies_[n] = 0.0;
		logPoint_[n] = 0;
		fitEnergyValues_[n] = 0.0;
		fitEnergyType_[n] = FALSE;
		fitEnergyWeights_[n] = 1.0;
	}
	distanceMatrix_ = NULL;
	vectorMatrix_ = NULL;
	atomPairs_ = NULL;
	nAtomPairs_ = 0;
	excludedPairs_ = NULL;
	nExcludedPairs_ = 0;
	nKVectors_ = 0;
	kVectors_ = NULL;
	rCos_ = NULL;
	rSin_ = NULL;
	atomCos_ = NULL;
	atomSin_ = NULL;
	typeCos_ = NULL;
	typeSin_ = NULL;
	chargeMap_ = NULL;
	Sij_ = NULL;

	// Public variables
	next = NULL;
	prev = NULL;
}

// Destructor
Configuration::~Configuration()
{
	deleteAtomArrays();
	deletePrecalculatedArrays();
}

/*
// Basic Data
*/

// Set name of configuration
void Configuration::setName(const char *name)
{
	name_ = name;
}

// Return name of this configuration
const char *Configuration::name()
{
	return name_.get();
}

/*
// Atomic Data
*/

// Set atom data
bool Configuration::setAtomData(int id, Vec3<double> coordinates, Vec3<double> forces, double weight)
{
	if (id < 0 || id >= nAtoms_)
	{
		printf("Error: Atom index in configuration is out of range.\n");
		return FALSE;
	}
	r_[id] = coordinates;
	f_[id] = forces;
	fWeight_[id] = weight;
	return TRUE;
}

// Set atom position data
bool Configuration::setAtomCoordinates(int id, Vec3<double> coordinates)
{
	if (id < 0 || id >= nAtoms_)
	{
		printf("Error: Atom index in configuration is out of range.\n");
		return FALSE;
	}
	r_[id] = coordinates;
	return TRUE;
}

// Set atom forces data
bool Configuration::setAtomForces(int id, Vec3<double> forces, double weight)
{
	if (id < 0 || id >= nAtoms_)
	{
		printf("Error: Atom index in configuration is out of range.\n");
		return FALSE;
	}
	f_[id] = forces;
	fWeight_[id] = weight;
	return TRUE;
}

// Return number of atoms in configuration
int Configuration::nAtoms()
{
	return nAtoms_;
}

// Return coordinates of atom specified
Vec3<double> Configuration::r(int id)
{
	return r_[id];
}

// Return array of atomic positions
Vec3<double> *Configuration::r()
{
	return r_;
}

/*
// Unit Cell
*/

// Setup unit cell
void Configuration::setCell(Matrix &cell, int cellType)
{
	// Set cell information
	cellType_ = cellType;
	cell_ = cell;
	inverseCell_ = cell;
	inverseCell_.invert();
	reciprocalCell_.setColumn(0, cell_.columnAsVec3(1) * cell_.columnAsVec3(2), 0.0);
	reciprocalCell_.setColumn(1, cell_.columnAsVec3(0) * cell_.columnAsVec3(2), 0.0);
	reciprocalCell_.setColumn(2, cell_.columnAsVec3(0) * cell_.columnAsVec3(1), 0.0);
	reciprocalVolume_ = 1.0 / fabs( cell_[0]*reciprocalCell_[0] + cell_[5]*reciprocalCell_[5] + cell_[10]*reciprocalCell_[10]);
	reciprocalCell_.columnMultiply(0, reciprocalVolume_);
	reciprocalCell_.columnMultiply(1, reciprocalVolume_);
	reciprocalCell_.columnMultiply(2, reciprocalVolume_);
}

// Return minimum image vector between supplied atoms
Vec3<double> Configuration::mimd(int i, int j)
{
	Vec3<double> result, half;
	switch (cellType_)
	{
		// No cell
		case (0):
			result = r_[j]-r_[i];
			break;
		// Cubic / orthorhombic cell
		case (1):
		case (2):
			result = r_[j]-r_[i];
			result.x -= floor(result.x/cell_[0] + 0.5)*cell_[0];
			result.y -= floor(result.y/cell_[5] + 0.5)*cell_[5];
			result.z -= floor(result.z/cell_[10] + 0.5)*cell_[10];
			break;
		// Parallelepiped cell
		case (3):
			result = inverseCell_.transform(r_[j]-r_[i]);
			result.x -= floor(result.x + 0.5);
			result.y -= floor(result.y + 0.5);
			result.z -= floor(result.z + 0.5);
			result = cell_.transform(result);
			break;
	}
	return result;
}

// Return reciprocal cell of configuration
Matrix &Configuration::reciprocalCell()
{
	return reciprocalCell_;
}

// Return reciprocal cell volume
double Configuration::reciprocalVolume()
{
	return reciprocalVolume_;
}

// Return reciprocal cell cutoff for Ewald sum
double Configuration::reciprocalCutoff()
{
	return reciprocalCutoff_;
}

/*
// Calculated Energy / Forces
*/

// Increase specified energy term
void Configuration::addEnergy(ContributionType type, double e)
{
	energies_[type] += e;
}

// Return specified energy term
double Configuration::energy(ContributionType type)
{
	return energies_[type];
}

// Return total energy
double Configuration::totalEnergy()
{
	double total = 0.0;
	for (int n=0; n<Configuration::nContributionTypes; ++n) total += energies_[n];
	return total;
}


// Add calculated force to specified atom
void Configuration::addForce(Configuration::ContributionType type, Vec3<double> &fvec, int i)
{
	if (i < 0 || i >= nAtoms_) printf("BUG: Index is out of bounds in addForce.\n");
	else fContributions_[type][i] += fvec;
}

// Subtract calculated force from specified atom
void Configuration::subtractForce(Configuration::ContributionType type, Vec3<double> &fvec, int i)
{
	if (i < 0 || i >= nAtoms_) printf("BUG: Index is out of bounds in subtractForce.\n");
	else fContributions_[type][i] -= fvec;
}

// Calculate total forces from individual contributions
void Configuration::sumForceContributions()
{
	int type;
	for (int n=0; n<nAtoms_; ++n)
	{
		fCalc_[n] = 0.0;
		for (type = 0; type<Configuration::nContributionTypes; ++type) fCalc_[n] += fContributions_[type][n];
	}
}

// Calculate and return sum of squares error in forces and any other specified quantities
double Configuration::totalSOSE()
{
	double sose = 0.0;
	double d;
	int n;

	// Atomic forces
	Vec3<double> fdelta;
	for (n=0; n<nAtoms_; ++n)
	{
		fdelta = f_[n] - fCalc_[n];
		sose += (fdelta.x*fdelta.x + fdelta.y*fdelta.y + fdelta.z*fdelta.z) * fWeight_[n];
	}

	// Energy components
	for (n = 0; n<Configuration::nContributionTypes; ++n)
	{
		if (!fitEnergyType_[n]) continue;
		d = fitEnergyValues_[n] - energies_[n];
		d *= d;
		sose += d * fitEnergyWeights_[n];
	}
	
	// ESP Points
	for (ESPPoint *esp = espPoints_.first(); esp != NULL; esp = esp->next)
	{
		// TODO
	}
	
	return sose;
}

// Return calculated force array
Vec3<double> *Configuration::fCalc()
{
	return fCalc_;
}

// Return specified logpoint
unsigned int Configuration::logPoint(ContributionType type)
{
	return logPoint_[type];
}

// Update specified logpoint
void Configuration::updateLogPoint(Forcefield *ff, ContributionType type)
{
	logPoint_[type] = ff->logPoint(type);
}

/*
// Precalculated Quantities
*/

// Update/create distance matrix
void Configuration::updateDistanceMatrix()
{
	if (distanceMatrix_ == NULL)
	{
		distanceMatrix_ = new double*[nAtoms_];
		vectorMatrix_ = new Vec3<double>*[nAtoms_];
		for (int i=0; i<nAtoms_; ++i)
		{
			distanceMatrix_[i] = new double[nAtoms_];
			vectorMatrix_[i] = new Vec3<double>[nAtoms_];
		}
	}
	// Recalculate mim vectors for all atoms
	int i, j;
	for (i=0; i<nAtoms_-1; ++i)
	{
		for (j=i+1; j<nAtoms_; ++j)
		{
			vectorMatrix_[i][j] = mimd(i, j);
			distanceMatrix_[i][j] = vectorMatrix_[i][j].magnitude();
			vectorMatrix_[j][i] = -vectorMatrix_[i][j];
			distanceMatrix_[j][i] = distanceMatrix_[i][j];
		}
	}
}

// Return distance for supplied pair
double Configuration::distance(int i, int j)
{
	return distanceMatrix_[i][j];
}

// Return mim vector for supplied pair
Vec3<double> Configuration::vector(int i, int j)
{
	return vectorMatrix_[i][j];
}

// (Re)create atompairs arrays
void Configuration::createPairArrays(int natompairs, int nexcludedpairs)
{
	nAtomPairs_ = natompairs;
	nExcludedPairs_ = nexcludedpairs;
	if (atomPairs_ != NULL) delete[] atomPairs_;
	atomPairs_ = new AtomPair[nAtomPairs_];
	if (excludedPairs_ != NULL) delete[] excludedPairs_;
	excludedPairs_ = new AtomPair[nExcludedPairs_];
}

// Return number of interacting atom pairs
int Configuration::nAtomPairs()
{
	return nAtomPairs_;
}

// Return interacting pair list
AtomPair *Configuration::atomPairs()
{
	return atomPairs_;
}

// Return number of excluded atom pairs
int Configuration::nExcludedPairs()
{
	return nExcludedPairs_;
}

// Return excluded atom pair list
AtomPair *Configuration::excludedPairs()
{
	return excludedPairs_;
}

// Create k-vector arrays
void Configuration::createKVectorArrays(Vec3<int> kmax, double alpha, double elecConvert)
{
	int i, n;
	if (rCos_ != NULL)
	{
		for (i=0; i<kMax_.max()+1; ++i) delete[] rCos_[i];
		for (i=0; i<2*kMax_.max()+1; ++i) delete[] rSin_[i];
		delete[] rCos_;
		delete[] rSin_;
	}
	kMax_ = kmax;
	rCos_ = new Vec3<double>*[kMax_.max()+1];
	for (i=0; i<kMax_.max()+1; i++) rCos_[i] = new Vec3<double>[nAtoms_];
	rSin_ = new Vec3<double>*[2*kMax_.max()+1];
	for (i=0; i<2*kMax_.max()+1; i++) rSin_[i] = new Vec3<double>[nAtoms_];

	// Calculate list of k-vectors within cutoff
	Matrix rcell = reciprocalCell_ * TWOPI;
	Vec3<double> rvec;
	double magsq, factor, cutoffsq = reciprocalCutoff_ * reciprocalCutoff_, alphasq = alpha*alpha;
	int kx, ky, kz;
	for (int pass=0; pass < 2; ++pass)
	{
		nKVectors_ = 0;
		for (kx=0; kx<=kMax_.x; ++kx)
			for (ky=-kMax_.y; ky<=kMax_.y; ++ky)
				for (kz=-kMax_.z; kz<=kMax_.z; ++kz)
				{
					if ((kx == 0) && (ky == 0) && (kz == 0)) continue;
					rvec.set(kx,ky,kz);
					rvec = rcell * rvec;
					magsq = rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z;
					if (magsq > cutoffsq) continue;
					factor = reciprocalVolume_ * TWOPI * elecConvert;
					factor *= exp(-magsq/(4.0*alphasq))/magsq * (kx == 0 ? 1.0 : 2.0);
					// Add to list
					if (pass == 1) kVectors_[nKVectors_].set(kx, ky, kz, rvec, magsq, factor);
					++nKVectors_;
				}
		// End of first pass? If so, create array...
		if (pass == 0)
		{
			if (kVectors_ != NULL) delete[] kVectors_;
			kVectors_ = new KVector[nKVectors_];
		}
	}
	
	// Create initial arrays
	updateKVectorArrays();
}

// Update k-vector arrays
void Configuration::updateKVectorArrays()
{
	// Check stored kmax value
	if (kMax_.max() < 1)
	{
		printf("Error: KVector arrays not created yet?\n");
		return;
	}
	
	// Generate reciprocal space coordinates for atoms at each cartesian k-vector.
	// Only create positive k-vector positions for cos since is an even function. For sin calculate negative
	// also (odd function), where rSin[k] runs from k=0,2*kMax_+1, with k=kMax_ the central (zero) vector (=complex(1,0)).
	// Thus, rSin[kMax_+i] gives the i'th positive kvector and rsun[kMax_-i] gives the i'th negative kvector.
	// Grab reciprocal cell and multiply by TWOPI
	Matrix rcell = reciprocalCell_ * TWOPI;
	Vec3<double> pos;
	int i, n, kmax = kMax_.max();
	for (i=0; i<nAtoms_; ++i)
	{
		// Set central (zero) vector elements to be cmplx(1.0,0.0)
		rCos_[0][i].set(1.0,1.0,1.0);
		rSin_[0][i].set(0.0,0.0,0.0);
		// Calculate first vector in the positive k-direction
		pos.x = rcell.columnAsVec3(0).dp(r_[i]);
		pos.y = rcell.columnAsVec3(1).dp(r_[i]);
		pos.z = rcell.columnAsVec3(2).dp(r_[i]);
		rCos_[1][i].x = cos(pos.x);
		rCos_[1][i].y = cos(pos.y);
		rCos_[1][i].z = cos(pos.z);
		rSin_[kmax+1][i].x = sin(pos.x);
		rSin_[kmax+1][i].y = sin(pos.y);
		rSin_[kmax+1][i].z = sin(pos.z);
		// Calculate vector in the negative k-direction for sin terms
		rSin_[kmax-1][i].x = -rSin_[kmax+1][i].x;
		rSin_[kmax-1][i].y = -rSin_[kmax+1][i].y;
		rSin_[kmax-1][i].z = -rSin_[kmax+1][i].z;
		// Calculate the extended reciprocal space position vectors (power expansion of first vectors).
		// TODO Timewaste: Assume, for now, that we have a cubic box so we do all three vector arrays at once, and up to kmax. 
		// Build up the kvector positions by multiplying by rCos_[+-1]/rSin_[+-1] each time
		int firstsin = kmax+1;
		for (n=1; n<kmax; n++)
		{
			// Complex multiplication: (a + bi)(c + di) = (ac - bd) + (ad + bc)i
			// Positive k-direction (cos and sin)
			int k = n+1;
			int sinpos = kmax + n;
			rCos_[k][i].x = rCos_[1][i].x * rCos_[n][i].x - rSin_[firstsin][i].x * rSin_[sinpos][i].x;
			rCos_[k][i].y = rCos_[1][i].y * rCos_[n][i].y - rSin_[firstsin][i].y * rSin_[sinpos][i].y;
			rCos_[k][i].z = rCos_[1][i].z * rCos_[n][i].z - rSin_[firstsin][i].z * rSin_[sinpos][i].z;
			rSin_[kmax+k][i].x = rCos_[1][i].x * rSin_[sinpos][i].x + rSin_[firstsin][i].x * rCos_[n][i].x;
			rSin_[kmax+k][i].y = rCos_[1][i].y * rSin_[sinpos][i].y + rSin_[firstsin][i].y * rCos_[n][i].y;
			rSin_[kmax+k][i].z = rCos_[1][i].z * rSin_[sinpos][i].z + rSin_[firstsin][i].z * rCos_[n][i].z;
			// Negative k-direction (sin)
			rSin_[kmax-k][i].x = -rSin_[kmax+k][i].x;
			rSin_[kmax-k][i].y = -rSin_[kmax+k][i].y;
			rSin_[kmax-k][i].z = -rSin_[kmax+k][i].z;
		}
	}
}

// Create static kvector/position sums for atom types
void Configuration::createStaticKVectorSums(Species *firstSp, int kMax, int nChargeTypes)
{
	// Create arrays first
	int k;
	if (atomCos_ != NULL)
	{
		for (k=0; k<nKVectors_; ++k)
		{
			delete[] atomCos_[k];
			delete[] atomSin_[k];
			delete[] typeCos_[k];
			delete[] typeSin_[k];
		}
		delete[] atomCos_;
		delete[] atomSin_;
		delete[] typeCos_;
		delete[] typeSin_;
		delete[] chargeMap_;
	}
	atomCos_ = new double*[nKVectors_];
	atomSin_ = new double*[nKVectors_];
	typeCos_ = new double*[nKVectors_];
	typeSin_ = new double*[nKVectors_];
	chargeMap_ = new int[nAtoms_];
	for (k=0; k<nKVectors_; ++k)
	{
		atomCos_[k] = new double[nAtoms_];
		atomSin_[k] = new double[nAtoms_];
		typeCos_[k] = new double[nChargeTypes];
		typeSin_[k] = new double[nChargeTypes];
	}
// 	printf("NVECTORS = %i, size (Mb) = %f\n", nKVectors_, nKVectors_*2*(nAtoms_+nChargeTypes)*8/1024.0/1024.0);

	int ct, n, m, i, kx, ky, kz;
	double xycos, xysin, xyzcos, xyzsin;
	for (k=0; k<nKVectors_; ++k)
	{
		// Zero arrays
		for (ct = 0; ct < nChargeTypes; ++ct)
		{
			typeCos_[k][ct] = 0.0;
			typeSin_[k][ct] = 0.0;
		}

		// Now sum contributions over atoms
		i = 0;
		kx = kVectors_[k].x();
		ky = kVectors_[k].y();
		kz = kVectors_[k].z();
		for (Species *sp = firstSp; sp != NULL; sp = sp->next)
		{
			for (m=0; m<sp->nMolecules(); ++m)
			{
				for (n=0; n<sp->nAtoms(); ++n)
				{
					// Calculate k-vector (x*y)
					xycos = rCos_[abs(kx)][i].x * rCos_[abs(ky)][i].y - rSin_[kMax+kx][i].x * rSin_[kMax+ky][i].y;
					xysin = rCos_[abs(kx)][i].x * rSin_[kMax+ky][i].y + rSin_[kMax+kx][i].x * rCos_[abs(ky)][i].y;
					// Calculate k-vector (xy*z);
					xyzcos = xycos * rCos_[abs(kz)][i].z - xysin * rSin_[kMax+kz][i].z;
					xyzsin = xycos * rSin_[kMax+kz][i].z + xysin * rCos_[abs(kz)][i].z;
					atomCos_[k][i] = xyzcos;
					atomSin_[k][i] = xyzsin;
					chargeMap_[i] = sp->atom(n)->chargeParameter()->parent()->id();
					if (chargeMap_[i] < 0) printf("Hideous error! ID in charge parameter interaction has not been set.\n");
					typeCos_[k][chargeMap_[i]] += atomCos_[k][i];
					typeSin_[k][chargeMap_[i]] += atomSin_[k][i];
					++i;
				}
			}
		}
	}
}

// Return rCos array
Vec3<double> **Configuration::rCos()
{
	return rCos_;
}

// Return rSin array
Vec3<double> **Configuration::rSin()
{
	return rSin_;
}

// Return atomCos array
double **Configuration::atomCos()
{
	return atomCos_;
}

// Return atomSin array
double **Configuration::atomSin()
{
	return atomSin_;
}

// Return typeCos array
double **Configuration::typeCos()
{
	return typeCos_;
}

// Return typeSin array
double **Configuration::typeSin()
{
	return typeSin_;
}

// Return atom charge map
int *Configuration::chargeMap()
{
	return chargeMap_;
}

// Return kVector array
KVector *Configuration::kVectors()
{
	return kVectors_;
}

// Return number of kVectors in array
int Configuration::nKVectors()
{
	return nKVectors_;
}

// Return whether configuration has precalculated atom/type arrays for Ewald sum
bool Configuration::hasStaticKVectorSums()
{
	return (atomCos_ != NULL);
}

// Free precalculated data arrays
void Configuration::deletePrecalculatedArrays()
{
	if (atomPairs_ != NULL) delete[] atomPairs_;
	atomPairs_ = NULL;
	if (excludedPairs_ != NULL) delete[] excludedPairs_;
	excludedPairs_ = NULL;
	if (kVectors_ != NULL) delete[] kVectors_;
	kVectors_ = NULL;
	if (rCos_ != NULL)
	{
		for (int i=0; i<kMax_.max()+1; ++i) delete[] rCos_[i];
		for (int i=0; i<2*kMax_.max()+1; ++i) delete[] rSin_[i];
		delete[] rCos_;
		delete[] rSin_;
		rCos_ = NULL;
		rSin_ = NULL;
	}
	if (atomCos_ != NULL)
	{
		for (int i=0; i<nKVectors_; ++i)
		{
			delete[] atomCos_[i];
			delete[] atomSin_[i];
			delete[] typeCos_[i];
			delete[] typeSin_[i];
		}
		delete[] atomCos_;
		delete[] atomSin_;
		delete[] typeCos_;
		delete[] typeSin_;
		atomCos_ = NULL;
		atomSin_ = NULL;
		typeCos_ = NULL;
		typeSin_ = NULL;
	}
}

/*
// Methods
*/

// Initialise all necessary arrays
bool Configuration::initialiseAtomArrays(int natoms)
{
	// Check that we haven't already initialised the arrays...
	if (r_ != NULL)
	{
		printf("Error: Atom-related arrays have already been initialised for this configuration.\n");
		return FALSE;
	}
	nAtoms_ = natoms;

	// Create arrays
	r_ = new Vec3<double>[nAtoms_];
	f_ = new Vec3<double>[nAtoms_];
	fCalc_ = new Vec3<double>[nAtoms_];
	fWeight_ = new double[nAtoms_];
	fContributions_ = new Vec3<double>*[Configuration::nContributionTypes];
	for (int i=0; i<Configuration::nContributionTypes; ++i) fContributions_[i] = new Vec3<double>[nAtoms_];
	for (int n=0; n<nAtoms_; ++n)
	{
		r_[n].zero();
		f_[n].zero();
		fWeight_[n] = 0.0;
		for (int i=0; i<Configuration::nContributionTypes; ++i) fContributions_[i][n].zero();
	}
	return TRUE;
}

// Delete existing atom arrays
void Configuration::deleteAtomArrays()
{
	if (r_ != NULL) delete[] r_;
	if (f_ != NULL) delete[] f_;
	if (fCalc_ != NULL) delete[] fCalc_;
	if (fWeight_ != NULL) delete[] fWeight_;
	if (fContributions_ != NULL)
	{
		for (int n=0; n<Configuration::nContributionTypes; ++n) delete[] fContributions_[n];
		delete[] fContributions_;
	}
	r_ = NULL;
	f_ = NULL;
	fCalc_ = NULL;
	fWeight_ = NULL;
	fContributions_ = NULL;
}

// Replicate atom/force data
bool Configuration::replicate(Vec3<int> rep, Species *firstsp)
{
	// Create temporary arrays
	int multiplier = rep.x*rep.y*rep.z;
	int newNAtoms = nAtoms_*multiplier;
	Vec3<double> *newR = new Vec3<double>[newNAtoms];
	Vec3<double> *newF = new Vec3<double>[newNAtoms];
	double *newFWeight = new double[newNAtoms];
	Vec3<double> transvec;

	// Duplicate atoms one species at a time
	int i = 0, j = 0, n, m, x, y, z, spstart = 0;
	Vec3<double> r;
	for (Species *sp = firstsp; sp != NULL; sp = sp->next)
	{
		if (sp->nAtoms() != 1)
		{
			printf("Error: Replicate cannot be used for non-monoatomic species, since a molecular fold has not yet been implemented.\n");
			return FALSE;
		}
		Vec3<int> disp;
		for (disp.x = 0; disp.x < rep.x; ++disp.x)
			for (disp.y = 0; disp.y < rep.y; ++disp.y)
				for (disp.z = 0; disp.z < rep.z; ++disp.z)
				{
					// Calculate displacement vector
					r.set(disp.x*1.0, disp.y*1.0, disp.z*1.0);
					r = cell_*r;
					j = spstart;
					for (m=0; m<sp->nMolecules(); ++m)
					{
						for (n=0; n<sp->nAtoms(); ++n)
						{
							newR[i] = r_[j+n]+r;
							newF[i] = f_[j+n];
							newFWeight[i] = fWeight_[j+n];
							++i;
						}
						j += sp->nAtoms();
					}
				}
		// Increase species starting atom counter
		spstart += sp->nAtoms()*sp->nMolecules();
	}
	
	// Remove old arrays, create new ones of the new size, and store new data
	deleteAtomArrays();
	initialiseAtomArrays(newNAtoms);
	for (n=0; n<newNAtoms; ++n)
	{
		r_[n] = newR[n];
		f_[n] = newF[n];
		fWeight_[n] = newFWeight[n];
	}
	
	// Adjust cell size
	cell_.columnMultiply(Vec3<double>(rep.x, rep.y, rep.z));
	
	// Delete temporary arrays
	delete[] newR;
	delete[] newF;
	delete[] newFWeight;
	
	return TRUE;
}

// Calculate Ewald sum cutoff
void Configuration::calculateEwaldCutoff(Vec3<int> kmax)
{
	// Reciprocal cutoff for Ewald sum is the shortest component of kVec * perpendicular reciprocal cell lengths
	Vec3<double> cross_ab = reciprocalCell_.columnAsVec3(0) * reciprocalCell_.columnAsVec3(1);
	Vec3<double> cross_bc = reciprocalCell_.columnAsVec3(1) * reciprocalCell_.columnAsVec3(2);
	Vec3<double> cross_ca = reciprocalCell_.columnAsVec3(2) * reciprocalCell_.columnAsVec3(0);
	Vec3<double> perpl(reciprocalVolume_ / cross_ab.magnitude(), reciprocalVolume_ / cross_bc.magnitude(), reciprocalVolume_ / cross_ca.magnitude());
	perpl.x *= kmax.x;
	perpl.y *= kmax.y;
	perpl.z *= kmax.z;
	reciprocalCutoff_ = perpl.min() * 1.05 * 2.0 * PI;
}

// Reset calculated forces (and energy)
void Configuration::reset(Configuration::ContributionType type)
{
	if (fContributions_ == NULL) printf("BUG: Force array does not exist to be zeroed.\n");
	else for (int n=0; n<nAtoms_; ++n) fContributions_[type][n].zero();
	energies_[type] = 0.0;
	logPoint_[type] = 0;
}

// Print configuration information, including calculated forces/deltas and energies
void Configuration::print(TargetSystem *targetSystem)
{
	Vec3<double> fdelta, fcalc;
	printf("Configuration '%s':\n", name_.get());
	printf("\t          vdW : %f\n", energies_[Configuration::VdwContribution]);
	printf("\tElectrostatic : %f\n", energies_[Configuration::ChargeContribution]);
	printf("\t         Bond : %f\n", energies_[Configuration::BondContribution]);
	printf("\t        Angle : %f\n", energies_[Configuration::AngleContribution]);
	printf("\t      Torsion : %f\n", energies_[Configuration::TorsionContribution]);
	printf("Sample of configuration coordinates and forces:\n");
	printf("Atom    RX       RY       RZ        FX           FY           FZ           CX           CY           CZ     \n");
	int skip = (nAtoms_/10 == 0 ? 1 : nAtoms_/10);
	skip = 1;
	for (int n=0; n<nAtoms_; n += skip)
	{
		fcalc = 0.0;
		for (int type = 0; type<Configuration::nContributionTypes; ++type) fcalc += fContributions_[type][n];
		fdelta = f_[n] - fcalc;
		printf("%-4i %-5s %8.3f %8.3f %8.3f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", n+1, targetSystem->atomTypeName(n), r_[n].x, r_[n].y, r_[n].z, f_[n].x, f_[n].y, f_[n].z, fcalc.x, fcalc.y, fcalc.z);
		printf("                                                                             %12.5e %12.5e %12.5e\n", fdelta.x, fdelta.y, fdelta.z);
	}
}

// Print energy information only
void Configuration::printEnergy()
{
	printf("\nConfiguration '%s':\n", name_.get());
	printf("\t          vdW : %f\n", energies_[Configuration::VdwContribution]);
	printf("\tElectrostatic : %f\n", energies_[Configuration::ChargeContribution]);
	printf("\t         Bond : %f\n", energies_[Configuration::BondContribution]);
	printf("\t        Angle : %f\n", energies_[Configuration::AngleContribution]);
	printf("\t      Torsion : %f\n", energies_[Configuration::TorsionContribution]);
}
