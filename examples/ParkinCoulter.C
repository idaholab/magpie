//
// clang++ -g --std=c++11 -o ParkinCoulter ParkinCoulter.C
//

#include <iostream>
#include <vector>

typedef double Real;
Real sqr(Real a) { return a * a; }

const unsigned int nelem = 2;
std::vector<std::vector<std::vector<Real>>> g(nelem + 1);

// Lattice A (final entry is PKA)
std::vector<Real> A = { 30.0, 70.0, 20.0 };

// atomic fraction
std::vector<Real> f = { 0.33333, 0.66666 };

// binding energy
std::vector<Real> Eb = { 0.0, 0.0 };

// capture energy (diagonals are displacement energies)
std::vector<std::vector<Real>> Ecap = { { 30.0,  0.0, 0.0 },
                                        {  0.0, 30.0, 0.0 },
                                        {  0.0,  0.0, 0.0 } };

// Eq. (2a)
Real Rho(unsigned int k, Real T) { return T > Ecap[k][k]; }

// Eq. (2b)
Real lambda(unsigned int i, unsigned int k, Real T) { return T < Ecap[i][k]; }


// electronic stopping (MyTRIM provides this)
Real s(Real E) { return 1.0; }

// differential crossection (MyTRIM probably provides this, too)
Real dRhoikdT(Real E, Real T) { return 1.0; };


// maximum energy and binning
const Real Emax = 1000.0;
const unsigned int nbin = 1000;
// bin size (function of Emax an nbin)
const Real dE = Emax / Real(nbin);

Real gij(unsigned int i, unsigned int j, Real E) { return E > g[i][j][E/dE] ? : 0.0; }

int main()
{
  // initialize dg and g tables
  for (unsigned int i = 0; i <= nelem; ++i)
  {
    g[i].resize(nelem);
    for (unsigned int j = 0; j < nelem; ++j)
      g[i][j].assign(nbin + 1, 0.0);
  }

  // Tabulate energy
  for (unsigned int bin = 0; bin < nbin; ++bin)
  {
    // energy
    const Real E = (bin * Emax) / nbin;

    // potential knock on atioms (0..nelem-1 are lattice atoms, nelem is the external PKA)
    for (unsigned int i = 0; i <= nelem; ++i)
    {
      for (unsigned int j = 0; j < nelem; ++j)
      {
        Real sum = 0.0;
        for (unsigned int k = 0; k < nelem; ++k)
        {
          Real integral = 0.0;
          const Real MikE = E * (2.0 * A[i] * A[j]) / sqr(A[i] + A[j]);
          for (Real T = 0; T < MikE; T += dE)
            integral += dE * dRhoikdT(E,T) * (Rho(k,T) * gij(k,j,T - Eb[k]) + (1.0 - Rho(k,T) * lambda(i,k,E-T)) * gij(i,j, E-T) - gij(i,j,E));

          sum += f[k] * integral;
        }

        // integrate gij
        g[i][j][bin+1] = g[i][j][bin] + dE * sum / s(E);
      }
    }
  }

  // output PKA data
  for (unsigned int j = 0; j < nelem; ++j)
  {
    if (j > 0)
      std::cout << '\n';

    for (unsigned int bin = 0; bin < nbin; ++bin)
    {
      // energy
      const Real E = (bin * Emax) / nbin;
      std::cout << E << ' ' << g[nelem][j][bin] << '\n';
    }
  }
}
