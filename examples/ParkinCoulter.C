//
// clang++ -g --std=c++11 -o ParkinCoulter ParkinCoulter.C
//

#include <cmath>
#include <iostream>
#include <vector>

#define _console std:: cout

typedef double Real;
Real
sqr(Real a)
{
  return a * a;
}

const unsigned int nelem = 2;
std::vector<std::vector<std::vector<Real>>> g(nelem + 1);

// Lattice A (final entry is PKA)
std::vector<Real> A = {30.0, 70.0, 20.0};
std::vector<Real> Z = {15.0, 30.0, 10.0};

// atomic fraction
std::vector<Real> f = {0.33333, 0.66666};

// binding energy
std::vector<Real> Eb = {0.0, 0.0};

// capture energy (diagonals are displacement energies)
std::vector<std::vector<Real>> Ecap = {{30.0, 0.0, 0.0},
                                       {0.0, 30.0, 0.0},
                                       {0.0, 0.0, 0.0}};

// Eq. (2a)
Real
Rho(unsigned int k, Real T)
{
  return T > Ecap[k][k];
}

// Eq. (2b)
Real
lambda(unsigned int i, unsigned int k, Real T)
{
  return T < Ecap[i][k];
}

// electronic stopping (MyTRIM provides this)
Real
s(Real E)
{
  return 1.0;
}

// differential crossection (Huang/Ghoniem)
Real
dRhoikdT(unsigned int i, unsigned int k, Real E, Real T)
{
  const Real ep1 = 15.0, ep2 = 0.369, ep3 = 0.0234;

  // charge in the right units!
  const Real e = 1.0;

  // tabulate Z^2/3!
  const Real a = 0.4683 / std::sqrt(std::pow(Z[i], 2.0 / 3.0) + std::pow(Z[k], 2.0 / 3.0));
  const Real e0 = Z[i] * Z[k] * e * e / a;
  const Real ep = E / e0;

  Real m, lm;
  if (ep < ep3)
  {
    m = 0;
    lm = 24.0;
  }
  else if (ep <= ep2)
  {
    m = 1.0 / 3.0;
    lm = 1.309;
  }
  else if (ep <= ep1)
  {
    m = 0.5;
    lm = 0.327;
  }
  else
  {
    m = 1.0;
    lm = 0.5;
  }

  const Real Cm = 3.14159265359 / 2.0 * lm * a * a * A[i] / A[k] * std::pow(2.0 * Z[i] * Z[k] * e * e / a, 2.0 * m);

  return Cm * std::pow(E, -m) * (-1.0 - m) * std::pow(T, -2.0 - m); // d/dT of (5) in Huang/Ghoniem 1992
};

// maximum energy and binning
const Real Emax = 1000.0;
const unsigned int nbin = 1000;
// bin size (function of Emax an nbin)
const Real dE = Emax / Real(nbin);

Real
gij(unsigned int i, unsigned int j, Real E)
{
  return E > g[i][j][E / dE] ?: 0.0;
}

int
main()
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
            integral += dE * dRhoikdT(i, k, E, T) * (Rho(k, T) * gij(k, j, T - Eb[k]) + (1.0 - Rho(k, T) * lambda(i, k, E - T)) * gij(i, j, E - T) - gij(i, j, E));

          sum += f[k] * integral;
        }

        // integrate gij
        g[i][j][bin + 1] = g[i][j][bin] + dE * sum / s(E);
      }
    }
  }

  // output PKA data
  for (unsigned int j = 0; j < nelem; ++j)
  {
    if (j > 0)
      _console << '\n';

    for (unsigned int bin = 0; bin < nbin; ++bin)
    {
      // energy
      const Real E = (bin * Emax) / nbin;
      _console << E << ' ' << g[nelem][j][bin] << '\n';
    }
  }
}
