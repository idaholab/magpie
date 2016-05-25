/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "DiscretePDFTest.h"

//Magpie includes
#include "DiscretePKAPDF.h"
#include "MultiIndex.h"

#include <cmath>

CPPUNIT_TEST_SUITE_REGISTRATION( DiscretePDFTest );

void
DiscretePDFTest::samplePKA()
{
  // Construct a very simple pdf: 3 ZAIDs, 2 energies, 4 azimuthal regions, 5 polar regions
  // For simplicity it'll be a product of 1-D pdfs P{ZAID} * P{E} * P{phi} * P{mu}
  std::vector<unsigned int> zaid(2);
  std::vector<Real> pzaid(2);
  zaid[0] = 10010;
  zaid[1] = 160320;
  pzaid[0] = 0.3;
  pzaid[1] = 0.7;

  std::vector<Real> energies(3);
  std::vector<Real> penergies(3);
  energies[0] = 0.0;
  energies[1] = 1.0;
  energies[2] = 9.0;
  penergies[0] = 0.4;
  penergies[1] = 0.6 / 8.0;

  Real pphi = 1.0 / (2.0 * libMesh::pi);
  Real dphi = 2.0 * libMesh::pi / 4.0;
  std::vector<Real> pmu(3);
  Real dmu = 2.0 / 3.0;
  pmu[0] = 0.2;
  pmu[1] = 0.5;
  pmu[2] = 0.8;

  std::vector<unsigned int> shape(4);
  shape[0] = 2;
  shape[1] = 2;
  shape[2] = 4;
  shape[3] = 3;
  MultiIndex<Real> mindex = MultiIndex<Real>(shape);

  std::vector<unsigned int> index(4);
  for (unsigned int j0 = 0; j0 < shape[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape[1]; ++j1)
      for (unsigned int j2 = 0; j2 < shape[2]; ++j2)
        for (unsigned int j3 = 0; j3 < shape[3]; ++j3)
        {
          index[0] = j0;
          index[1] = j1;
          index[2] = j2;
          index[3] = j3;
          mindex(index) = pzaid[j0] * penergies[j1] * pphi * pmu[j3];
        }

  DiscretePKAPDF pdf = DiscretePKAPDF(1.0, zaid, energies, 4, 3, mindex);

  // set up a frequency counter
  unsigned int j1, j2, j3, j4;
  MultiIndex<Real> frequency = MultiIndex<Real>(shape);

  // test it by sampling
  Real counter = 0.0, max = 2.0e5;
  while (counter < max)
  {
    DiscretePKAPDFBase::initialPKAState i_state;
    pdf.drawSample(i_state);
    counter += 1.0;

    // find the correct "bin"
    if (i_state._Z == 1)
      j1 = 0;
    else if (i_state._Z == 16)
      j1 = 1;
    else
      mooseError("Incorrect Z in initialPKAState.");

    if (i_state._energy <= 1.0)
      j2 = 0;
    else if (i_state._energy > 1.0 && i_state._energy <= 9.0)
      j2 = 1;
    else
      mooseError("Incorrect energy in initialPKAState.");

    Real mu = i_state._direction(0);
    j4 = std::floor((mu + 1.0) / dmu);

    Real phi = std::atan2(i_state._direction(2), i_state._direction(1));
    j3 = std::floor((phi + libMesh::pi) / dphi);

    index[0] = j1;
    index[1] = j2;
    index[2] = j3;
    index[3] = j4;
    frequency(index) += 1.0 / max;
  }

  for (MultiIndex<Real>::iterator it = frequency.begin(); it !=  frequency.end() ; ++it)
  {
    index = it.indices();
    Real width = (energies[index[1] + 1] - energies[index[1]]) * dmu * dphi;
    CPPUNIT_ASSERT( std::abs(1.0 - mindex(index) * width / frequency(index)) < 0.1);
  }
}
