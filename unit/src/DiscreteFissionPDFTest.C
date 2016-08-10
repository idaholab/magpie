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

#include "DiscreteFissionPDFTest.h"

//Magpie includes
#include "DiscreteFissionPKAPDF.h"
#include "PKAGeneratorBase.h"
#include "MultiIndex.h"
#include "mytrim/ion.h"

#include <cmath>
#include <iostream>

CPPUNIT_TEST_SUITE_REGISTRATION( DiscreteFissionPDFTest );

void
DiscreteFissionPDFTest::sampleFissionPKA()
{
  //add environment variable for path of sum yield data files
  setenv("ENDF_FP_DIR","../data/fission_yield/",1);

  // Define vector of zaids and probabilites
  std::vector<unsigned int> zaid(2);
  std::vector<Real> pzaid(2);
  zaid[0] = 922350;
  zaid[1] = 942390;
  pzaid[0] = 0.0;
  pzaid[1] = 1.0;

  //Define energies and probabilites
  std::vector<Real> energies(3);
  std::vector<Real> penergies(2);
  energies[0] = 1.0e-6;
  energies[1] = 0.5;
  energies[2] = 9.0;
  penergies[0] = 1.0 / (0.5-1.0e-6);
  penergies[1] = 0.0 / 8.5;

  //Define MultiIndex probabilites
  std::vector<unsigned int> shape(2);
  shape[0] = 2;
  shape[1] = 2;
  MultiIndex<Real> mindex = MultiIndex<Real>(shape);

  //"Build" MultiIndex
  std::vector<unsigned int> index(2);
  for (unsigned int j0 = 0; j0 < shape[0]; ++j0)
    for (unsigned int j1 = 0; j1 < shape[1]; ++j1)
    {
      index[0] = j0;
      index[1] = j1;
      mindex(index) = pzaid[j0] * penergies[j1];
    }

  DiscreteFissionPKAPDF pdf = DiscreteFissionPKAPDF(zaid, energies, mindex);

  // set up a frequency counter
  unsigned int j1, j2;
  MultiIndex<Real> frequency = MultiIndex<Real>(shape);

  // test it by sampling
  Real counter = 0.0, max = 5.0e5;
  Real total_energy = 0.0;
  Real total_mass = 0.0;
  Real total_Z = 0.0;
  Real mubin1 = 0.0, mubin2 = 0.0, mubin3 = 0.0, mubin4 = 0.0;
  Real phibin1 = 0.0, phibin2 = 0.0, phibin3 = 0.0, phibin4 = 0.0;
  Real true_energy = 173.5265049499354998486; //for Pu-239

  //create std::vector for A values of fission products, corresponding tallies, and true pdf
  std::vector<unsigned int> fproduct_A;
  std::vector<Real> fproduct_pdf;

  //read in pdf for Pu-239 Epithermal
  std::string filename = "./data/Pu239_A.txt";
  std::ifstream infile(filename);
  int aa;
  Real bb;
  while (infile >> aa >> bb)
  {
    fproduct_A.push_back(aa);
    fproduct_pdf.push_back(bb);
  }

  std::vector<Real> fproduct_tally(fproduct_pdf.size());
  unsigned int tally_index;

  while (counter < max)
  {
    std::vector<MyTRIM_NS::IonBase> i_state;
    pdf.drawSample(i_state);
    counter += 1.0;

    //find matching bin to the sample A
    for (unsigned int i = 0; i < fproduct_A.size(); ++i)
    {
      if ( int(i_state[0]._m) == fproduct_A[i])
        {
          tally_index = i;
          break;
        }
    }
    //tally mass number in corresponding bin
    fproduct_tally[tally_index] += 1.0;

    total_energy += i_state[0]._E + i_state[1]._E;
    total_mass += i_state[0]._m + i_state[1]._m + 2.8836;
    total_Z += i_state[0]._Z + i_state[1]._Z;

    DiscreteFissionPDFTest::setRandomDirection(i_state[0]);
    i_state[1]._dir = -i_state[0]._dir;

    //sampled_mu (int)
    int sampled_mu = int((i_state[0]._dir(2) + 1.0) * 2.0);

    //mu bins
    if ( sampled_mu >= 0 && sampled_mu < 1)
      mubin1 += 1.0;
    else if (sampled_mu >= 1 && sampled_mu < 2)
      mubin2 += 1.0;
    else if (sampled_mu >= 2 && sampled_mu < 3)
      mubin3 += 1.0;
    else if (sampled_mu >= 3 && sampled_mu <= 4)
      mubin4 += 1.0;
    else
      mooseError("mu is not between -1 and 1");

    //retrieve phi
    Real sampled_phi = std::atan2(i_state[0]._dir(1),i_state[0]._dir(0));

    //phi bins
    if (sampled_phi >= -libMesh::pi && sampled_phi < -libMesh::pi / 2.0)
      phibin1 += 1.0;
    else if (sampled_phi >= -libMesh::pi / 2.0 && sampled_phi < 0.0)
      phibin2 += 1.0;
    else if (sampled_phi >= 0.0 && sampled_phi < libMesh::pi / 2.0)
      phibin3 += 1.0;
    else if (sampled_phi >= libMesh::pi / 2.0 && sampled_phi <= libMesh::pi)
      phibin4 += 1.0;
    else
      mooseError("phi is not between -pi and pi");

    index[0] = j1;
    index[1] = j2;
    frequency(index) += 1.0 / max;
  }

  //normalize the talley to pdf
  Real sum_tally = 0.0;
  for (auto n : fproduct_tally)
    sum_tally += n;

  sum_tally *= 0.5;

  for (unsigned int j = 0; j < fproduct_tally.size(); ++j)
    fproduct_tally[j] = fproduct_tally[j]/sum_tally;

  sum_tally = 0.0;
  for (auto n : fproduct_tally)
    sum_tally += n;

  //calculate the error for each mass number bin
  std::vector<Real> error(fproduct_tally.size());
  for (unsigned int j = 0; j < error.size(); ++j)
    error[j] = std::abs(fproduct_tally[j] - fproduct_pdf[j]);

  //find the maximum error
  auto maxerror = std::max_element(error.begin(), error.end());

  //check if total energy is conserved
  CPPUNIT_ASSERT( std::abs(total_energy / max - true_energy) < 1.0e-6 );

  //check if total mass is conserved
  CPPUNIT_ASSERT( std::abs(total_mass / max - 239.0) < 1.0e-3);

  //check if the z number is conserved
  CPPUNIT_ASSERT( std::abs(total_Z / max - 94.0) < 1.0e-6);

  //check if mu bins are uniformly distributed
  CPPUNIT_ASSERT( std::abs(mubin1 / max - 0.25) < 1.0e-2);
  CPPUNIT_ASSERT( std::abs(mubin2 / max - 0.25) < 1.0e-2);
  CPPUNIT_ASSERT( std::abs(mubin3 / max - 0.25) < 1.0e-2);
  CPPUNIT_ASSERT( std::abs(mubin4 / max - 0.25) < 1.0e-2);

  //check if phi bins are uniformly distributed
  CPPUNIT_ASSERT( std::abs(phibin1 / max - 0.25) < 1.0e-2);
  CPPUNIT_ASSERT( std::abs(phibin2 / max - 0.25) < 1.0e-2);
  CPPUNIT_ASSERT( std::abs(phibin3 / max - 0.25) < 1.0e-2);
  CPPUNIT_ASSERT( std::abs(phibin4 / max - 0.25) < 1.0e-2);

  //check max error of sampled distribution vs true distribution
  CPPUNIT_ASSERT( *maxerror < 1.0e-2);
}

void
DiscreteFissionPDFTest::setRandomDirection(MyTRIM_NS::IonBase & ion)
{
  Real nsq, x1, x2;

  // Marsaglia's method for uniformly sampling the surface of the sphere
  do
  {
    x1 = 2 * MooseRandom::rand() - 1.0;
    x2 = 2 * MooseRandom::rand() - 1.0;
    nsq = x1 * x1 + x2 * x2;
  } while (nsq >= 1);

  // construct normalized direction vector
  ion._dir(0) = 2.0 * x1 * std::sqrt(1.0 - nsq);
  ion._dir(1) = 2.0 * x2 * std::sqrt(1.0 - nsq);
  ion._dir(2) = 1.0 - 2.0 * nsq;
}
