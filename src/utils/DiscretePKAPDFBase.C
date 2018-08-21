/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "DiscretePKAPDFBase.h"
#include "MooseError.h"
#include "MooseRandom.h"

DiscretePKAPDFBase::DiscretePKAPDFBase(const std::vector<unsigned int> & ZAID,
                                       const std::vector<Real> & energies)
  : _zaids(ZAID), _nZA(_zaids.size()), _energies(energies), _ng(_energies.size() - 1)
{
}

unsigned int
DiscretePKAPDFBase::sampleHelper(const MultiIndex<Real> & marginal_pdf,
                                 const MultiIndex<Real>::size_type indices) const
{
  if (marginal_pdf.dim() - indices.size() != 1)
    mooseError("For sampling the indices vector must reduce the marginal_pdf to a one-dimensional "
               "ladder function.");
  MultiIndex<Real>::size_type dimension(indices.size());
  for (unsigned int j = 0; j < indices.size(); ++j)
    dimension[j] = j;
  MultiIndex<Real> new_marginal_pdf = marginal_pdf.slice(dimension, indices);
  return sampleHelper(new_marginal_pdf);
}

unsigned int
DiscretePKAPDFBase::sampleHelper(const MultiIndex<Real> & marginal_pdf) const
{
  Real r = MooseRandom::rand();
  std::vector<unsigned int> index(1);
  unsigned int j = 0;
  for (; j < marginal_pdf.size()[0]; ++j)
  {
    index[0] = j;
    if (r <= marginal_pdf(index))
      break;
  }
  return j;
}

unsigned int
DiscretePKAPDFBase::sampleHelper(const std::vector<Real> & marginal_pdf) const
{
  Real r = MooseRandom::rand();
  std::vector<unsigned int> index(1);
  unsigned int j = 0;
  for (; j < marginal_pdf.size(); ++j)
  {
    if (r <= marginal_pdf[j])
      break;
  }
  return j;
}
