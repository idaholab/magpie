/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "ComplexTypes.h"
#include "MooseError.h"
#include "MooseUtils.h"

namespace MooseUtils
{

template <>
bool
absoluteFuzzyEqual(const std::complex<Real> & var1,
                   const std::complex<Real> & var2,
                   const std::complex<Real> & tol)
{
  return (std::abs(var1 - var2) <= std::abs(tol));
}

template <>
bool
absoluteFuzzyEqual(const std::complex<Real> & var1,
                   const Real & var2,
                   const std::complex<Real> & tol)
{
  return (std::abs(var1 - var2) <= std::abs(tol));
}

} // namespace MooseUtils

#include "RankTwoTensorImplementation.h"
#include "RankThreeTensorImplementation.h"
#include "RankFourTensorImplementation.h"

namespace libMesh
{
template class VectorValue<Complex>;
}

template class RankTwoTensorTempl<Complex>;
template class RankThreeTensorTempl<Complex>;
template class RankFourTensorTempl<Complex>;
