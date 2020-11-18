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

// This instantiation for complex variable simplifies logic and skips
// checks on user's input validity.
template <>
void
RankFourTensorTempl<Complex>::fillGeneralOrthotropicFromInputVector(
    const std::vector<Complex> & input)
{
  if (input.size() != 12)
    mooseError("To use fillGeneralOrhotropicFromInputVector, your input must have size 12. Yours "
               "has size ",
               input.size());

  const Complex & Ea = input[0];
  const Complex & Eb = input[1];
  const Complex & Ec = input[2];
  const Complex & Gab = input[3];
  const Complex & Gbc = input[4];
  const Complex & Gca = input[5];
  const Complex & nuba = input[6];
  const Complex & nuca = input[7];
  const Complex & nucb = input[8];
  const Complex & nuab = input[9];
  const Complex & nuac = input[10];
  const Complex & nubc = input[11];

  unsigned int ntens = N * (N + 1) / 2;

  std::vector<Complex> mat;
  mat.assign(ntens * ntens, 0);

  Complex k = 1 - nuab * nuba - nubc * nucb - nuca * nuac - nuab * nubc * nuca - nuba * nucb * nuac;

  mat[0] = Ea * (1 - nubc * nucb) / k;
  mat[1] = Ea * (nubc * nuca + nuba) / k;
  mat[2] = Ea * (nuba * nucb + nuca) / k;

  mat[6] = Eb * (nuac * nucb + nuab) / k;
  mat[7] = Eb * (1 - nuac * nuca) / k;
  mat[8] = Eb * (nuab * nuca + nucb) / k;

  mat[12] = Ec * (nuab * nubc + nuac) / k;
  mat[13] = Ec * (nuac * nuba + nubc) / k;
  mat[14] = Ec * (1 - nuab * nuba) / k;

  mat[21] = 2 * Gab;
  mat[28] = 2 * Gca;
  mat[35] = 2 * Gbc;

  int nskip = N - 1;

  unsigned int index = 0;
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      for (unsigned int k = 0; k < N; ++k)
        for (unsigned int l = 0; l < N; ++l)
        {
          if (i == j)
            (*this)._vals[index] =
                k == l ? mat[i * ntens + k] : mat[i * ntens + k + nskip + l] / 2.0;
          else
            (*this)._vals[index] = k == l ? mat[(nskip + i + j) * ntens + k]
                                          : mat[(nskip + i + j) * ntens + k + nskip + l] / 2.0;
          index++;
        }
}

template class RankTwoTensorTempl<Complex>;
template class RankThreeTensorTempl<Complex>;
template class RankFourTensorTempl<Complex>;
template RankTwoTensorTempl<Complex>
RankFourTensorTempl<Complex>::operator*(const RankTwoTensorTempl<Complex> & a) const;
