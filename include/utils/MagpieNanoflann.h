/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "libmesh/nanoflann.hpp"

// Make newer nanoflann API compatible with older nanoflann versions
#if NANOFLANN_VERSION < 0x150
namespace nanoflann
{
template <typename T, typename U>
using ResultItem = std::pair<T, U>;
}
namespace nanoflann
{
typedef SearchParams SearchParameters;
}
#endif
