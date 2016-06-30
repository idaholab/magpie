#ifndef MAGPIEUTILS_H
#define MAGPIEUTILS_H

#include "Moose.h"
#include "MooseEnum.h"
#include "MultiMooseEnum.h"
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"

namespace MagpieUtils
{

/**
 * Returns a point inside element el distributed with uniform
 * probability. For higher order elements the non-linearity of the
 * mapping from the reference element coordinates to the physical element
 * could introduce a non-uniform distribution.
 * @param el The target element to locate the random point in
 * @param rnd A set of random numbers in [0:1) to use for generating the projected point
 */
Point randomElementPoint(const Elem & el, const Point & rnd);

///@{ ZAID (ZZAAA) conversion helper
unsigned int getZFromZAID(unsigned int zaid);
unsigned int getAFromZAID(unsigned int  zaid);
///@}

///enum of different neutron energy groups (thermal, epithermal, fast, high)
enum neutronEnergyTypes{Thermal, Epithermal, Fast, High, enum_max};

/// vector of strings for each neutron energy type
static const std::vector<std::string> neutronEnergyString = {"Thermal", "Epithermal", "Fast", "High"};

///Determines neutron energy type given it's energy
neutronEnergyTypes determineNeutronType(Real energy);

} // namespace MagpieUtils

#endif //MAGPIEUTILS_H
