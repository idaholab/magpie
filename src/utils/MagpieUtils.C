/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MagpieUtils.h"
#include "MooseError.h"

#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"

#include <sstream>

namespace MagpieUtils
{

Point
randomElementPoint(const Elem & el, const Point & rnd)
{
  FEType fe_type(el.default_order());

  // generate point on the reference element
  Point ref(rnd);

  // we need to remap the random numbers according to the reference domain
  // as described in FEAbstract::on_reference_element
  switch (el.type())
  {
    case NODEELEM:
      // those are just single points
      ref = 0.0;
      break;

    case EDGE2:
    case EDGE3:
    case EDGE4:
      // one dimensional elements in [-1:1]
      ref = Point(ref(0) * 2 - 1.0, 0.0, 0.0);
      break;

    case TRI3:
    case TRI6:
      // flip the upper triangle onto the lower triangle
      if (ref(0) + ref(1) > 1)
        ref = Point(1.0 - ref(0), 1.0 - ref(1), 0.0);
      else
        ref(2) = 0.0;
      break;

    case QUAD4:
    case QUAD8:
    case QUAD9:
      // two dimensional elements in [-1:1]^2
      ref = Point(ref(0) * 2 - 1.0, ref(1) * 2 - 1.0, 0.0);
      break;

    case TET4:
    case TET10:
      // already in the reference volume
      if (ref(0) + ref(1) + ref(2) < 1)
        break;

      // three more tets in the corners of the cube
      else if (ref(0) + 1 - ref(1) + 1 - ref(2) < 1)
        ref = Point(ref(0), 1.0 - ref(1), 1.0 - ref(2));
      else if (1 - ref(0) + ref(1) + 1 - ref(2) < 1)
        ref = Point(1.0 - ref(0), ref(1), 1.0 - ref(2));
      else if (1 - ref(0) + 1 - ref(1) + ref(2) < 1)
        ref = Point(1.0 - ref(0), 1.0 - ref(1), ref(2));

      // internal tet pair
      {
        Real s = (ref.norm_sq() - 1) / 2.0;
        ref -= Point(s, s, s);
      }
      break;

    case HEX8:
    case HEX20:
    case HEX27:
      // three dimensional elements in [-1:1]^3
      ref = Point(ref(0) * 2 - 1.0, ref(1) * 2 - 1.0, ref(2) * 2 - 1.0);
      break;

    case PRISM6:
    case PRISM15:
    case PRISM18:
      // flip the upper triangle onto the lower triangle
      if (ref(0) + ref(1) > 1)
        ref = Point(1.0 - ref(0), 1.0 - ref(1), ref(2) * 2 - 1.0);
      else
        ref(2) = ref(2) * 2 - 1.0;
      break;

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
      mooseError("Random points on elements does not support pyramids");

    default:
      mooseError("Random points on elements does not support infinite elements");
  }

  return FEInterface::map(el.dim(), fe_type, &el, ref);
}

/**
 * ZAID is an identifier of an isotope of the form
 * ZZAAAm, where m is the state.
 */
unsigned int
getZFromZAID(unsigned int zaid)
{
  mooseAssert(zaid <= 999999 && zaid > 9999, "ZAID " << zaid << " is invalid.");
  return zaid / 10000;
}

unsigned int
getAFromZAID(unsigned int zaid)
{
  mooseAssert(zaid <= 999999 && zaid > 9999, "ZAID " << zaid << " is invalid.");
  return (zaid % 10000) / 10;
}

const std::string &
neutronEnergyName(unsigned int i)
{
  const static std::vector<std::string> names({"Thermal", "Epithermal", "Fast", "High"});
  return names[i];
}

NeutronEnergyType
determineNeutronType(Real energy)
{
  if (energy < 0.5)
    return Thermal;
  else if (energy >= 0.5 && energy < 0.75e6)
    return Epithermal;
  else if (energy >= 0.75e6 && energy < 7.0e6)
    return Fast;
  return High;
}

} // namespace MagpieUtils
