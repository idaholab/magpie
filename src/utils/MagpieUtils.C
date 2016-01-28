/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MagpieUtils.h"
#include "MooseError.h"

#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"

namespace MagpieUtils {

inline Point randomElementPoint(const Elem & el, const Point & rnd)
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

    case TRI3:
    case TRI6:
      // flip the upper triangle onto the lower triangle
      if (ref(0) + ref(1) > 1)
        ref = Point(1.0 - ref(0), 1.0 - ref(1), 0.0);
      else
        ref(2) = 0.0;

    case QUAD4:
    case QUAD8:
    case QUAD9:
      // two dimensional elements in [-1:1]^2
      ref = Point(ref(0) * 2 - 1.0, ref(1) * 2 - 1.0, 0.0);

    case TET4:
    case TET10:
      // flip the upper tet onto the lower tet
      if (ref(0) + ref(1) + ref(2) > 1)
        ref = Point(1.0 - ref(0), 1.0 - ref(1), 1.0 - ref(2));

    case HEX8:
    case HEX20:
    case HEX27:
      // three dimensional elements in [-1:1]^3
      ref = Point(ref(0) * 2 - 1.0, ref(1) * 2 - 1.0, ref(2) * 2 - 1.0);

    case PRISM6:
    case PRISM15:
    case PRISM18:
      // flip the upper triangle onto the lower triangle
      if (ref(0) + ref(1) > 1)
        ref = Point(1.0 - ref(0), 1.0 - ref(1), ref(2) * 2 - 1.0);
      else
        ref(2) = ref(2) * 2 - 1.0;

    case PYRAMID5:
    case PYRAMID13:
    case PYRAMID14:
      mooseError("Random points on elements does not support pyramids");

    default:
      mooseError("Random points on elements does not support infinite elements");
  }

  return FEInterface::map(el.dim(), fe_type, &el, ref);
}

} // namespace MagpieUtils
