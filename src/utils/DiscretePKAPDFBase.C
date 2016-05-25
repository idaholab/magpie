/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "DiscretePKAPDFBase.h"
#include "MooseError.h"

DiscretePKAPDFBase::DiscretePKAPDFBase(Real magnitude) :
    _magnitude(magnitude)
{
}
