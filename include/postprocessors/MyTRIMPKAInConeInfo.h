/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMPKAInfo.h"

class MyTRIMPKAInConeInfo : public MyTRIMPKAInfo
{
public:
  static InputParameters validParams();

  MyTRIMPKAInConeInfo(const InputParameters & parameters);

protected:
  virtual bool skipPKA(const MyTRIM_NS::IonBase & ion) const override;

  const RealVectorValue _direction;
  const Real _min_cosine;
};
