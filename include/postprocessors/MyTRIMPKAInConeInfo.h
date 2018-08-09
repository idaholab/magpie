/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMPKAINCONEINFO_H
#define MYTRIMPKAINCONEINFO_H

#include "MyTRIMPKAInfo.h"

// forward declarations
class MyTRIMPKAInConeInfo;

template <>
InputParameters validParams<MyTRIMPKAInConeInfo>();

class MyTRIMPKAInConeInfo : public MyTRIMPKAInfo
{
public:
  MyTRIMPKAInConeInfo(const InputParameters & parameters);

protected:
  virtual bool skipPKA(const MyTRIM_NS::IonBase & ion) const override;

  const RealVectorValue _direction;
  const Real _min_cosine;
};

#endif // MyTRIMPKAInConeInfo_H
