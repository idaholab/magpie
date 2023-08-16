/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralPostprocessor.h"
#include "mytrim/ion.h"

class MyTRIMRasterizer;

class MyTRIMPKAInfo : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  MyTRIMPKAInfo(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

  using Postprocessor::getValue;
  virtual Real getValue() const override;

protected:
  /// override this function to add conditions for the considered PKAs
  virtual bool skipPKA(const MyTRIM_NS::IonBase & /*ion*/) const { return false; }

  /// the rasterizer object to pull data from
  const MyTRIMRasterizer & _rasterizer;

  /// which quantity is to be aggregated;
  enum ValueType
  {
    TOTAL_MASS = 0,
    TOTAL_ENERGY,
    TOTAL_CHARGE,
    TOTAL_NUMBER
  } _value_type;

  /// the value this PP computes (determined by the value_type parameter)
  Real _value;
};
