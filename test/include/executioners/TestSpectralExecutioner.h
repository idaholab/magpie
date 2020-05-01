/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef TESTSPECTRALEXECUTIONER_H
#define TESTSPECTRALEXECUTIONER_H

#include "SpectralExecutionerBase.h"
#include "FFTWBufferBase.h"
#include "FFTProblem.h"

// System includes
#include <string>

// Forward declarations
class InputParameters;

/**
 * FFT Executioner base class.
 */
class TestSpectralExecutioner : public SpectralExecutionerBase
{
public:
  static InputParameters validParams();

  TestSpectralExecutioner(const InputParameters & parameters);
  virtual void execute() override;

};


#endif
