#ifndef MULTIAPPFISSIONPDFTRANSFER_H
#define MULTIAPPFISSIONPDFTRANSFER_H

#include "MultiAppTransfer.h"

class MultiAppFissionPDFTransfer;

template<>
InputParameters validParams<MultiAppFissionPDFTransfer>();

/**
 * Transfer a neutronics probability density to a neutronics based pka generator.
 */
class MultiAppFissionPDFTransfer : public MultiAppTransfer
{
public:
  MultiAppFissionPDFTransfer(const InputParameters & parameters);

  virtual void execute();

protected:
  UserObjectName _pka_generator_name;
  UserObjectName _neutronics_pdf_name;
};

#endif // MULTIAPPFISSIONPDFTRANSFER_H
