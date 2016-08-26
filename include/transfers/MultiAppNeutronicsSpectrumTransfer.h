#ifndef MULTIAPPNEUTRONICSSPECTRUMTRANSFER_H
#define MULTIAPPNEUTRONICSSPECTRUMTRANSFER_H

#include "MultiAppTransfer.h"

class MultiAppNeutronicsSpectrumTransfer;

template<>
InputParameters validParams<MultiAppNeutronicsSpectrumTransfer>();

/**
 * Transfer a neutronics probability density to a neutronics based pka generator.
 */
class MultiAppNeutronicsSpectrumTransfer : public MultiAppTransfer
{
public:
  MultiAppNeutronicsSpectrumTransfer(const InputParameters & parameters);

  virtual void execute();

protected:
  UserObjectName _pka_generator_name;
  UserObjectName _neutronics_pdf_name;
};

#endif //MULTIAPPNEUTRONICSSPECTRUMTRANSFER_H
