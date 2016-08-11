#ifndef MULTIAPPRADIATIONDAMAGETRANSFER_H
#define MULTIAPPRADIATIONDAMAGETRANSFER_H

#include "MultiAppTransfer.h"

class MultiAppRadiationDamageTransfer;

template<>
InputParameters validParams<MultiAppRadiationDamageTransfer>();

/**
 * Transfer a neutronics probability density to a neutronics based pka generator.
 */
class MultiAppRadiationDamageTransfer : public MultiAppTransfer
{
public:
  MultiAppRadiationDamageTransfer(const InputParameters & parameters);

  virtual void execute();

protected:
  UserObjectName _pka_generator_name;
  UserObjectName _neutronics_pdf_name;
};

#endif // MULTIAPPRADIATIONDAMAGETRANSFER_H
