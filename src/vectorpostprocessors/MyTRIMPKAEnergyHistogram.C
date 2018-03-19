/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MyTRIMPKAEnergyHistogram.h"
#include "MyTRIMRasterizer.h"
#include "mytrim/ion.h"

template<>
InputParameters validParams<MyTRIMPKAEnergyHistogram>()
{
  InputParameters params = validParams<GeneralVectorPostprocessor>();
  params.addClassDescription("Generate an energy histogram for the primary knock-on atom (PKA) list");
  params.addRequiredParam<UserObjectName>("rasterizer", "Name of the MyTRIMRasterizer userobject to pull data from");
  params.addParam<unsigned int>("channel_number", 50, "Number of energy channels");
  params.addRangeCheckedParam<Real>("channel_width", 5.0e6, "channel_width > 0", "Energy channel width in eV");
  return params;
}

MyTRIMPKAEnergyHistogram::MyTRIMPKAEnergyHistogram(const InputParameters & params) :
    GeneralVectorPostprocessor(params),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _nchannels(getParam<unsigned int>("channel_number")),
    _deltaE(getParam<Real>("channel_width")),
    _channel_center(declareVector("E")),
    _count(declareVector("n"))
{
  // initialize the channel energy vector
  _channel_center.resize(_nchannels);
  for (unsigned i = 0; i < _nchannels; ++i)
    _channel_center[i] = (i + 0.5) * _deltaE;
}

void
MyTRIMPKAEnergyHistogram::initialize()
{
  // reset the statistics
  _count.assign(_nchannels, 0.0);
}

void
MyTRIMPKAEnergyHistogram::execute()
{
  const std::vector<MyTRIM_NS::IonBase> & pka_list = _rasterizer.getPKAList();
  for (auto & pka : pka_list)
  {
    int channel = pka._E / _deltaE;

    if (channel >= 0 && channel < _nchannels)
      _count[channel]++;
  }
}

void
MyTRIMPKAEnergyHistogram::finalize()
{
  gatherSum(_count);
}
