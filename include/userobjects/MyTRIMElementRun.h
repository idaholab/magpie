/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MYTRIMRUN_H
#define MYTRIMRUN_H

#include "MyTRIMRunBase.h"

class MyTRIMElementRun;

template<>
InputParameters validParams<MyTRIMElementRun>();

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library
 */
class MyTRIMElementRun : public MyTRIMRunBase
{
public:
  MyTRIMElementRun(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute();
  virtual void finalize();

  /// @{ shorthand typedefs
  typedef ThreadedRecoilElementAveragedLoop::MyTRIMResult MyTRIMResult;
  typedef ThreadedRecoilElementAveragedLoop::MyTRIMResultMap MyTRIMResultMap;
  /// @}

  /// get the TRIM result data
  const MyTRIMResult & result(const Elem *) const;

protected:
  ///@{ pack/unpack the _result_map into a structure suitable for parallel communication
  void serialize(std::string & serialized_buffer);
  void deserialize(std::vector<std::string> & serialized_buffers);
  ///@}

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultMap _result_map;

private:
  /// zero result to return for elements that have not been touched by the cascades
  MyTRIMResult _zero;
};

#endif //MYTRIMRUN_H
