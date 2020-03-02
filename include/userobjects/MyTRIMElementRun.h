/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MyTRIMRunBase.h"
#include "ThreadedRecoilElementAveragedLoop.h"

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library
 */
class MyTRIMElementRun : public MyTRIMRunBase
{
public:
  static InputParameters validParams();

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

  ///@{ timers
  PerfID _perf_trim;
  PerfID _perf_finalize;
  ///@}

private:
  /// zero result to return for elements that have not been touched by the cascades
  MyTRIMResult _zero;
};
