/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MYTRIMDIRACRUN_H
#define MYTRIMDIRACRUN_H

#include "MyTRIMRunBase.h"
#include "ThreadedRecoilDiracSourceLoop.h"

class MyTRIMDiracRun;

template<>
InputParameters validParams<MyTRIMDiracRun>();

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library and
 * stores the resulting defect distributions as
 */
class MyTRIMDiracRun : public MyTRIMRunBase
{
public:
  MyTRIMDiracRun(const InputParameters & parameters);

  virtual void initialize() {}
  virtual void execute();
  virtual void finalize();

  /// @{ shorthand typedefs
  typedef ThreadedRecoilDiracSourceLoop::MyTRIMResult MyTRIMResult;
  typedef ThreadedRecoilDiracSourceLoop::MyTRIMResultList MyTRIMResultList;
  /// @}

  /// get the TRIM result data
  const MyTRIMResultList & result() const;

protected:
  ///@{ pack/unpack the _result_map into a structure suitable for parallel communication
  void serialize(std::string & serialized_buffer);
  void deserialize(std::vector<std::string> & serialized_buffers);
  ///@}

  /// data such as interstitials and vacancies produced will be stored here
  MyTRIMResultList _result_list;
};

#endif //MYTRIMDIRACRUN_H
