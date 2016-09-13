#ifndef MYTRIMDIRACRESULT_H
#define MYTRIMDIRACRESULT_H

#include "GeneralVectorPostprocessor.h"

class MyTRIMDiracResult;
class MyTRIMDiracRun;

template<>
InputParameters validParams<MyTRIMDiracResult>();

/**
 * Outputs the list of MyTRIM defects comopiled with the the MyTRIMDiracRunner
 */
class MyTRIMDiracResult : public GeneralVectorPostprocessor
{
public:
  MyTRIMDiracResult(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void threadJoin(const UserObject &);

protected:
  const MyTRIMDiracRun & _mytrim;
  const unsigned int _ivar;
  const unsigned int _defect;

  ///@{ coordinates of the defects
  VectorPostprocessorValue & _x;
  VectorPostprocessorValue & _y;
  VectorPostprocessorValue & _z;
  VectorPostprocessorValue & _elem_id;
  ///@}
};

#endif //MYTRIMDIRACRESULT_H
