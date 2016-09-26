#ifndef MOOSEMYTRIMCORE_H
#define MOOSEMYTRIMCORE_H

#include "mytrim/trim.h"
#include "MooseTypes.h"

#include <list>

class MooseMyTRIMSample;

/**
 * MyTRIM simulation core class
 */
class MooseMyTRIMCore : public MyTRIM_NS::TrimBase
{
public:
  MooseMyTRIMCore(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample, std::list<std::pair<Point, unsigned int> > & vac_list);

  virtual void vacancyCreation();

protected:
  /// list of vacancies generated during the recoil
  std::list<std::pair<Point, unsigned int> > & _vac_list;

  // dimension of the mesh
  const unsigned int _dim;
};

#endif //MOOSEMYTRIMCORE_H
