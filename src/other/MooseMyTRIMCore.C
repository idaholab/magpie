#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMSample.h"

MooseMyTRIMCore::MooseMyTRIMCore(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample,
                                 std::list<std::pair<Point, unsigned int> > & vac_list) :
    MyTRIM_NS::TrimBase(simconf, sample),
    _vac_list(vac_list),
    _dim(sample->getDim())
{
}

void
MooseMyTRIMCore::vacancyCreation()
{
  // called if both atoms in the recoil have sufficient energy to leave the site
  _vac_list.push_back(std::make_pair(Point(_recoil->_pos(0), _recoil->_pos(1), _dim == 2 ? 0.0 :  _recoil->_pos(2)), _recoil->_tag));
}
