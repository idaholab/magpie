#include "MooseMyTRIMCore.h"
#include "MooseMyTRIMSample.h"

MooseMyTRIMCore::MooseMyTRIMCore(MyTRIM_NS::SimconfType * simconf, MooseMyTRIMSample * sample, std::vector<std::pair<Point, unsigned int> > & vac) :
    MyTRIM_NS::TrimBase(simconf, sample),
    _vac(vac),
    _dim(sample->getDim())
{
}

void
MooseMyTRIMCore::vacancyCreation()
{
  // called if both atoms in the recoil have sufficient energy to leave the site
  _vac.push_back(std::make_pair(Point(_recoil->_pos(0), _recoil->_pos(1), _dim == 2 ? 0.0 :  _recoil->_pos(2)), _recoil->_tag));
}

void
MooseMyTRIMCore::checkPKAState()
{
}

void
MooseMyTRIMCore::dissipateRecoilEnergy()
{
}
