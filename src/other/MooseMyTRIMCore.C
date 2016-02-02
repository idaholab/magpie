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
  _vac.push_back(std::make_pair(Point(recoil->pos(0), recoil->pos(1), _dim == 2 ? 0.0 :  recoil->pos(2)), recoil->tag));
}

void
MooseMyTRIMCore::checkPKAState()
{
}

void
MooseMyTRIMCore::dissipateRecoilEnergy()
{
}
