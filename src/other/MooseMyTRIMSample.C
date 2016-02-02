#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

MooseMyTRIMSample::MooseMyTRIMSample(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh, MyTRIM_NS::SimconfType * simconf) :
    MyTRIM_NS::SampleBase(),
    _rasterizer(rasterizer),
    _nvars(_rasterizer.nVars()),
    _trim_mass(_rasterizer.mass()),
    _trim_charge(_rasterizer.charge()),
    _mesh(mesh),
    _dim(_mesh.dimension()),
    _pl(_mesh.getMesh().sub_point_locator()),
    _simconf(simconf)
{
  if (_dim < 2 || _dim > 3)
    mooseError("TRIM simulation works in 2D or 3D only.");
}

MooseMyTRIMSample::~MooseMyTRIMSample()
{
  // delete all elements of the materials in the cache
  for (MaterialsCache::iterator i = _materials_cache.begin(); i != _materials_cache.end(); ++i)
    for (unsigned int j = 0; j <  i->second.element.size(); ++j)
      delete i->second.element[j];
}

void
MooseMyTRIMSample::averages(const MyTRIM_NS::IonBase  * pka)
{
  _current_ion = pka;

  // apply averages for all cached materials
  for (MaterialsCache::iterator i = _materials_cache.begin(); i != _materials_cache.end(); ++i)
    i->second.average(pka);
}

MyTRIM_NS::MaterialBase *
MooseMyTRIMSample::lookupMaterial(Point & pos)
{
  // point to sample the material at
  Point p(pos(0), pos(1), _dim == 2 ? 0.0 : pos(2));

  // get element containing the point
  mooseAssert(_pl != NULL, "initialize() must be called on the MooseMyTRIMSample object.");
  const Elem * elem = (*_pl)(p);

  // no element found means we have left the mesh
  if (elem == NULL)
    return NULL;

  // try to find the element in the cache
  MaterialsCache::iterator i = _materials_cache.find(elem);
  if (i != _materials_cache.end())
    return &i->second;

  // otherwise prepare the material using data from the rasterizer
  const std::vector<Real> & material_data = _rasterizer.material(elem);
  MooseMyTRIMMaterial material(_simconf, 1.0);

  // set elements
  for (unsigned int i = 0; i < _nvars; ++i)
  {
    MyTRIM_NS::ElementBase * element = new MyTRIM_NS::ElementBase;
    element->_Z = _trim_mass[i];
    element->_m = _trim_charge[i];
    element->_t = material_data[i];
    material.element.push_back(element);
  }

  // prepare material
  material.prepare();
  material.average(_current_ion);

  // add the material to the cache and return pointer
  MaterialsCache::iterator mi = _materials_cache.insert(_materials_cache.begin(), std::make_pair(elem, material));
  return &mi->second;
}
