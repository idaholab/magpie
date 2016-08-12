#include "PointLocatorRegularOrthogonal.h"
#include "PointLocatorRegularOrthogonalData.h"
#include "MooseError.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

PointLocatorRegularOrthogonal::PointLocatorRegularOrthogonal(const MeshBase & mesh,
                                                             const PointLocatorBase * master) :
    PointLocatorBase(mesh, master),
    _out_of_mesh_mode(false),
    _data(nullptr)
{
  if (master != nullptr)
  {
    const PointLocatorRegularOrthogonal * my_master =
      cast_ptr<const PointLocatorRegularOrthogonal *>(master);

    if (!my_master->initialized())
      mooseError("Linking a slave point locator to an uninitialized master is not allowed.");

    _data = my_master->_data;
    _initialized = true;
  }
}

PointLocatorRegularOrthogonal::~PointLocatorRegularOrthogonal()
{
  this->clear ();
}

void
PointLocatorRegularOrthogonal::clear()
{
  // only delete the root element table when we are the master
  if (this->_data != nullptr)
  {
    if (this->_master == nullptr)
      // we own the root element table
      delete this->_data;
    else
      // someone else owns and therefore deletes the root element table
      this->_data = nullptr;
  }

  _initialized = false;
}

void
PointLocatorRegularOrthogonal::init()
{
  mooseError("PointLocatorRegularOrthogonal needs to be explicitly initialized with mesh meta data.");
}

void
PointLocatorRegularOrthogonal::init(const std::vector<unsigned int> & cell_count,
                                    const Point & min_corner,
                                    const Point & max_corner)
{
  // initialize root element table on the master point locator only
  if (this->_master == nullptr)
    _data = new PointLocatorRegularOrthogonalData(cell_count, min_corner, max_corner, _mesh);

  _initialized = true;
}

const Elem *
PointLocatorRegularOrthogonal::operator() (const Point & p,
                                           const std::set<subdomain_id_type> * /* allowed_subdomains */) const
{
  const Elem * el = _data->rootElement(p);
  mooseAssert(el != nullptr, "No element found at p");

  // element has no active children
  if (el->active())
    return el;

  // TODO: bisect the children
  mooseError("AMR is not yet supported for use with the MyTRIMMesh and its associated point locator.");
}

void
PointLocatorRegularOrthogonal::operator() (const Point & p,
                                           std::set<const Elem *> & candidate_elements,
                                           const std::set<subdomain_id_type> * /* allowed_subdomains */) const
{
  // return just one match, the exact one TODO: return all neighbors if we are fuzzily on a face/edge/node
  candidate_elements.insert((*this)(p));
  return;
}
