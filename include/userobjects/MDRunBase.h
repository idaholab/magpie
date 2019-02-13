/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef MDRUNBASE_H
#define MDRUNBASE_H

#include "GeneralUserObject.h"
#include "KDTree.h"
#include "libmesh/bounding_box.h"

class MDRunBase;
class MooseMesh;

template <>
InputParameters validParams<MDRunBase>();

/**
 * Base class for molecular dynamics runs in Magpie
 */
class MDRunBase : public GeneralUserObject
{
public:
  MDRunBase(const InputParameters & parameters);

  void initialSetup() override;
  void timestepSetup() override;

  class MDParticles
  {
  public:
    /// particle's position
    std::vector<Point> pos;

    /// particle's velocity
    std::vector<Point> vel;

    /// the id of particle in the MD calculation
    std::vector<unsigned int> id;

    /// the mesh element the particles are in
    std::vector<unique_id_type> elem_id;
  };

  /// access to the MDParticles
  const MDParticles & particles() const { return _md_particles; }

  /// access to the element to particle map
  const std::vector<unsigned int> elemParticles(unique_id_type elem_id) const;

  /// List of quantities to get from MD simulation
  MultiMooseEnum getMDQuantities() const;

protected:
  /// call back function to update the particle list
  virtual void updateParticleList() = 0;

  /// updates the KDTree object
  void updateKDTree();

  /// map MDParticles to elements
  void mapMDParticles();

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;

  /// dimension of the mesh
  const unsigned int _nproc;

  /// stores bounding boxes of all other processors
  std::vector<BoundingBox> _bbox;

  /// total number of particles
  unsigned int _n_particles;

  /// number of local particles
  unsigned int _n_local_particles;

  /// stores the location of
  MDParticles _md_particles;

  /// a map from elem pointer to particles in this element
  std::map<unique_id_type, std::vector<unsigned int>> _elem_particles;

  /// A KDTree object to handle md_particles
  std::unique_ptr<KDTree> _kd_tree;
};

#endif // MDRUNBASE_H
