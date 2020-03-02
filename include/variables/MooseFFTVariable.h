/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MooseVariable.h"

class MooseFFTVariable : public MooseVariable
{
public:
  static InputParameters validParams();

  MooseFFTVariable(const InputParameters & parameters);

  virtual bool isNodal() const { return false; }

  /**
   * Clear out the dof indices.  We do this in case this variable is not going to be prepared at
   * all...
   */
  virtual void clearDofIndices() {}

  /**
   * Prepare the elemental degrees of freedom
   */
  virtual void prepare() {}

  /**
   * Prepare the neighbor element degrees of freedom
   */
  virtual void prepareNeighbor() { mooseError("Neighbor FFT variables are not supported yet."); }

  /**
   * Prepare a lower dimensional element's degrees of freedom
   */
  virtual void prepareLowerD()
  {
    mooseError("Lower dimensional FFT variables are not supported yet.");
  }

  virtual void prepareAux() {}

  virtual void reinitNode() {}
  virtual void reinitAux() {}
  virtual void reinitAuxNeighbor() { mooseError("Neighbor FFT variables are not supported yet."); }

  virtual void reinitNodes(const std::vector<dof_id_type> &)
  {
    mooseError("Nodal FFT variables are not supported.");
  }
  virtual void reinitNodesNeighbor(const std::vector<dof_id_type> &)
  {
    mooseError("Nodal FFT variables are not supported.");
  }

  /**
   * Field type of this variable
   */
  virtual Moose::VarFieldType fieldType() const { return Moose::VarFieldType::VAR_FIELD_STANDARD; }

  /**
   * @returns true if this is a vector-valued element, false otherwise.
   */
  virtual bool isVector() const { return false; };

  /**
   * Is this variable defined at nodes
   * @return true if it the variable is defined at nodes, otherwise false
   */
  virtual bool isNodalDefined() const { return false; }

  virtual const dof_id_type & nodalDofIndex() const
  {
    mooseError("Nodal FFT variables are not supported.");
  }

  virtual const dof_id_type & nodalDofIndexNeighbor() const
  {
    mooseError("Nodal FFT variables are not supported.");
  }

  /**
   * Current element this variable is evaluated at
   */
  virtual const Elem * const & currentElem() const
  {
    mooseError("Current element access not supported.");
  }

  /**
   * The subdomains the variable is active on
   */
  virtual const std::set<SubdomainID> & activeSubdomains() const { return _all_subdomains; }
  /**
   * Is the variable active on the subdomain?
   * @param subdomain The subdomain id in question
   * @return true if active on subdomain, false otherwise
   */
  virtual bool activeOnSubdomain(SubdomainID /*subdomain*/) const { return true; }

  /**
   * Prepare the initial condition
   */
  virtual void prepareIC() {}

  /**
   * Compute values at interior quadrature points
   */
  virtual void computeElemValues()
  { // meat!
  }

  /**
   * Compute values at facial quadrature points
   */
  virtual void computeElemValuesFace() {}
  /**
   * Compute values at facial quadrature points for the neighbor
   */
  virtual void computeNeighborValuesFace() {}
  /**
   * Compute values at quadrature points for the neighbor
   */
  virtual void computeNeighborValues() {}
  /**
   * compute values at quadrature points on the lower dimensional element
   */
  virtual void computeLowerDValues() {}
  /**
   * Compute nodal values of this variable in the neighbor
   */
  virtual void computeNodalNeighborValues() {}
  /**
   * Compute nodal values of this variable
   */
  virtual void computeNodalValues() {}

  virtual void getDofIndices(const Elem * /*elem*/,
                             std::vector<dof_id_type> & /*dof_indices*/) const
  {
  }

  virtual unsigned int numberOfDofsNeighbor() { return 0; }

  virtual void insert(NumericVector<Number> & /*residual*/) {}
  virtual void add(NumericVector<Number> & /*residual*/) {}

protected:
  std::vector<dof_id_type> _no_dofs;
  std::set<SubdomainID> _all_subdomains;
};
