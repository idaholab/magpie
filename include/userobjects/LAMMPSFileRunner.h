/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "MDRunBase.h"

class LAMMPSFileRunner;
class Function;

template <>
InputParameters validParams<LAMMPSFileRunner>();

/**
 * Reads lammps dump files to emulate MD simulation
 */
class LAMMPSFileRunner : public MDRunBase
{
public:
  LAMMPSFileRunner(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual void finalize() override {}
  virtual void initialSetup() override;

  virtual void updateParticleList() override;

protected:
  /// reads a LAMMPS file
  void readLAMMPSFile(FileName filename);

  /// reads two LAMMPS files and interpolates times
  void readLAMMPSFileHistory(std::pair<FileName, FileName> filenames,
                             std::pair<Real, Real> timestamps,
                             Real current_time);

  /// helper function that finds two files, one right before and one right after md_time
  void findBracketingLAMMPSFiles(Real md_time,
                                 std::pair<std::string, std::string> & filenames,
                                 std::pair<Real, Real> & timestamps);

  /// whether a sequence of input files or a single input file is read
  bool _time_sequence;

  /// name of LAMMPS file or file base if _time_sequence == true
  FileName _lammps_file;

  /// column of x, y, z coordinate in LAMMPS files
  std::vector<unsigned int> _pos_columns;

  /// Conversion from FEM time to MD time_stamp
  const Function * _time_conversion;
};

