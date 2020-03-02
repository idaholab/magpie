/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "LAMMPSFileRunner.h"
#include "MooseUtils.h"
#include "Function.h"
#include <sstream>
#include <cstring>
#include <cctype>

registerMooseObject("MagpieApp", LAMMPSFileRunner);

InputParameters
LAMMPSFileRunner::validParams()
{
  InputParameters params = MDRunBase::validParams();
  params.addClassDescription("Reads a LAMMPS dump file or sequence of dump files and "
                             "maps LAMMPS particles to MOOSE FEM mesh.");
  params.addParam<bool>("time_sequence",
                        false,
                        "If true, a sequence of dump files is read, "
                        "else a single dump file is read.");
  params.addRequiredParam<FileName>("lammps_file",
                                    "Name of LAMMPS dump file or file base of LAMMPS file for "
                                    "a sequence.");
  params.addRequiredParam<std::vector<unsigned int>>(
      "xyz_columns", "Column ids of the x, y, and z coordinates of particles.");
  std::vector<unsigned int> empty = {};
  params.addParam<std::vector<unsigned int>>(
      "property_columns", empty, "Column ids of the properties.");
  params.addParam<FunctionName>("time_conversion",
                                "A conversion from FEM simulation time to MD time stamps.");
  params.addClassDescription("Allows coupling FEM calculations to LAMMPS dump files.");
  return params;
}

LAMMPSFileRunner::LAMMPSFileRunner(const InputParameters & parameters)
  : MDRunBase(parameters),
    _time_sequence(getParam<bool>("time_sequence")),
    _lammps_file(getParam<FileName>("lammps_file")),
    _pos_columns(getParam<std::vector<unsigned int>>("xyz_columns")),
    _prop_columns(getParam<std::vector<unsigned int>>("property_columns")),
    _time_conversion(nullptr)
{
  ;
  if (_pos_columns.size() != _dim)
    mooseError(
        "Dimensionality is ", _dim, " but xyz_columns has ", _pos_columns.size(), " entries.");

  if (isParamValid("time_conversion"))
    _time_conversion = &getFunction("time_conversion");

  if (_properties.size() != _prop_columns.size())
    mooseError("properties array must have same length as property_columns");
}

void
LAMMPSFileRunner::initialSetup()
{
  // call MDRunBase initialSetup
  MDRunBase::initialSetup();

  if (_time_sequence)
    return;

  readLAMMPSFile(_lammps_file);
  updateKDTree();
}

void
LAMMPSFileRunner::updateParticleList()
{
  // read LAMMPS file if a time sequence is read
  if (_time_sequence)
  {
    // map to MD time
    Real md_time = _t;
    if (_time_conversion)
      md_time = _time_conversion->value(_t, Point());

    // find files w/ timestamps "bracketing" md_time
    std::pair<std::string, std::string> files;
    std::pair<Real, Real> timestamps;
    findBracketingLAMMPSFiles(md_time, files, timestamps);
    readLAMMPSFileHistory(files, timestamps, md_time);
    updateKDTree();
  }

  // update the mapping to mesh if mesh has changed or ...
  mapMDParticles();

  // update the granular candidates as well if necessary
  if (_granular)
    updateElementGranularVolumes();
}

void
LAMMPSFileRunner::findBracketingLAMMPSFiles(Real md_time,
                                            std::pair<std::string, std::string> & filenames,
                                            std::pair<Real, Real> & timestamps)
{
  /// the return values
  std::string before;
  std::string after;

  /// extract path and file_base from provided lammps file entry
  std::pair<std::string, std::string> p = MooseUtils::splitFileName(_lammps_file);

  /// get all files in the directory
  std::list<std::string> files = MooseUtils::getFilesInDirs({p.first});

  /// loop over all files, retain the ones that begin with base and extract their timestamps
  std::vector<Real> times;
  times.reserve(files.size());

  unsigned int lammps_files_found = 0;
  Real closest_after_time = std::numeric_limits<Real>::max();
  Real closest_before_time = -std::numeric_limits<Real>::max();
  for (auto & file : files)
  {
    std::pair<std::string, std::string> f = MooseUtils::splitFileName(file);
    if (f.second.compare(0, p.second.size(), p.second))
      continue;

    // split the file at dots and check the second item ==> timestamp
    Real timestamp;
    std::vector<std::string> elements;
    MooseUtils::tokenize(f.second, elements, 1, ".");
    if (elements.size() < 2)
      mooseError("Tokenizing filename ",
                 f.second,
                 " failed. LAMMPS filename must be base.<t>.xyz where <t> is the timestamp.");
    std::stringstream sstr(elements[1]);

    if (!isTimestamp(elements[1]))
      continue;
    sstr >> timestamp;

    // if (sstr >> timestamp || sstr.eof())
    //  if (sstr.fail())
    //    continue;

    // increase the counter & check if this file is a candidate for before/after
    ++lammps_files_found;
    if (timestamp > md_time && timestamp < closest_after_time)
    {
      closest_after_time = timestamp;
      after = file;
    }
    else if (timestamp <= md_time && timestamp > closest_before_time)
    {
      closest_before_time = timestamp;
      before = file;
    }
  }

  if (lammps_files_found == 0)
    mooseError("No LAMMPS file with name base ", p.second, " in folder ", p.first);

  // guard against the case that only before or after were found => set the two equal
  if (before.size() == 0)
    before = after;
  else if (after.size() == 0)
    after = before;

  filenames = std::pair<std::string, std::string>(before, after);
  timestamps = std::pair<Real, Real>(closest_before_time, closest_after_time);
}

void
LAMMPSFileRunner::readLAMMPSFile(FileName filename)
{
  // check if this file is readable
  MooseUtils::checkFileReadable(filename);

  // read LAMMPS file
  std::ifstream file(filename);

  // skip first three lines
  for (unsigned int j = 0; j < 3; ++j)
    file.ignore(10000, '\n');

  // get number of particles
  file >> _n_particles;

  // skip another five lines
  for (unsigned int j = 0; j < 6; ++j)
    file.ignore(10000, '\n');

  // attach the right _md_particles.pos
  std::vector<Point> & position = _md_particles.pos;

  // zero out the position & id vector
  position.clear();
  _md_particles.id.clear();
  _md_particles.properties.clear();

  // size the md particle vector considering not all processors have all particles
  position.reserve(std::ceil(2.0 * _n_particles / _nproc));
  _md_particles.id.reserve(std::ceil(2.0 * _n_particles / _nproc));
  _md_particles.properties.reserve(std::ceil(2.0 * _n_particles / _nproc));

  // read particle entries
  _n_local_particles = 0;
  for (unsigned int j = 0; j < _n_particles; ++j)
  {
    // check if this particle is in the local bounding box
    std::string line;
    std::getline(file, line);
    std::vector<std::string> elements;
    MooseUtils::tokenize<std::string>(line, elements, 1, " ");

#if DEBUG
    // check that enough columns exist in debug
    // find largest requested column
    unsigned int largest_column = 0;
    for (unsigned int i = 0; i < _dim; ++i)
      if (_pos_columns[i] > largest_column)
        largest_column = _pos_columns[i];

    for (unsigned int i = 0; i < _prop_columns.size(); ++i)
      if (_prop_columns[i] > largest_column)
        largest_column = _prop_columns[i];

    if (largest_column >= elements.size())
      mooseError("Error reading ", filename, " on line ", line);
#endif

    // read location
    Point pos;
    for (unsigned int i = 0; i < _dim; ++i)
    {
      std::stringstream sstr(elements[_pos_columns[i]]);
      sstr >> pos(i);
    }

    // read properties
    std::vector<Real> props(_md_particles._prop_size);
    for (unsigned int i = 0; i < _prop_columns.size(); ++i)
    {
      unsigned int prop_id = _properties.get(i);
      std::stringstream sstr(elements[_prop_columns[i]]);
      sstr >> props[propIndex(prop_id)];
    }

    // check if this particle is in this processor BBox
    if (!_bbox[processor_id()].contains_point(pos))
      continue;

    ++_n_local_particles;
    position.push_back(pos);
    _md_particles.id.push_back(j);
    _md_particles.properties.push_back(props);
  }
  file.close();

  // size back
  position.shrink_to_fit();
  _md_particles.id.shrink_to_fit();
}

bool
LAMMPSFileRunner::isTimestamp(std::string ts_candidate) const
{
  for (auto & s : ts_candidate)
    if (!isdigit(s))
      return false;
  return true;
}

void
LAMMPSFileRunner::readLAMMPSFileHistory(std::pair<FileName, FileName> filenames,
                                        std::pair<Real, Real> timestamps,
                                        Real current_time)
{
  // check if this file is readable
  MooseUtils::checkFileReadable(filenames.first);
  MooseUtils::checkFileReadable(filenames.second);

  // open LAMMPS file
  std::ifstream file_before(filenames.first);
  std::ifstream file_after(filenames.second);

  // compute weights for particle
  Real weight = (current_time - timestamps.first) / (timestamps.second - timestamps.first);

  // skip first three lines
  for (unsigned int j = 0; j < 3; ++j)
  {
    file_before.ignore(10000, '\n');
    file_after.ignore(10000, '\n');
  }

  // get number of particles
  unsigned int np;
  file_before >> _n_particles;
  file_after >> np;
  if (np != _n_particles)
    mooseError("Files ",
               filenames.first,
               " : ",
               _n_particles,
               " and ",
               filenames.second,
               " : ",
               np,
               " have different number of particles.");

  // skip another five lines
  for (unsigned int j = 0; j < 6; ++j)
  {
    file_before.ignore(10000, '\n');
    file_after.ignore(10000, '\n');
  }

  // attach the right _md_particles.pos
  std::vector<Point> & position = _md_particles.pos;

  // zero out the position & id vector
  position.clear();
  _md_particles.id.clear();
  _md_particles.properties.clear();

  // size the md particle vector considering not all processors have all particles
  position.reserve(std::ceil(2.0 * _n_particles / _nproc));
  _md_particles.id.reserve(std::ceil(2.0 * _n_particles / _nproc));
  _md_particles.properties.reserve(std::ceil(2.0 * _n_particles / _nproc));

  // read particle entries
  _n_local_particles = 0;
  for (unsigned int j = 0; j < _n_particles; ++j)
  {
    // check if this particle is in the local bounding box
    // either in before or after state

    // get the MD particle from the before file
    std::vector<std::string> elements;
    std::string line_before;
    std::getline(file_before, line_before);
    MooseUtils::tokenize<std::string>(line_before, elements, 1, " ");

#if DEBUG
    // check that enough columns exist in debug
    // find largest requested column
    unsigned int largest_column = 0;
    for (unsigned int i = 0; i < _dim; ++i)
      if (_pos_columns[i] > largest_column)
        largest_column = _pos_columns[i];

    for (unsigned int i = 0; i < _prop_columns.size(); ++i)
      if (_prop_columns[i] > largest_column)
        largest_column = _prop_columns[i];

    if (largest_column >= elements.size())
      mooseError("Error reading ", filenames.first, " on line ", line_before);
#endif

    Point pos_before;
    for (unsigned int i = 0; i < _dim; ++i)
    {
      std::stringstream sstr(elements[_pos_columns[i]]);
      sstr >> pos_before(i);
    }

    // read properties from before file
    std::vector<Real> props_before(_md_particles._prop_size);
    for (unsigned int i = 0; i < _prop_columns.size(); ++i)
    {
      unsigned int prop_id = _properties.get(i);
      std::stringstream sstr(elements[_prop_columns[i]]);
      sstr >> props_before[propIndex(prop_id)];
    }

    // get the MD particle from the after file
    elements.clear();
    std::string line_after;
    std::getline(file_after, line_after);

    MooseUtils::tokenize<std::string>(line_after, elements, 1, " ");

    // we have determined largest_column in DEBUG so we can use it in mooseAssert
    mooseAssert(largest_column < elements.size(),
                "Error reading " << filenames.second << " on line " << line_after);

    Point pos_after;
    for (unsigned int i = 0; i < _dim; ++i)
    {
      std::stringstream sstr(elements[_pos_columns[i]]);
      sstr >> pos_after(i);
    }

    // read properties from before file
    std::vector<Real> props_after(_md_particles._prop_size);
    for (unsigned int i = 0; i < _prop_columns.size(); ++i)
    {
      unsigned int prop_id = _properties.get(i);
      std::stringstream sstr(elements[_prop_columns[i]]);
      sstr >> props_after[propIndex(prop_id)];
    }

    // check if this particle is in this processor BBox
    if (!_bbox[processor_id()].contains_point(pos_after) &&
        !_bbox[processor_id()].contains_point(pos_before))
      continue;

    ++_n_local_particles;
    position.push_back((1 - weight) * pos_before + weight * pos_after);
    _md_particles.id.push_back(j);
    _md_particles.properties.push_back(std::vector<Real>(_md_particles._prop_size));
    for (unsigned int i = 0; i < _md_particles._prop_size; ++i)
      _md_particles.properties.back()[i] = (1 - weight) * props_before[i] + weight * props_after[i];
  }
  file_before.close();
  file_after.close();

  // size back
  position.shrink_to_fit();
  _md_particles.id.shrink_to_fit();
}
