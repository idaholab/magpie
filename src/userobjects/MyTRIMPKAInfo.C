#include "MyTRIMPKAInfo.h"
#include "MyTRIMRasterizer.h"
#include "MagpieParallel.h"
#include "mytrim/ion.h"

template<>
InputParameters validParams<MyTRIMPKAInfo>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<UserObjectName>("rasterizer", "Name of the MyTRIMRasterizer userobject to pull data from");
  params.addParam<unsigned int>("number_channels", 50, "Number of energy channels");
  params.addParam<Real>("channel_width", 5.0e6, "Energy channel width in eV");
  return params;
}

MyTRIMPKAInfo::MyTRIMPKAInfo(const InputParameters & params) :
    GeneralUserObject(params),
    _rasterizer(getUserObject<MyTRIMRasterizer>("rasterizer")),
    _nchannels(getParam<unsigned int>("number_channels")),
    _deltaE(getParam<Real>("channel_width"))
{
}

void
MyTRIMPKAInfo::initialize()
{
  // reset the statistics
  // _npka_per_element.clear();
  _npka_per_zaid.clear();
  _npka_per_mass.clear();
  _npka_per_energy.clear();
  _npka_per_energy.resize(_nchannels);
}

void
MyTRIMPKAInfo::threadJoin(const UserObject &y)
{
  const MyTRIMPKAInfo & uo = static_cast<const MyTRIMPKAInfo &>(y);

  for (unsigned int j = 0; j < uo._npka_per_energy.size(); ++j)
    _npka_per_energy[j] += uo._npka_per_energy[j];

  for (auto & entry : uo._npka_per_mass)
    if (_npka_per_mass.find(entry.first) != _npka_per_mass.end())
      _npka_per_mass[entry.first] += entry.second;
    else
      _npka_per_mass[entry.first] = entry.second;

  for (auto & entry : uo._npka_per_zaid)
    if (_npka_per_zaid.find(entry.first) != _npka_per_zaid.end())
      _npka_per_zaid[entry.first] += entry.second;
    else
      _npka_per_zaid[entry.first] = entry.second;

  //for (auto & entry : uo._npka_per_element)
  //  if (_npka_per_element.find(entry.first) != _npka_per_element.end())
  //    _npka_per_element[entry.first] += entry.second;
  //  else
  //    _npka_per_element[entry.first] = entry.second;
}

void
MyTRIMPKAInfo::execute()
{
  const std::vector<MyTRIM_NS::IonBase> & pka_list = _rasterizer.getPKAList();
  for (auto & entry : pka_list)
  {
    unsigned int M = round(entry._m);
    unsigned int Z = entry._Z;
    Real energy = entry._E;

    unsigned int channel = energy / _deltaE;
    if (channel >= _nchannels)
      channel = _nchannels - 1;
    _npka_per_energy[channel] += 1;

    if (_npka_per_mass.find(M) != _npka_per_mass.end())
      _npka_per_mass[M] += 1;
    else
      _npka_per_mass[M] = 1;

    // zaids
    unsigned int zaid = Z * 1000 + M;
    if (_npka_per_zaid.find(zaid) != _npka_per_zaid.end())
      _npka_per_zaid[zaid] += 1;
    else
      _npka_per_zaid[zaid] = 1;
  }
}

void
MyTRIMPKAInfo::finalize()
{
  // for single processor runs we do not need to do anything here
  if (_app.n_processors() > 1)
  {
    // create send buffer
    std::string send_buffer;

    // create byte buffers for the streams received from all processors
    std::vector<std::string> recv_buffers;

    // pack the complex datastructures into the string stream
    serialize(send_buffer);

    // broadcast serialized data to and receive from all processors
    MagpieUtils::allgatherStringBuffers(_communicator, send_buffer, recv_buffers);

    // unpack the received data and merge it into the local data structures
    deserialize(recv_buffers);
  }

  // print statistics on processor 0
  if (processor_id() == 0)
    printPKAStats();
}

void
MyTRIMPKAInfo::serialize(std::string & serialized_buffer)
{
  // stream for serializing the _npka_per<mass,zaid,energy> data structures to byte streams
  std::ostringstream oss;
  //dataStore(oss, _npka_per_element, this);
  dataStore(oss, _npka_per_mass, this);
  dataStore(oss, _npka_per_zaid, this);
  dataStore(oss, _npka_per_energy, this);

  // Populate the passed in string pointer with the string stream's buffer contents
  serialized_buffer.assign(oss.str());
}

void
MyTRIMPKAInfo::deserialize(std::vector<std::string> & serialized_buffers)
{
  mooseAssert(serialized_buffers.size() == _app.n_processors(), "Unexpected size of serialized_buffers: " << serialized_buffers.size());

  // The input string stream used for deserialization
  std::istringstream iss;

  // Loop over all datastructures for all procerssors to perfrom the gather operation
  for (unsigned int rank = 0; rank < serialized_buffers.size(); ++rank)
  {
    // skip the current processor (its data is already in the structures)
    if (rank == processor_id())
      continue;

    // populate the stream with a new buffer and reset stream state
    iss.str(serialized_buffers[rank]);
    iss.clear();

    // Load the communicated data into temporary structures
    // std::map<dof_id_type, unsigned int> other_npka_per_element;
    // dataLoad(iss, other_npka_per_element, this);
    std::map<unsigned int, unsigned int> other_npka_per_mass;
    dataLoad(iss, other_npka_per_mass, this);
    std::map<unsigned int, unsigned int> other_npka_per_zaid;
    dataLoad(iss, other_npka_per_zaid, this);
    std::vector<unsigned int> other_npka_per_energy;
    dataLoad(iss, other_npka_per_energy, this);

    // merge the data in with the current processor's data
    // _npka_per_element.insert(other_npka_per_element.begin(), other_npka_per_element.end());
    _npka_per_mass.insert(other_npka_per_mass.begin(), other_npka_per_mass.end());
    _npka_per_zaid.insert(other_npka_per_zaid.begin(), other_npka_per_zaid.end());

    // merging the PKA lists
    _npka_per_energy.insert(_npka_per_energy.begin(), other_npka_per_energy.begin(), other_npka_per_energy.end());
  }
}

void
MyTRIMPKAInfo::printPKAStats() const
{
  // Print by mass
  _console << "Printing PKA information:" << std::endl << "1. PKAs by mass" << std::endl;
  for (auto & mass_entry : _npka_per_mass)
    _console << "M: " << mass_entry.first << " Number: " << mass_entry.second << std::endl;

  // print by ZAID
  _console << std::endl << "2. PKAs by ZAID" << std::endl;
  for (auto & zaid_entry : _npka_per_zaid)
    _console << "ZAID: " << zaid_entry.first << " Number: " << zaid_entry.second << std::endl;

  // Print by energy
  _console << std::endl << "3. PKAs by energy" << std::endl;
  for (unsigned int j = 0; j < _npka_per_energy.size(); ++j)
    _console << "Energy range (MeV): " << std::scientific << std::setprecision(4) << j * _deltaE << " - " << (j + 1) * _deltaE << ": "<< _npka_per_energy[j] << std::endl;
}
