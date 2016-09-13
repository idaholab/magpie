#ifndef MAGPIEPARALLEL_H
#define MAGPIEPARALLEL_H

#include "Moose.h"
#include "libmesh/parallel.h"
#include <vector>

// Let's not allow non-MPI builds for now. In the future a simple workaround would be a plain bufer copy.
#ifndef LIBMESH_HAVE_MPI
#error Magpie needs libMesh compiled with MPI.
#endif

namespace MagpieUtils
{
/**
 * Specialized replacement for the libMesh packed range communication that allows unlimited
 * string buffer sizes.
 */
void allgatherStringBuffers(const Parallel::Communicator & comm, std::string & send_buffer, std::vector<std::string> & recv_buffers);

}

#endif //MAGPIEPARALLEL_H
