#ifndef UNSIGNEDINTEGERPOINT_H
#define UNSIGNEDINTEGERPOINT_H

#include "libmesh/type_vector.h"

class UnsignedIntegerPoint : public TypeVector<unsigned int>
{
public:
  UnsignedIntegerPoint(const unsigned int x = 0,
                       const unsigned int y = 0,
                       const unsigned int z = 0) :
    TypeVector<unsigned int>(x, y, z) {}
};

#endif //UNSIGNEDINTEGERPOINT_H
