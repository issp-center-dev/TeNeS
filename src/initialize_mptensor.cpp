#include "tenes.hpp"

#include <mptensor/mptensor.hpp>

namespace tenes{
void initialize_mptensor(){
#ifndef _NO_MPI
  mptensor::scalapack::BlacsGrid::init();
#endif
}
}
