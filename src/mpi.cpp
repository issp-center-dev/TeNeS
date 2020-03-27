#include <cstring>

#include "mpi.hpp"

#ifdef _NO_MPI

int MPI_Init(int*, char***) { return 0; }
int MPI_Barrier(MPI_Comm) { return 0; }
int MPI_Finalize() { return 0; }

int MPI_Comm_size(MPI_Comm, int* size) {
  *size = 1;
  return 0;
}
int MPI_Comm_rank(MPI_Comm, int* rank) {
  *rank = 0;
  return 0;
}

int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm) { return 0; }
#endif

namespace tenes{

template <>
MPI_Datatype get_MPI_Datatype<int>() { return MPI_INT; }
template <>
MPI_Datatype get_MPI_Datatype<double>() { return MPI_DOUBLE; }
template <>
MPI_Datatype get_MPI_Datatype<bool>() { return MPI_INT; }

int bcast(bool &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  int v = val?0:1;
  ret=MPI_Bcast(&v, 1, MPI_INT, root, comm);
  val = v==0;
#endif
  return ret;
}

int bcast(std::string &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  int mpirank;
  MPI_Comm_rank(comm, &mpirank);
  int sz = val.size()+1;
  bcast(sz, root, comm);
  char* buf = new char[sz];
  if(mpirank==root){
    std::strncpy(buf, val.c_str(), sz);
  }
  ret = MPI_Bcast(buf, sz, MPI_CHAR, root, comm);
  if(mpirank!=root){
    val = buf;
  }
  delete buf;
#endif
  return ret;
}

int bcast(std::vector<std::string> &val, int root, MPI_Comm comm){
  int ret=0;
#ifndef _NO_MPI
  int mpirank;
  MPI_Comm_rank(comm, &mpirank);
  int sz = val.size();
  bcast(sz, root, comm);
  if(mpirank!=root){
    val.resize(sz);
  }
  for(int i=0; i<sz; ++i){
    int r=bcast(val[i], root, comm);
    if(r!=0){
      ret = r;
    }
  }
#endif
  return ret;
}

} // end of namespace tenes
