/* TeNeS - Massively parallel tensor network solver /
/ Copyright (C) 2019- The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify /
/ it under the terms of the GNU General Public License as published by /
/ the Free Software Foundation, either version 3 of the License, or /
/ (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, /
/ but WITHOUT ANY WARRANTY; without even the implied warranty of /
/ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
/ GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License /
/ along with this program. If not, see http://www.gnu.org/licenses/. */

#include "iTPS.hpp"
#include <cmath>

namespace tenes {
namespace itps {
template <class tensor>
void iTPS<tensor>::make_LGS() {
  using mptensor::Shape, mptensor::Axes, mptensor::Index;

  // initialize iTPS for LGS //  
  std::vector<tensor> LGS;
  LGS.clear(); 
  for (int i = 0; i < N_UNIT; ++i) {
    const auto pdim = lattice.physical_dims[i];
    if (pdim > 3){
      std::stringstream ss;
      ss << "ERROR: LGS support only S=1/2 and S=1" << std::endl;
      throw std::runtime_error(ss.str());
    }
    LGS.push_back(tensor(comm, Shape(1, 1, 1, 1, pdim)));
  }

  // assuming all pdim s are identical
  int pdim;
  pdim = lattice.physical_dims[0];
  std::vector<double> real_part,img_part;
  double factor;
  factor = (std::sqrt(3.0) - 1.0) * 0.5;
  Index index;

  if (pdim == 2){
    // S = 1/2
    real_part.push_back(1.0);
    real_part.push_back(factor);

    img_part.push_back(0.0);
    img_part.push_back(factor);
      
    for (int i = 0; i < lattice.N_UNIT; ++i){
      for (int n = 0; n < LGS[i].local_size(); ++n){
	index = LGS[i].global_index(n);
	// in the present case, local_size == pdim
	auto v = std::complex<double>(real_part[index[4]], img_part[index[4]]);
	LGS[i].set_value(index, to_tensor_type(v));
      }
    }

    // projectors
    tensor Q = tensor(comm, Shape(2, 2, 2, pdim, pdim));

    // Identity
    Q.set_value(Index(0 ,0 ,0, 0, 0), to_tensor_type(1.0));
    Q.set_value(Index(0 ,0 ,0, 1, 1), to_tensor_type(1.0));
    
    // sigma_x
    Q.set_value(Index(0 ,1 ,1, 0, 1), to_tensor_type(1.0));
    Q.set_value(Index(0 ,1 ,1, 1, 0), to_tensor_type(1.0));
    
    // sigma_y
    Q.set_value(Index(1 ,0 ,1, 0, 1), to_tensor_type(std::complex<double>(0.0, -1.0)));
    Q.set_value(Index(1 ,0 ,1, 1, 0), to_tensor_type(std::complex<double>(0.0, 1.0)));
    
    // sigma_z
    Q.set_value(Index(1 ,1 ,0, 0, 0), to_tensor_type(1.0));
    Q.set_value(Index(1 ,1 ,0, 1, 1), to_tensor_type(-1.0));

    if (lattice.initial_dirs[0][0] == lattice.initial_dirs[1][0]){
      // ferro LGS
      std::cout<< "Creating Ferro LGS for S = 1/2"<<std::endl;
      for (int i = 0; i < lattice.N_UNIT; ++i){
	int ix = i % LX;
	int iy = i / LX;
	if ((ix + iy) % 2 == 0){
	  // sublattice 1
	  LGS[i] = mptensor::reshape(mptensor::tensordot(Q, LGS[i], Axes(4), Axes(4)).transpose(Axes(1, 0, 2, 3, 4, 5, 6, 7)), Shape(2, 1, 2, 2, 2));
	} else{
	  // sublattice 2
	  LGS[i] = mptensor::reshape(mptensor::tensordot(Q, LGS[i], Axes(4), Axes(4)).transpose(Axes(0, 2, 1, 3, 4, 5, 6, 7)), Shape(2, 2, 2, 1, 2));
	}
      }
    } else {
      // antiferro LGS
      std::cout<< "Creating Antiferro LGS for S = 1/2"<<std::endl;
      tensor v = tensor(comm, Shape(2, 2));
      tensor vd = tensor(comm, Shape(2, 2));

      v.set_value(Index(0,1), to_tensor_type(std::complex<double>(0.0, 1.0)));
      v.set_value(Index(1,0), to_tensor_type(std::complex<double>(1.0, 0.0)));

      vd.set_value(Index(0,1), to_tensor_type(std::complex<double>(1.0, 0.0)));
      vd.set_value(Index(1,0), to_tensor_type(std::complex<double>(0.0, -1.0)));
      
      for (int i = 0; i < lattice.N_UNIT; ++i){
	int ix = i % LX;
	int iy = i / LX;
	if ((ix -  iy + LX) % 4 == 0){
	  // sublattice 1
	  LGS[i] = mptensor::reshape(mptensor::tensordot(
							 v, mptensor::tensordot(
									     vd, mptensor::tensordot(
												 v, mptensor::tensordot(
														     Q, LGS[i], Axes(4), Axes(4)), 
												 Axes(1), Axes(2)),
									     Axes(0), Axes(1)),
							 Axes(0), Axes(2)),
				     Shape(2, 1, 2, 2, 2));
	  
	} else if (((ix - iy + LX) % 4) % 2 == 1){
	  // sublattice 2 or 4
	  LGS[i] = mptensor::reshape(mptensor::tensordot(Q, LGS[i], Axes(4), Axes(4)).transpose(Axes(0, 2, 1, 3, 4, 5, 6, 7)), Shape(2, 2, 2, 1, 2));
	} else if ((ix - iy + LX) % 4 == 2){
	  // sublattice 3
	  LGS[i] = mptensor::reshape(mptensor::tensordot(
							 vd, mptensor::tensordot(
									     v, mptensor::tensordot(
												 vd, mptensor::tensordot(
														     Q, LGS[i], Axes(4), Axes(4)), 
												 Axes(1), Axes(2)),
									     Axes(1), Axes(1)),
							 Axes(0), Axes(2)),
				     Shape(2, 1, 2, 2, 2));
	}
      }
    }
    
  } else {
    // S = 1
    real_part.push_back(1.0);
    real_part.push_back(std::sqrt(2.0) * factor);
    real_part.push_back(0.0);

    img_part.push_back(0.0);
    img_part.push_back(std::sqrt(2.0) * factor);
    img_part.push_back(2.0 * factor * factor);
      
    for (int i = 0; i < lattice.N_UNIT; ++i){
      for (int n = 0; n < LGS[i].local_size(); ++n){
	index = LGS[i].global_index(n);
	// in the present case, local_size == pdim
	auto v = std::complex<double>(real_part[index[4]], img_part[index[4]]);
	LGS[i].set_value(index, to_tensor_type(v));
      }
    }
    // projectors
    tensor Q = tensor(comm, Shape(2, 2, 2, pdim, pdim));

    // Identity
    Q.set_value(Index(0 ,0 ,0, 0, 0), to_tensor_type(1.0));
    Q.set_value(Index(0 ,0 ,0, 1, 1), to_tensor_type(1.0));
    Q.set_value(Index(0 ,0 ,0, 2, 2), to_tensor_type(1.0));
    
    // sigma_x
    Q.set_value(Index(0 ,1 ,1, 0, 2), to_tensor_type(-1.0));
    Q.set_value(Index(0 ,1 ,1, 1, 1), to_tensor_type(-1.0));
    Q.set_value(Index(0 ,1 ,1, 2, 0), to_tensor_type(-1.0));
    
    // sigma_y
    Q.set_value(Index(1 ,0 ,1, 0, 2), to_tensor_type(1.0));
    Q.set_value(Index(1 ,0 ,1, 1, 1), to_tensor_type(-1.0));
    Q.set_value(Index(1 ,0 ,1, 2, 0), to_tensor_type(1.0));
    
    // sigma_z
    Q.set_value(Index(1 ,1 ,0, 0, 0), to_tensor_type(-1.0));
    Q.set_value(Index(1 ,1 ,0, 1, 1), to_tensor_type(1.0));
    Q.set_value(Index(1 ,1 ,0, 2, 2), to_tensor_type(-1.0));


    if (lattice.initial_dirs[0][0] == lattice.initial_dirs[1][0]){
      // ferro LGS
      std::cout<< "Creating Ferro LGS for S = 1"<<std::endl;

      for (int i = 0; i < lattice.N_UNIT; ++i){
	int ix = i % LX;
	int iy = i / LX;
	if ((ix + iy) % 2 == 0){
	  // sublattice 1
	  LGS[i] = mptensor::reshape(mptensor::tensordot(Q, LGS[i], Axes(4), Axes(4)).transpose(Axes(1, 0, 2, 3, 4, 5, 6, 7)), Shape(2, 1, 2, 2, 3));
	} else{
	  // sublattice 2
	  LGS[i] = mptensor::reshape(mptensor::tensordot(Q, LGS[i], Axes(4), Axes(4)).transpose(Axes(0, 2, 1, 3, 4, 5, 6, 7)), Shape(2, 2, 2, 1, 3));
	}
      }
    } else {
      // antiferro LGS
      std::cout<< "Creating Antiferro LGS for S = 1"<<std::endl;
      tensor x = tensor(comm, Shape(2, 2));

      x.set_value(Index(0,1), to_tensor_type(1.0));
      x.set_value(Index(1,0), to_tensor_type(1.0));
      
      for (int i = 0; i < lattice.N_UNIT; ++i){
	int ix = i % LX;
	int iy = i / LX;
	if ((ix + iy) % 2 == 0){
	  // sublattice 1
	  LGS[i] = mptensor::reshape(mptensor::tensordot(
							 x, mptensor::tensordot(
									     x, mptensor::tensordot(
												 x, mptensor::tensordot(
														     Q, LGS[i], Axes(4), Axes(4)), 
												 Axes(1), Axes(2)),
									     Axes(1), Axes(1)),
							 Axes(1), Axes(2)),
				     Shape(2, 1, 2, 2, 3));
	} else{
	  // sublattice 2
	  LGS[i] = mptensor::reshape(mptensor::tensordot(Q, LGS[i], Axes(4), Axes(4)).transpose(Axes(0, 2, 1, 3, 4, 5, 6, 7)), Shape(2, 2, 2, 1, 3));
	}
      }
    }
  }

  // copy to Tn

  Tn.clear();
  for (int i = 0; i < N_UNIT; ++i) {
    const auto pdim = lattice.physical_dims[i];
    const auto vdim = lattice.virtual_dims[i];

    Tn.push_back(tensor(comm, Shape(vdim[0], vdim[1], vdim[2], vdim[3], pdim)));
  }

  for (int i = 0; i < lattice.N_UNIT; ++i) {
    const auto pdim = lattice.physical_dims[i];
    const auto vdim = lattice.virtual_dims[i];
    

    for (int n = 0; n < Tn[i].local_size(); ++n) {
      index = Tn[i].global_index(n);
      if (index[0] < 2 && index[1] < 2 && index[2] < 2 && index[3] < 2) {
	typename tensor::value_type v;
	bool test;
	test = LGS[i].get_value(index, v);
	Tn[i].set_value(index, to_tensor_type(v));
      } else {
	Tn[i].set_value(index, to_tensor_type(0.0));
      }
    }
  }
  // output tensors
  
}
// template specialization
template class iTPS<real_tensor>;
template class iTPS<complex_tensor>;

}  // namespace itps
}  // namespace tenes
