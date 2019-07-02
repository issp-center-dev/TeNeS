#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>
#include <complex>

#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "Parameters.hpp"
#include "Square_lattice_CTM.hpp"
#include "edge.hpp"
#include "hamiltonian.hpp"

using namespace mptensor;
typedef Tensor<scalapack::Matrix, double> ptensor;

void Initialize_Tensors(std::vector<ptensor> &Tn, Lattice lattice) {
  // ferro
  int D = Tn[0].shape()[0];

  std::vector<double> ran(D * D * D * D * 2);
  std::mt19937 gen(11);
  std::uniform_real_distribution<double> dist(-0.01, 0.01);
  for (int i = 0; i < D * D * D * D * 2; i++) {
    ran[i] = dist(gen);
  }
  int nr;

  Index index;
  for (int i = 0; i < lattice.N_UNIT; ++i) {
    for (int n = 0; n < Tn[i].local_size(); ++n) {
      index = Tn[i].global_index(n);
      if (index == Index(0, 0, 0, 0, 0)) {
        Tn[i].set_value(index, 1.0);
      } else {
        nr = index[0] + index[1] * D + index[2] * D * D + index[3] * D * D * D +
             index[4] * D * D * D * D;
        Tn[i].set_value(index, ran[nr]);
      }
    }
  }
}

template <class ptensor>
int tnsolve(MPI_Comm comm,
            PEPS_Parameters peps_parameters,
            Parameters local_parameters,
            Lattice lattice,
            Edges simple_edges,
            Edges full_edges,
            std::vector<ptensor> hams
    ){
  int mpisize, mpirank;

  /* MPI initialization */
  // MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &mpirank);
  MPI_Comm_size(comm, &mpisize);

  // for measure time
  double time_simple_update = 0.0;
  double time_full_update = 0.0;
  double time_env = 0.0;
  double time_obs = 0.0;
  double start_time;

  std::cout << std::setprecision(12);

  constexpr int nleg = 4;

  local_parameters.Bcast(comm);
  peps_parameters.Bcast(comm);
  lattice.Bcast(comm);

  // output debug or warning info only from process 0
  if (mpirank != 0) {
    peps_parameters.Debug_flag = false;
    peps_parameters.Warning_flag = false;
  }

  if (mpirank == 0) {
    // folder check
    struct stat status;
    if (stat("output_data", &status) != 0) {
      mkdir("output_data", 0755);
    }

    peps_parameters.save("output_data/output_params.dat");
    lattice.save_append("output_data/output_params.dat");
    local_parameters.save_append("output_data/output_params.dat");
    std::ofstream ofs;
    ofs.open("output_data/output_params.dat", std::ios::out | std::ios::app);
    ofs << std::endl;
    ofs.close();
  }

  // set seed for randomized svd
  random_tensor::set_seed(11 + mpirank);

  // for convenience//
  int D = peps_parameters.D;
  int CHI = peps_parameters.CHI;

  int LX = lattice.LX;
  int LY = lattice.LY;
  int N_UNIT = lattice.N_UNIT;

  // Tensors
  std::vector<ptensor> Tn(N_UNIT, ptensor(Shape(D, D, D, D, 2)));
  std::vector<ptensor>
      eTt(N_UNIT, ptensor(Shape(CHI, CHI, D, D))),
      eTr(N_UNIT, ptensor(Shape(CHI, CHI, D, D))),
      eTb(N_UNIT, ptensor(Shape(CHI, CHI, D, D))),
      eTl(N_UNIT, ptensor(Shape(CHI, CHI, D, D)));
  std::vector<ptensor> 
      C1(N_UNIT, ptensor(Shape(CHI, CHI))),
      C2(N_UNIT, ptensor(Shape(CHI, CHI))),
      C3(N_UNIT, ptensor(Shape(CHI, CHI))),
      C4(N_UNIT, ptensor(Shape(CHI, CHI)));

  std::vector<std::vector<std::vector<double> > > lambda_tensor(
      N_UNIT, std::vector<std::vector<double> >(nleg, std::vector<double>(D)));

  Initialize_Tensors(Tn, lattice);
  for (int i1 = 0; i1 < N_UNIT; ++i1) {
    for (int i2 = 0; i2 < nleg; ++i2) {
      for (int i3 = 0; i3 < D; ++i3) {
        lambda_tensor[i1][i2][i3] = 1.0;
      }
    }
  }


  ptensor Tn1_new, Tn2_new;
  std::vector<double> lambda_c;

  std::vector<ptensor> ops;
  for(auto Ham: hams){
    ptensor U = EvolutionaryTensor(Ham, local_parameters.tau_simple);
    ptensor op12 = transpose(reshape(U, Shape(2, 2, 2, 2)), Axes(2, 3, 0, 1));
    ops.push_back(op12);
  }

  // simple update
  start_time = MPI_Wtime();
  for (int int_tau = 0; int_tau < local_parameters.num_simple_step; ++int_tau) {
    for(auto ed: simple_edges){
      const int source = ed.source_site;
      const int target = ed.target_site;
      const int source_leg = ed.source_leg;
      const int target_leg = ed.target_leg;
      Simple_update_bond(Tn[source], Tn[target],
          lambda_tensor[source], lambda_tensor[target],
          ops[ed.op_id], source_leg, peps_parameters,
          Tn1_new, Tn2_new, lambda_c
          );
      lambda_tensor[source][source_leg] = lambda_c;
      lambda_tensor[target][target_leg] = lambda_c;
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;
    }
  }
  time_simple_update += MPI_Wtime() - start_time;
  // done simple update

  // Start full update
  if (local_parameters.num_full_step > 0) {
    ops.clear();
    for(auto Ham: hams){
      ptensor U = EvolutionaryTensor(Ham, local_parameters.tau_full);
      ptensor op12 = transpose(reshape(U, Shape(2, 2, 2, 2)), Axes(2, 3, 0, 1));
      ops.push_back(op12);
    }
    // Environment
    start_time = MPI_Wtime();
    Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
        peps_parameters, lattice);
    time_env += MPI_Wtime() - start_time;
  }

  start_time = MPI_Wtime();
  for (int int_tau = 0; int_tau < local_parameters.num_full_step; ++int_tau) {

    for(auto ed: full_edges){
      const int source = ed.source_site;
      const int target = ed.target_site;

      if(ed.is_horizontal()){
        /*
         *  C1 t t' C2'
         *  l  T T' r'
         *  C4 b b' C3'
         */
        Full_update_bond(C1[source], C2[target], C3[target], C4[source],
            eTt[source], eTt[target], eTr[target], // t  t' r'
            eTb[target], eTb[source], eTl[source], // b' b  l
            Tn[source], Tn[target],
            ops[ed.op_id], ed.source_leg, peps_parameters,
            Tn1_new, Tn2_new);
      }else{
        /*
         * C1' t' C2'
         *  l' T' r'
         *  l  T  r
         * C4  b  C3
         *
         *   |
         *   | rotate
         *   V
         * 
         *  C4 l l' C1'
         *  b  T T' t'
         *  C3 r r' C2'
         */
        Full_update_bond(C4[source], C1[target], C2[target], C3[source],
            eTl[source], eTl[target], eTt[target], // l  l' t'
            eTr[target], eTr[source], eTb[source], // r' r  b
            Tn[source], Tn[target],
            ops[ed.op_id], ed.source_leg, peps_parameters,
            Tn1_new, Tn2_new);
      }
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;

      if(peps_parameters.Full_Use_FFU){
        if(ed.is_horizontal()){
          const int source_x = source % LX;
          const int target_x = target % LX;
          Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
              source_x, peps_parameters, lattice);
          Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
              target_x, peps_parameters, lattice);
        }else{
          const int source_y = source / LX;
          const int target_y = target / LX;
          Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
              source_y, peps_parameters, lattice);
          Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
              target_y, peps_parameters, lattice);
        }
      }else{
        Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
            peps_parameters, lattice);
      }
    }
  }
  time_full_update += MPI_Wtime() - start_time;
  // done full update

  // Calc physical quantities

  start_time = MPI_Wtime();
  Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
      peps_parameters, lattice);
  time_env += MPI_Wtime() - start_time;

  ptensor op_identity(Shape(2, 2)), op_mz(Shape(2, 2)), op_mx(Shape(2, 2));

  op_identity.set_value(Index(0, 0), 1.0);
  op_identity.set_value(Index(0, 1), 0.0);
  op_identity.set_value(Index(1, 0), 0.0);
  op_identity.set_value(Index(1, 1), 1.0);

  op_mz.set_value(Index(0, 0), 0.5);
  op_mz.set_value(Index(0, 1), 0.0);
  op_mz.set_value(Index(1, 0), 0.0);
  op_mz.set_value(Index(1, 1), -0.5);

  op_mx.set_value(Index(0, 0), 0.0);
  op_mx.set_value(Index(0, 1), 0.5);
  op_mx.set_value(Index(1, 0), 0.5);
  op_mx.set_value(Index(1, 1), 0.0);

  std::vector<double> mz(N_UNIT), mx(N_UNIT);
  std::vector<std::vector<double> > zz(N_UNIT, std::vector<double>(2));
  for (int i = 0; i < N_UNIT; ++i) {
    mx[i] = 0.0;
    mz[i] = 0.0;
    zz[i][0] = 0.0;
    zz[i][1] = 0.0;
  }

  double norm;
  int num_j;
  start_time = MPI_Wtime();
  for (int i = 0; i < N_UNIT; ++i) {
    norm = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
        eTb[i], eTl[i], Tn[i], op_identity);
    mz[i] = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
        eTb[i], eTl[i], Tn[i], op_mz)
    / norm;
    mx[i] = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
        eTb[i], eTl[i], Tn[i], op_mx)
    / norm;
    if (peps_parameters.Debug_flag) {
      std::cout << i << " "
      << " " << norm << " " << mz[i] << " " << mx[i] << std::endl;
    }
  }
  double norm_x, norm_y;
  for (int num = 0; num < N_UNIT; ++num) {
    num_j = lattice.NN_Tensor[num][2];

    // x direction
    norm_x = Contract_two_sites_holizontal(
        C1[num], C2[num_j], C3[num_j], C4[num], eTt[num], eTt[num_j],
        eTr[num_j], eTb[num_j], eTb[num], eTl[num], Tn[num], Tn[num_j],
        op_identity, op_identity);
    zz[num][0] = Contract_two_sites_holizontal(
        C1[num], C2[num_j], C3[num_j], C4[num], eTt[num],
        eTt[num_j], eTr[num_j], eTb[num_j], eTb[num], eTl[num],
        Tn[num], Tn[num_j], op_mz, op_mz) /
    norm_x;

    // y direction
    num_j = lattice.NN_Tensor[num][3];

    norm_y = Contract_two_sites_vertical(
        C1[num], C2[num], C3[num_j], C4[num_j], eTt[num], eTr[num],
        eTr[num_j], eTb[num_j], eTl[num_j], eTl[num], Tn[num], Tn[num_j],
        op_identity, op_identity);
    zz[num][1] = Contract_two_sites_vertical(
        C1[num], C2[num], C3[num_j], C4[num_j], eTt[num],
        eTr[num], eTr[num_j], eTb[num_j], eTl[num_j], eTl[num],
        Tn[num], Tn[num_j], op_mz, op_mz) /
    norm_y;

    if (peps_parameters.Debug_flag) {
      std::cout <<  num << " " << norm_x << " " << norm_y << " "
      << zz[num][0] << " " << zz[num][1] << std::endl;
    }
  }
  time_obs += MPI_Wtime() - start_time;
  double sum_mz, sum_mx, sum_zz;
  sum_zz = 0.0;
  sum_mx = 0.0;
  sum_mz = 0.0;
  for (int i = 0; i < N_UNIT; ++i) {
    sum_mx += mx[i];
    sum_mz += mz[i];
    sum_zz += zz[i][0] + zz[i][1];
  }
  if (mpirank == 0) {
    // TODO 
    const double hx = 1.4;
    std::cout <<  -sum_zz / N_UNIT - hx * sum_mx / N_UNIT << " "
    << sum_mz / N_UNIT << " " << sum_mx / N_UNIT << " "
    << sum_zz / N_UNIT << std::endl;
  }

  if (mpirank == 0) {
    std::cout << "##time simple update= " << time_simple_update << std::endl;
    std::cout << "##time full update= " << time_full_update << std::endl;
    std::cout << "##time environmnent= " << time_env << std::endl;
    std::cout << "##time observable= " << time_obs << std::endl;
  }
  return 0;
}

using d_tensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;
using c_tensor = mptensor::Tensor<mptensor::scalapack::Matrix, std::complex<double>>;

template
int tnsolve<d_tensor>(MPI_Comm comm,
                      PEPS_Parameters peps_parameters,
                      Parameters local_parameters,
                      Lattice lattice,
                      Edges simple_edges,
                      Edges full_edges,
                      std::vector<d_tensor> hams
    );
/*
template
int tnsolve<c_tensor>(MPI_Comm comm,
                      PEPS_Parameters peps_parameters,
                      Parameters local_parameters,
                      Lattice lattice,
                      Edges simple_edges,
                      Edges full_edges,
                      std::vector<c_tensor> hams
    );
    */
