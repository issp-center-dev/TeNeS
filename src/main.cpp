#define _USE_MATH_DEFINES
#include <random>
#include <sys/stat.h>

#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>

#include <toml11/toml.hpp>

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "Parameters.hpp"
#include "Square_lattice_CTM.hpp"


using namespace mptensor;
typedef Tensor<scalapack::Matrix, double> ptensor;
typedef Tensor<lapack::Matrix, double> tensor;

/*
 * axis:
 *
 *  y
 *  ^
 *  |
 *  .->x
 *
 * edge index:
 *
 *   1
 *  0.2
 *   3
 *
 * corner and edge tensor
 *
 * C1 t C2
 * l  .  r
 * C4 b C3
 *
 */

struct Edge{
  enum direction{ horizontal, vertical };
  int source_site;
  int source_leg;
  int target_site;
  int target_leg;
  direction dir;
  Edge():source_site(-1), source_leg(-1), target_site(-1), target_leg(-1), dir(horizontal){}

  /*
   * horizontal:
   *
   *  s-t
   *
   * vertical:
   *
   *  t
   *  |
   *  s
   */
  Edge(int source_site, int target_site, direction dir)
      : source_site(source_site), source_leg(dir==horizontal?2:1),
        target_site(target_site), target_leg(dir==horizontal?0:3),
        dir(dir){}

  bool is_horizontal() const { return dir==horizontal; }
  bool is_vertical() const { return !is_horizontal(); }
};


/* for MPI */
int mpirank;
int mpisize;

ptensor Set_Hamiltonian(double hx) {
  ptensor Ham(Shape(4, 4));
  Ham.set_value(Index(0, 0), -0.25);
  Ham.set_value(Index(0, 1), -0.125 * hx);
  Ham.set_value(Index(0, 2), -0.125 * hx);
  Ham.set_value(Index(0, 3), 0.0);

  Ham.set_value(Index(1, 0), -0.125 * hx);
  Ham.set_value(Index(1, 1), 0.25);
  Ham.set_value(Index(1, 2), 0.0);
  Ham.set_value(Index(1, 3), -0.125 * hx);

  Ham.set_value(Index(2, 0), -0.125 * hx);
  Ham.set_value(Index(2, 1), 0.0);
  Ham.set_value(Index(2, 2), 0.25);
  Ham.set_value(Index(2, 3), -0.125 * hx);

  Ham.set_value(Index(3, 0), 0.0);
  Ham.set_value(Index(3, 1), -0.125 * hx);
  Ham.set_value(Index(3, 2), -0.125 * hx);
  Ham.set_value(Index(3, 3), -0.25);
  return Ham;
}

void Initialize_Tensors(std::vector<ptensor> &Tn, Lattice lattice) {
  // ferro
  int D = Tn[0].shape()[0];

  std::vector<double> ran(D * D * D * D * 2);
#ifdef CPP11
  std::mt19937 gen(11);
  std::uniform_real_distribution<double> dist(-0.01, 0.01);
  for (int i = 0; i < D * D * D * D * 2; i++) {
    ran[i] = dist(gen);
  }
#elif BOOST
  boost::mt19937 gen(11);
  boost::uniform_real<double> dist(-0.01, 0.01);
  for (int i = 0; i < D * D * D * D * 2; i++) {
    ran[i] = dist(gen);
  }
#else
  for (int i = 0; i < D * D * D * D * 2; i++) {
    ran[i] = 0.0;
  }
#endif
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

int main(int argc, char **argv) {

  /* MPI initialization */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  // for measure time
  double time_simple_update = 0.0;
  double time_full_update = 0.0;
  double time_env = 0.0;
  double time_obs = 0.0;
  double start_time;

  std::cout << std::setprecision(12);

  // Parameters
  PEPS_Parameters peps_parameters;
  Parameters local_parameters;

  constexpr int nleg = 4;
  Lattice lattice(2,2);

  if (mpirank == 0) {
    local_parameters.set("input.toml");
    peps_parameters.set("input.toml");
  }

  local_parameters.Bcast_parameters(MPI_COMM_WORLD);
  peps_parameters.Bcast_parameters(MPI_COMM_WORLD);

  lattice.Bcast_parameters(MPI_COMM_WORLD);

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

    peps_parameters.output_parameters("output_data/output_params.dat");
    lattice.output_parameters_append("output_data/output_params.dat");
    local_parameters.output_parameters_append("output_data/output_params.dat");
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
  const double hx = 1.4;

  ptensor Ham = Set_Hamiltonian(hx);
  ptensor U = EvolutionaryTensor(Ham, local_parameters.tau_simple);
  ptensor op12 = transpose(reshape(U, Shape(2, 2, 2, 2)), Axes(2, 3, 0, 1));

  std::vector<Edge> simple_edges;
  // x-bond A sub-lattice
  simple_edges.push_back(Edge(0, 1, Edge::horizontal));
  simple_edges.push_back(Edge(3, 2, Edge::horizontal));
  // x-bond B sub-lattice
  simple_edges.push_back(Edge(1, 0, Edge::horizontal));
  simple_edges.push_back(Edge(2, 3, Edge::horizontal));
  // y-bond A sub-lattice
  simple_edges.push_back(Edge(0, 2, Edge::vertical));
  simple_edges.push_back(Edge(3, 1, Edge::vertical));
  // y-bond B sub-lattice
  simple_edges.push_back(Edge(1, 3, Edge::vertical));
  simple_edges.push_back(Edge(2, 0, Edge::vertical));

  ptensor Tn1_new, Tn2_new;
  std::vector<double> lambda_c;

  // simple update
  start_time = MPI_Wtime();
  for (int int_tau = 0; int_tau < local_parameters.tau_simple_step; ++int_tau) {
    for(auto ed: simple_edges){
      const int source = ed.source_site;
      const int target = ed.target_site;
      const int source_leg = ed.source_leg;
      const int target_leg = ed.target_leg;
      Simple_update_bond(Tn[source], Tn[target],
          lambda_tensor[source], lambda_tensor[target],
          op12, source_leg, peps_parameters,
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
  if (local_parameters.tau_full_step > 0) {
    Ham = Set_Hamiltonian(hx);
    EvolutionaryTensor(U, Ham, local_parameters.tau_full_step);
    op12 = transpose(reshape(U, Shape(2, 2, 2, 2)), Axes(2, 3, 0, 1));

    // Environment
    start_time = MPI_Wtime();
    Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn,
        peps_parameters, lattice);
    time_env += MPI_Wtime() - start_time;
  }

  std::vector<Edge> full_edges;
  // x-bond A sub-lattice
  full_edges.push_back(Edge(0, 1, Edge::horizontal));
  full_edges.push_back(Edge(3, 2, Edge::horizontal));
  // x-bond B sub-lattice
  full_edges.push_back(Edge(1, 0, Edge::horizontal));
  full_edges.push_back(Edge(2, 3, Edge::horizontal));
  // y-bond A sub-lattice
  full_edges.push_back(Edge(0, 2, Edge::vertical));
  full_edges.push_back(Edge(3, 1, Edge::vertical));
  // y-bond B sub-lattice
  full_edges.push_back(Edge(1, 3, Edge::vertical));
  full_edges.push_back(Edge(2, 0, Edge::vertical));

  start_time = MPI_Wtime();
  for (int int_tau = 0; int_tau < local_parameters.tau_full_step; ++int_tau) {

    for(auto ed: full_edges){
      const int source = ed.source_site;
      const int target = ed.target_site;

      if(ed.is_horizontal()){
        Full_update_bond(C1[source], C2[target], C3[target], C4[source],
            eTt[source], eTt[target], eTr[target],
            eTb[target], eTb[source], eTl[source],
            Tn[source], Tn[target],
            op12, ed.source_leg, peps_parameters,
            Tn1_new, Tn2_new);
      }else{
        Full_update_bond(C4[source], C1[target], C2[target], C3[source],
            eTl[source], eTl[target], eTt[target],
            eTr[target], eTr[source], eTb[source],
            Tn[source], Tn[target],
            op12, ed.source_leg, peps_parameters,
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
      std::cout << hx << " " << i << " "
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
      std::cout << hx << " " << num << " " << norm_x << " " << norm_y << " "
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
    std::cout << hx << " " << -sum_zz / N_UNIT - hx * sum_mx / N_UNIT << " "
    << sum_mz / N_UNIT << " " << sum_mx / N_UNIT << " "
    << sum_zz / N_UNIT << std::endl;
  }

  if (mpirank == 0) {
    std::cout << "##time simple update= " << time_simple_update << std::endl;
    std::cout << "##time full update= " << time_full_update << std::endl;
    std::cout << "##time environmnent= " << time_env << std::endl;
    std::cout << "##time observable= " << time_obs << std::endl;
  }
  MPI_Finalize();
  return 0;
}
