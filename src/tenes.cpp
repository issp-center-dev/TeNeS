#define _USE_MATH_DEFINES
#include <algorithm>
#include <complex>
#include <random>
#include <sys/stat.h>

#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>

#include "type.hpp"

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "Square_lattice_CTM.hpp"
#include "correlation.hpp"
#include "edge.hpp"
#include "timer.hpp"

#include "tenes.hpp"

namespace tenes {

using namespace mptensor;

template <class ptensor> class TeNeS {
public:
  TeNeS(MPI_Comm comm_, PEPS_Parameters peps_parameters_, Lattice lattice_,
        Operators<ptensor> simple_updates_, Operators<ptensor> full_updates_,
        Edges ham_edges_, std::vector<ptensor> hams_,
        std::vector<ptensor> lops_, CorrelationParameter corparam_);

  void initialize_tensors();
  void update_CTM();
  void simple_update();
  void full_update();

  void measure();
  std::vector<std::vector<double>> measure_local(bool save);
  double measure_energy(bool save);
  std::vector<std::vector<std::vector<double>>> measure_NN(bool save);
  std::vector<Correlation> measure_correlation(bool save);
  void optimize();
  void save_tensors() const;

private:
  static constexpr int nleg = 4;

  MPI_Comm comm;
  int mpisize, mpirank;

  PEPS_Parameters peps_parameters;
  Lattice lattice;

  Operators<ptensor> simple_updates;
  Operators<ptensor> full_updates;
  Edges ham_edges;
  std::vector<ptensor> hams;
  std::vector<ptensor> lops;
  ptensor op_identity;

  CorrelationParameter corparam;

  std::vector<ptensor> Tn;
  std::vector<ptensor> eTt, eTr, eTb, eTl;
  std::vector<ptensor> C1, C2, C3, C4;
  std::vector<std::vector<std::vector<double>>> lambda_tensor;

  int D;
  int CHI;
  int LX;
  int LY;
  int N_UNIT;
  int ldof;

  std::string outdir;

  double time_simple_update;
  double time_full_update;
  double time_environment;
  double time_observable;
};

template <class ptensor>
TeNeS<ptensor>::TeNeS(MPI_Comm comm_, PEPS_Parameters peps_parameters_,
                      Lattice lattice_, Operators<ptensor> simple_updates_,
                      Operators<ptensor> full_updates_, Edges ham_edges_,
                      std::vector<ptensor> hams_, std::vector<ptensor> lops_,
                      CorrelationParameter corparam_)
    : comm(comm_), peps_parameters(peps_parameters_), lattice(lattice_),
      simple_updates(simple_updates_), full_updates(full_updates_),
      ham_edges(ham_edges_), lops(lops_), corparam(corparam_), outdir("output"),
      time_simple_update(), time_full_update(), time_environment(),
      time_observable() {
  MPI_Comm_size(comm, &mpisize);
  MPI_Comm_rank(comm, &mpirank);

  peps_parameters.Bcast(comm);
  // output debug or warning info only from process 0
  if (mpirank != 0) {
    peps_parameters.print_level = PEPS_Parameters::PrintLevel::none;
  }

  D = peps_parameters.D;
  CHI = peps_parameters.CHI;

  lattice.Bcast(comm);

  LX = lattice.LX;
  LY = lattice.LY;
  N_UNIT = lattice.N_UNIT;
  ldof = lops.begin()->shape()[0];

  // set seed for randomized svd
  int seed = peps_parameters.seed;
  random_tensor::set_seed(seed + mpirank);

  if (outdir.empty()) {
    outdir += ".";
  }

  if (mpirank == 0) {
    // folder check
    struct stat status;
    if (stat(outdir.c_str(), &status) != 0) {
      mkdir(outdir.c_str(), 0755);
    }

    std::string param_file = outdir + "/parameters.dat";

    peps_parameters.save(param_file.c_str());
    lattice.save_append(param_file.c_str());
  }

  op_identity = ptensor(Shape(ldof, ldof));
  for (int i = 0; i < ldof; ++i) {
    for (int j = 0; j < i; ++j)
      op_identity.set_value(Index(i, j), 0.0);
    op_identity.set_value(Index(i, i), 1.0);
    for (int j = i + 1; j < ldof; ++j)
      op_identity.set_value(Index(i, j), 0.0);
  }

  const int LDOF = ldof;

  std::transform(hams_.begin(), hams_.end(), std::back_inserter(hams),
                 [LDOF](ptensor const &A) {
                   return transpose(reshape(A, Shape(LDOF, LDOF, LDOF, LDOF)),
                                    Axes(2, 3, 0, 1));
                 });

  initialize_tensors();
}

template <class ptensor> void TeNeS<ptensor>::initialize_tensors() {
  Tn.clear();
  eTt.clear();
  eTr.clear();
  eTb.clear();
  eTl.clear();
  C1.clear();
  C2.clear();
  C3.clear();
  C4.clear();
  lambda_tensor.clear();

  for (int i = 0; i < N_UNIT; ++i) {
    Tn.push_back(ptensor(Shape(D, D, D, D, ldof)));
    eTt.push_back(ptensor(Shape(CHI, CHI, D, D)));
    eTr.push_back(ptensor(Shape(CHI, CHI, D, D)));
    eTb.push_back(ptensor(Shape(CHI, CHI, D, D)));
    eTl.push_back(ptensor(Shape(CHI, CHI, D, D)));
    C1.push_back(ptensor(Shape(CHI, CHI)));
    C2.push_back(ptensor(Shape(CHI, CHI)));
    C3.push_back(ptensor(Shape(CHI, CHI)));
    C4.push_back(ptensor(Shape(CHI, CHI)));
    lambda_tensor.push_back(
        std::vector<std::vector<double>>(nleg, std::vector<double>(D)));
    for (int j = 0; j < nleg; ++j) {
      for (int k = 0; k < D; ++k) {
        lambda_tensor[i][j][k] = 1.0;
      }
    }
  }

  std::vector<double> ran(D * D * D * D * ldof);
  std::mt19937 gen(peps_parameters.seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  int nr;

  std::string const &load_dir = peps_parameters.tensor_load_dir;
  if (load_dir.empty()) {
    Index index;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      for (int i = 0; i < D * D * D * D * ldof; i++) {
        ran[i] = dist(gen);
      }
      auto &dir = lattice.initial_dirs[i];
      if (std::all_of(dir.begin(), dir.end(),
                      [=](double x) { return x == 0.0; })) {
        // random
        dir.resize(ldof);
        for (int j = 0; j < ldof; ++j) {
          dir[j] = dist(gen);
        }
      }
      for (int n = 0; n < Tn[i].local_size(); ++n) {
        index = Tn[i].global_index(n);
        if (index[0] == 0 && index[1] == 0 && index[2] == 0 && index[3] == 0) {
          Tn[i].set_value(index, dir[index[4]]);
        } else {
          nr = index[0] + index[1] * D + index[2] * D * D +
               index[3] * D * D * D + index[4] * D * D * D * D;
          Tn[i].set_value(index, lattice.noises[i] * ran[nr]);
        }
      }
    }
  } else {
    struct stat status;
    if (stat(load_dir.c_str(), &status) != 0) {
      std::string msg = load_dir + " does not exists.";
      throw std::runtime_error(msg);
    }
    for (int i = 0; i < N_UNIT; ++i) {
      std::string filename = load_dir + "/";
      std::string suffix = "_" + std::to_string(i) + ".dat";
      Tn[i].load((filename + "T" + suffix).c_str());
      eTt[i].load((filename + "Et" + suffix).c_str());
      eTr[i].load((filename + "Er" + suffix).c_str());
      eTb[i].load((filename + "Eb" + suffix).c_str());
      eTl[i].load((filename + "El" + suffix).c_str());
      C1[i].load((filename + "C1" + suffix).c_str());
      C2[i].load((filename + "C2" + suffix).c_str());
      C3[i].load((filename + "C3" + suffix).c_str());
      C4[i].load((filename + "C4" + suffix).c_str());
    }
    int nl = N_UNIT * nleg * D;
    std::vector<double> ls;
    if (mpirank == 0) {
      for (int i = 0; i < N_UNIT; ++i) {
        std::ifstream ifs(load_dir + "/lambda_" + std::to_string(i) + ".dat");
        for (int j = 0; j < nleg; ++j) {
          for (int k = 0; k < D; ++k) {
            double temp = 0.0;
            ifs >> temp;
            ls.push_back(temp);
          }
        }
      }
    }
    MPI_Bcast(&(ls[0]), nl, MPI_DOUBLE, 0, comm);
    int index = 0;
    for (int i = 0; i < N_UNIT; ++i) {
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < D; ++k) {
          lambda_tensor[i][j][k] = ls[index];
          ++index;
        }
      }
    }
  } // end of else part of if(load_dir.empty())
}

template <class ptensor> inline void TeNeS<ptensor>::update_CTM() {
  Timer<> timer;
  Calc_CTM_Environment(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, peps_parameters,
                       lattice);
  time_environment += timer.elapsed();
}

template <class ptensor> void TeNeS<ptensor>::simple_update() {
  Timer<> timer;
  ptensor Tn1_new;
  ptensor Tn2_new;
  std::vector<double> lambda_c;
  const int nsteps = peps_parameters.num_simple_step;

  int ireport = 1;
  int next_report_step = 0.1 * nsteps - 1;

  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : simple_updates) {
      const int source = up.source_site;
      const int source_leg = up.source_leg;
      const int target = lattice.neighbor(source, source_leg);
      const int target_leg = (source_leg + 2) % 4;
      Simple_update_bond(Tn[source], Tn[target], lambda_tensor[source],
                         lambda_tensor[target], up.op, source_leg,
                         peps_parameters, Tn1_new, Tn2_new, lambda_c);
      lambda_tensor[source][source_leg] = lambda_c;
      lambda_tensor[target][target_leg] = lambda_c;
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;
    }

    if (mpirank == 0 &&
        peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      if (int_tau == next_report_step) {
        std::cout << 100.0 * (int_tau + 1) / nsteps << "% done" << std::endl;
        ++ireport;
        next_report_step = 0.1 * ireport * nsteps - 1;
      }
    }
  }
  time_simple_update += timer.elapsed();
}

template <class ptensor> void TeNeS<ptensor>::full_update() {
  Timer<> timer;

  ptensor Tn1_new, Tn2_new;
  if (peps_parameters.num_full_step > 0) {
    update_CTM();
  }

  const int nsteps = peps_parameters.num_full_step;

  int ireport = 1;
  int next_report_step = 0.1 * nsteps - 1;

  timer.reset();
  for (int int_tau = 0; int_tau < nsteps; ++int_tau) {
    for (auto up : full_updates) {
      const int source = up.source_site;
      const int source_leg = up.source_leg;
      const int target = lattice.neighbor(source, source_leg);
      // const int target_leg = (source_leg + 2) % 4;

      switch (source_leg) {
      case 0:
        /*
         *  C1' t' t C3
         *  l'  T' T r
         *  C2' b' b C4
         *
         *   |
         *   | rotate
         *   V
         *
         *  C4 b b' C2'
         *  r  T T' l'
         *  C3 t t' C1'
         */
        Full_update_bond(C4[source], C2[target], C1[target], C3[source],
                         eTb[source], eTb[target], eTl[target], eTt[target],
                         eTt[source], eTr[source], Tn[source], Tn[target],
                         up.op, source_leg, peps_parameters, Tn1_new, Tn2_new);
        break;
      case 1:
        /*
         * C1' t' C2'
         *  l' T'  r'
         *  l  T   r
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
                         eTl[source], eTl[target], eTt[target], eTr[target],
                         eTr[source], eTb[source], Tn[source], Tn[target],
                         up.op, source_leg, peps_parameters, Tn1_new, Tn2_new);
        break;
      case 2:
        /*
         *  C1 t t' C2'
         *  l  T T' r'
         *  C4 b b' C3'
         */
        Full_update_bond(C1[source], C2[target], C3[target], C4[source],
                         eTt[source], eTt[target], eTr[target], // t  t' r'
                         eTb[target], eTb[source], eTl[source], // b' b  l
                         Tn[source], Tn[target], up.op, source_leg,
                         peps_parameters, Tn1_new, Tn2_new);
        break;
      case 3:
        /*
         * C1  t C2
         *  l  T  r
         *  l' T' r'
         * C4' b C3'
         *
         *   |
         *   | rotate
         *   V
         *
         *  C2 r r' C3'
         *  t  T T' b'
         *  C1 l l' C4'
         */
        Full_update_bond(C2[source], C3[target], C4[target], C1[source],
                         eTr[source], eTr[target], eTb[target], eTl[target],
                         eTl[source], eTt[source], Tn[source], Tn[target],
                         up.op, source_leg, peps_parameters, Tn1_new, Tn2_new);
        break;

      default:
        break;
      } // end of switch
      Tn[source] = Tn1_new;
      Tn[target] = Tn2_new;

      if (peps_parameters.Full_Use_FastFullUpdate) {
        if (up.is_horizontal()) {
          const int source_x = source % LX;
          const int target_x = target % LX;
          Left_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_x,
                    peps_parameters, lattice);
          Right_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_x,
                     peps_parameters, lattice);
        } else {
          const int source_y = source / LX;
          const int target_y = target / LX;
          Top_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, source_y,
                   peps_parameters, lattice);
          Bottom_move(C1, C2, C3, C4, eTt, eTr, eTb, eTl, Tn, target_y,
                      peps_parameters, lattice);
        }
      } else {
        update_CTM();
      }
    }

    if (mpirank == 0 &&
        peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      if (int_tau == next_report_step) {
        std::cout << 100.0 * (int_tau + 1) / nsteps << "% done" << std::endl;
        ++ireport;
        next_report_step = 0.1 * ireport * nsteps - 1;
      }
    }
  }
  time_full_update += timer.elapsed();
}

template <class ptensor> void TeNeS<ptensor>::optimize() {
  // for measure time
  std::cout << std::setprecision(12);

  ptensor Tn1_new, Tn2_new;
  std::vector<double> lambda_c;

  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "Start simple update" << std::endl;
  }
  simple_update();

  if (peps_parameters.num_full_step > 0) {
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "Start full update" << std::endl;
    }
    full_update();
  }
}

template <class ptensor>
std::vector<std::vector<double>> TeNeS<ptensor>::measure_local(bool save) {
  Timer<> timer;
  const int nlops = lops.size();
  std::vector<std::vector<double>> local_obs(nlops,
                                             std::vector<double>(N_UNIT, 0));

  for (int i = 0; i < N_UNIT; ++i) {
    double norm = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                    eTb[i], eTl[i], Tn[i], op_identity);
    for (int ilops = 0; ilops < nlops; ++ilops) {
      double val = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                     eTb[i], eTl[i], Tn[i], lops[ilops]) /
                   norm;
      local_obs[ilops][i] = val;
    }
  }

  time_observable += timer.elapsed();

  if (save && mpirank == 0) {
    std::string filename = outdir + "/site_obs.dat";
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save site observables to " << filename << std::endl;
    }
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << "# $1: op_index\n";
    ofs << "# $2: site_index\n";
    ofs << "# $3: real\n";
    ofs << "# $4: imag\n";
    ofs << std::endl;

    for (int ilops = 0; ilops < nlops; ++ilops) {
      double sum = 0.0;
      for (int i = 0; i < N_UNIT; ++i) {
        sum += local_obs[ilops][i];
        ofs << ilops << " " << i << " " << local_obs[ilops][i] << " " << 0.0
            << std::endl;
      }
    }
  }
  return local_obs;
}

template <class ptensor> double TeNeS<ptensor>::measure_energy(bool save) {
  Timer<> timer;
  double energy = 0.0;

  for (auto ed : ham_edges) {
    const int source = ed.source_site;
    const int target = ed.target_site;
    if (ed.is_horizontal()) {
      const double local_norm = Contract_two_sites_horizontal(
          C1[source], C2[target], C3[target], C4[source], eTt[source],
          eTt[target], eTr[target], eTb[target], eTb[source], eTl[source],
          Tn[source], Tn[target], op_identity, op_identity);
      energy += Contract_two_sites_horizontal_op12(
                    C1[source], C2[target], C3[target], C4[source], eTt[source],
                    eTt[target], eTr[target], eTb[target], eTb[source],
                    eTl[source], Tn[source], Tn[target], hams[ed.op_id]) /
                local_norm;
    } else {
      const double local_norm = Contract_two_sites_vertical(
          C1[source], C2[source], C3[target], C4[target], eTt[source],
          eTr[source], eTr[target], eTb[target], eTl[target], eTl[source],
          Tn[source], Tn[target], op_identity, op_identity);
      energy += Contract_two_sites_vertical_op12(
                    C1[source], C2[source], C3[target], C4[target], eTt[source],
                    eTr[source], eTr[target], eTb[target], eTl[target],
                    eTl[source], Tn[source], Tn[target], hams[ed.op_id]) /
                local_norm;
    }
  }
  energy /= N_UNIT;

  time_observable += timer.elapsed();

  if (save && mpirank == 0) {
    std::string filename = outdir + "/energy.dat";
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save energy to " << filename << std::endl;
    }
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << energy << std::endl;
  }

  return energy;
}

template <class ptensor>
std::vector<std::vector<std::vector<double>>>
TeNeS<ptensor>::measure_NN(bool save) {
  Timer<> timer;
  const int nlops = lops.size();
  std::vector<std::vector<std::vector<double>>> neighbor_obs(
      nlops,
      std::vector<std::vector<double>>(N_UNIT, std::vector<double>(2, 0.0)));

  for (int source = 0; source < N_UNIT; ++source) {
    { // horizontal
      int target = lattice.right(source);
      const double local_norm = Contract_two_sites_horizontal(
          C1[source], C2[target], C3[target], C4[source], eTt[source],
          eTt[target], eTr[target], eTb[target], eTb[source], eTl[source],
          Tn[source], Tn[target], op_identity, op_identity);
      for (int ilops = 0; ilops < nlops; ++ilops) {
        double val =
            Contract_two_sites_horizontal(
                C1[source], C2[target], C3[target], C4[source], eTt[source],
                eTt[target], eTr[target], eTb[target], eTb[source], eTl[source],
                Tn[source], Tn[target], lops[ilops], lops[ilops]) /
            local_norm;
        neighbor_obs[ilops][source][0] = val;
      }
    }
    { // vertical
      int target = lattice.top(source);
      const double local_norm = Contract_two_sites_vertical(
          C1[target], C2[target], C3[source], C4[source], eTt[target],
          eTr[target], eTr[source], eTb[source], eTl[source], eTl[target],
          Tn[target], Tn[source], op_identity, op_identity);
      for (int ilops = 0; ilops < nlops; ++ilops) {
        double val =
            Contract_two_sites_vertical(
                C1[target], C2[target], C3[source], C4[source], eTt[target],
                eTr[target], eTr[source], eTb[source], eTl[source], eTl[target],
                Tn[target], Tn[source], lops[ilops], lops[ilops]) /
            local_norm;
        neighbor_obs[ilops][source][1] = val;
      }
    }
  }
  time_observable += timer.elapsed();

  if (save && mpirank == 0) {
    std::string filename = outdir + "/neighbor_obs.dat";
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save NN correlation to " << filename << std::endl;
    }
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << "# $1: op_index\n";
    ofs << "# $2: source_site\n";
    ofs << "# $3: target_site\n";
    ofs << "# $4: real\n";
    ofs << "# $5: imag\n";
    ofs << std::endl;
    for (int ilops = 0; ilops < nlops; ++ilops) {
      for (int source = 0; source < N_UNIT; ++source) {
        int target = lattice.right(source);
        ofs << ilops << " " << source << " " << target << " "
            << neighbor_obs[ilops][source][0] << " " << 0.0 << std::endl;
        target = lattice.top(source);
        ofs << ilops << " " << source << " " << target << " "
            << neighbor_obs[ilops][source][1] << " " << 0.0 << std::endl;
      }
    }
  }

  return neighbor_obs;
}

template <class ptensor>
std::vector<Correlation> TeNeS<ptensor>::measure_correlation(bool save) {
  Timer<> timer;

  const int nlops = lops.size();
  const int r_max = corparam.r_max;
  std::vector<std::vector<int>> r_ops(nlops);
  for (auto ops : corparam.operators) {
    // TODO: range check
    r_ops[std::get<0>(ops)].push_back(std::get<1>(ops));
  }

  ptensor correlation_T(Shape(CHI, CHI, D, D));
  ptensor correlation_norm(Shape(CHI, CHI, D, D));
  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    for (int left_ilop = 0; left_ilop < nlops; ++left_ilop) {
      if (r_ops[left_ilop].empty()) {
        continue;
      }

      const int left_x = lattice.x(left_index);
      const int left_y = lattice.y(left_index);
      { // horizontal
        StartCorrelation(correlation_T, C1[left_index], C4[left_index],
                         eTt[left_index], eTb[left_index], eTl[left_index],
                         Tn[left_index], lops[left_ilop]);
        StartCorrelation(correlation_norm, C1[left_index], C4[left_index],
                         eTt[left_index], eTb[left_index], eTl[left_index],
                         Tn[left_index], op_identity);

        for (int r = 0; r < r_max; ++r) {
          const int right_x = (left_x + r + 1) % LX;
          const int right_y = left_y;
          const int offset_x = (left_x + r + 1) / LX;
          const int offset_y = 0;
          const int right_index = lattice.index(right_x, right_y);
          double norm = FinishCorrelation(correlation_norm, C2[right_index],
                                          C3[right_index], eTt[right_index],
                                          eTr[right_index], eTb[right_index],
                                          Tn[right_index], op_identity);
          for (auto right_ilop : r_ops[left_ilop]) {
            double val = FinishCorrelation(correlation_T, C2[right_index],
                                           C3[right_index], eTt[right_index],
                                           eTr[right_index], eTb[right_index],
                                           Tn[right_index], lops[left_ilop]) /
                         norm;
            correlations.push_back(Correlation{left_index, right_index,
                                               offset_x, offset_y, left_ilop,
                                               right_ilop, val, 0.0});
          }

          Transfer(correlation_T, eTt[right_index], eTb[right_index],
                   Tn[right_index]);
          Transfer(correlation_norm, eTt[right_index], eTb[right_index],
                   Tn[right_index]);
        }
      }
      { // vertical
        ptensor tn = transpose(Tn[left_index], Axes(3, 0, 1, 2, 4));
        StartCorrelation(correlation_T, C4[left_index], C3[left_index],
                         eTl[left_index], eTr[left_index], eTb[left_index], tn,
                         lops[left_ilop]);
        StartCorrelation(correlation_norm, C4[left_index], C3[left_index],
                         eTl[left_index], eTr[left_index], eTb[left_index], tn,
                         op_identity);

        for (int r = 0; r < r_max; ++r) {
          const int right_x = left_x;
          const int right_y = (left_y + r + 1) % LY;
          const int offset_x = 0;
          const int offset_y = (left_y + r + 1) / LY;
          const int right_index = lattice.index(right_x, right_y);
          tn = transpose(Tn[right_index], Axes(3, 0, 1, 2, 4));
          double norm = FinishCorrelation(correlation_norm, C1[right_index],
                                          C2[right_index], eTl[right_index],
                                          eTt[right_index], eTr[right_index],
                                          tn, op_identity);
          for (auto right_ilop : r_ops[left_ilop]) {
            double val = FinishCorrelation(correlation_T, C1[right_index],
                                           C2[right_index], eTl[right_index],
                                           eTt[right_index], eTr[right_index],
                                           tn, lops[left_ilop]) /
                         norm;
            correlations.push_back(Correlation{left_index, right_index,
                                               offset_x, offset_y, left_ilop,
                                               right_ilop, val, 0.0});
          }

          Transfer(correlation_T, eTl[right_index], eTr[right_index], tn);
          Transfer(correlation_norm, eTl[right_index], eTr[right_index], tn);
        }
      }
    }
  }

  time_observable += timer.elapsed();

  if (save && mpirank == 0) {
    std::string filename = outdir + "/correlation.dat";
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save long-range correlations to " << filename
                << std::endl;
    }
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << "# $1: left_op\n";
    ofs << "# $2: left_site\n";
    ofs << "# $3: right_op\n";
    ofs << "# $4: right_site\n";
    ofs << "# $5: offset_x\n";
    ofs << "# $6: offset_y\n";
    ofs << "# $7: real\n";
    ofs << "# $8: imag\n";
    ofs << std::endl;
    for (auto const &cor : correlations) {
      ofs << cor.left_op << " " << cor.left_index << " " << cor.right_op << " "
          << cor.right_index << " " << cor.offset_x << " " << cor.offset_y
          << " " << cor.real << " " << cor.imag << " " << std::endl;
    }
  }
  return correlations;
}

template <class ptensor> void TeNeS<ptensor>::measure() {
  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "Start calculating observables" << std::endl;
    std::clog << "  Start updating environment" << std::endl;
  }
  update_CTM();

  const int nlops = lops.size();

  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "  Start calculating local operators" << std::endl;
  }
  auto local_obs = measure_local(true);

  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "  Start calculating energy" << std::endl;
  }
  auto energy = measure_energy(true);

  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "  Start calculating NN correlation" << std::endl;
  }
  auto NN_obs = measure_NN(true);

  if (corparam.r_max > 0) {
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "  Start calculating long range correlation" << std::endl;
    }
    auto correlations = measure_correlation(true);
  }

  if (mpirank == 0) {
    std::string filename = outdir + "/time.dat";
    std::ofstream ofs(filename.c_str());
    ofs << "time simple update = " << time_simple_update << std::endl;
    ofs << "time full update   = " << time_full_update << std::endl;
    ofs << "time environmnent  = " << time_environment << std::endl;
    ofs << "time observable    = " << time_observable << std::endl;
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save elapsed times to " << filename << std::endl;
    }

    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::cout << std::endl;

      std::cout << "Energy = " << energy << std::endl;

      for (int ilops = 0; ilops < nlops; ++ilops) {
        double sum = 0.0;
        for (int i = 0; i < N_UNIT; ++i) {
          sum += local_obs[ilops][i];
        }
        std::cout << "Local operator " << ilops << " = " << sum / N_UNIT
                  << std::endl;
      }
      std::cout << std::endl;

      std::cout << "time simple update = " << time_simple_update << std::endl;
      std::cout << "time full update   = " << time_full_update << std::endl;
      std::cout << "time environmnent  = " << time_environment << std::endl;
      std::cout << "time observable    = " << time_observable << std::endl;
    }
  }
}

template <class ptensor> void TeNeS<ptensor>::save_tensors() const {
  std::string const &save_dir = peps_parameters.tensor_save_dir;
  if (save_dir.empty()) {
    return;
  }
  if (mpirank == 0) {
    struct stat status;
    if (stat(save_dir.c_str(), &status) != 0) {
      mkdir(save_dir.c_str(), 0755);
    }
  }
  for (int i = 0; i < N_UNIT; ++i) {
    std::string filename = save_dir + "/";
    std::string suffix = "_" + std::to_string(i) + ".dat";
    Tn[i].save((filename + "T" + suffix).c_str());
    eTt[i].save((filename + "Et" + suffix).c_str());
    eTr[i].save((filename + "Er" + suffix).c_str());
    eTb[i].save((filename + "Eb" + suffix).c_str());
    eTl[i].save((filename + "El" + suffix).c_str());
    C1[i].save((filename + "C1" + suffix).c_str());
    C2[i].save((filename + "C2" + suffix).c_str());
    C3[i].save((filename + "C3" + suffix).c_str());
    C4[i].save((filename + "C4" + suffix).c_str());
  }
  if (mpirank == 0) {
    for (int i = 0; i < N_UNIT; ++i) {
      std::ofstream ofs(save_dir + "/lambda_" + std::to_string(i) + ".dat");
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < D; ++k) {
          ofs << lambda_tensor[i][j][k] << std::endl;
        }
      }
    }
  }
}

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
          Operators<tensor> simple_updates, Operators<tensor> full_updates,
          Edges ham_edges, std::vector<tensor> hamiltonians,
          std::vector<tensor> local_operators, CorrelationParameter corparam){
  TeNeS<tensor> tns(comm, peps_parameters, lattice, simple_updates, full_updates,
                     ham_edges, hamiltonians, local_operators, corparam);
  tns.optimize();
  tns.save_tensors();
  tns.measure();
  return 0;
}

// template specialization
using d_tensor = mptensor::Tensor<mptensor_matrix_type, double>;
template int tenes<d_tensor>(MPI_Comm comm, PEPS_Parameters peps_parameters,
                             Lattice lattice,
                             Operators<d_tensor> simple_updates,
                             Operators<d_tensor> full_updates, Edges ham_edges,
                             std::vector<d_tensor> hams,
                             std::vector<d_tensor> lops,
                             CorrelationParameter corparam);
/*
using c_tensor = mptensor::Tensor<mptensor_matrix_type,
std::complex<double>>;
template int tenes<c_tensor>(MPI_Comm comm,
                      PEPS_Parameters peps_parameters,
                      Lattice lattice,
                      Edges simple_edges,
                      Edges full_edges,
                      std::vector<c_tensor> hams,
                      std::vector<d_tensor> lops,
                      CorrelationParameter corparam
    );
    */

} // end of namespace tenes
