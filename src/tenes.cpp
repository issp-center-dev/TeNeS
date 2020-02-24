#define _USE_MATH_DEFINES
#include <algorithm>
#include <complex>
#include <limits>
#include <map>
#include <random>
#include <sys/stat.h>
#include <tuple>

#include <mptensor/complex.hpp>
#include <mptensor/rsvd.hpp>
#include <mptensor/tensor.hpp>

#include "type.hpp"

#include "Lattice.hpp"
#include "PEPS_Basics.hpp"
#include "PEPS_Parameters.hpp"
#include "Square_lattice_CTM.hpp"
#include "correlation.hpp"
#include "timer.hpp"

#include "tenes.hpp"

namespace tenes {

using namespace mptensor;

struct Bond {
  int source_site;
  int target_site;
  int offset_x;
  int offset_y;
};

bool operator<(const Bond &a, const Bond &b) {
  return std::tie(a.source_site, a.target_site, a.offset_x, a.offset_y) <
         std::tie(b.source_site, b.target_site, b.offset_x, b.offset_y);
}

template <class ptensor> class TeNeS {
public:
  TeNeS(MPI_Comm comm_, PEPS_Parameters peps_parameters_, Lattice lattice_,
        NNOperators<ptensor> simple_updates_,
        NNOperators<ptensor> full_updates_,
        Operators<ptensor> onesite_operators,
        Operators<ptensor> twosite_operators, CorrelationParameter corparam_);

  void initialize_tensors();
  void update_CTM();
  void simple_update();
  void full_update();

  void measure();
  std::vector<std::vector<double>> measure_onesite(bool save);
  std::vector<std::map<Bond, double>> measure_twosite(bool save);
  // std::vector<Correlation> measure_correlation(bool save);
  void optimize();
  void save_tensors() const;

private:
  int siteoperator_index(int site, int group) const {
    return site_ops_indices[site][group];
  }

  static constexpr int nleg = 4;

  MPI_Comm comm;
  int mpisize, mpirank;

  PEPS_Parameters peps_parameters;
  Lattice lattice;

  NNOperators<ptensor> simple_updates;
  NNOperators<ptensor> full_updates;
  Operators<ptensor> onesite_operators;
  Operators<ptensor> twosite_operators;
  std::vector<std::vector<int>> site_ops_indices;
  int num_onesite_operators;
  int num_twosite_operators;

  std::vector<ptensor> op_identity;

  CorrelationParameter corparam;

  std::vector<ptensor> Tn;
  std::vector<ptensor> eTt, eTr, eTb, eTl;
  std::vector<ptensor> C1, C2, C3, C4;
  std::vector<std::vector<std::vector<double>>> lambda_tensor;

  int CHI;
  int LX;
  int LY;
  int N_UNIT;

  std::string outdir;

  double time_simple_update;
  double time_full_update;
  double time_environment;
  double time_observable;
};

template <class ptensor>
TeNeS<ptensor>::TeNeS(MPI_Comm comm_, PEPS_Parameters peps_parameters_,
                      Lattice lattice_, NNOperators<ptensor> simple_updates_,
                      NNOperators<ptensor> full_updates_,
                      Operators<ptensor> onesite_operators_,
                      Operators<ptensor> twosite_operators_,
                      CorrelationParameter corparam_)
    : comm(comm_), peps_parameters(peps_parameters_), lattice(lattice_),
      simple_updates(simple_updates_), full_updates(full_updates_),
      onesite_operators(onesite_operators_),
      twosite_operators(twosite_operators_), corparam(corparam_),
      outdir("output"), time_simple_update(), time_full_update(),
      time_environment(), time_observable() {

  MPI_Comm_size(comm, &mpisize);
  MPI_Comm_rank(comm, &mpirank);

  peps_parameters.Bcast(comm);
  // output debug or warning info only from process 0
  if (mpirank != 0) {
    peps_parameters.print_level = PEPS_Parameters::PrintLevel::none;
  }

  CHI = peps_parameters.CHI;

  lattice.Bcast(comm);

  LX = lattice.LX;
  LY = lattice.LY;
  N_UNIT = lattice.N_UNIT;

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

  initialize_tensors();

  int maxops = -1;
  for (auto const &op : onesite_operators) {
    maxops = std::max(op.group, maxops);
  }
  num_onesite_operators = maxops + 1;

  maxops = -1;
  for (auto const &op : twosite_operators) {
    maxops = std::max(op.group, maxops);
  }
  num_twosite_operators = maxops + 1;

  site_ops_indices.resize(N_UNIT, std::vector<int>(num_onesite_operators, -1));
  for(int i=0; i<onesite_operators.size(); ++i){
    auto const& op = onesite_operators[i];
    site_ops_indices[op.source_site][op.group] = i;
  }
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
    const auto pdim = lattice.physical_dims[i];
    const auto vdim = lattice.virtual_dims[i];

    Tn.push_back(ptensor(Shape(vdim[0], vdim[1], vdim[2], vdim[3], pdim)));
    eTt.push_back(ptensor(Shape(CHI, CHI, vdim[1], vdim[1])));
    eTr.push_back(ptensor(Shape(CHI, CHI, vdim[2], vdim[2])));
    eTb.push_back(ptensor(Shape(CHI, CHI, vdim[3], vdim[3])));
    eTl.push_back(ptensor(Shape(CHI, CHI, vdim[0], vdim[0])));
    C1.push_back(ptensor(Shape(CHI, CHI)));
    C2.push_back(ptensor(Shape(CHI, CHI)));
    C3.push_back(ptensor(Shape(CHI, CHI)));
    C4.push_back(ptensor(Shape(CHI, CHI)));

    std::vector<std::vector<double>> lambda(nleg);
    for (int j = 0; j < nleg; ++j) {
      lambda[j] = std::vector<double>(vdim[j], 1.0);
    }
    lambda_tensor.push_back(lambda);

    ptensor id(mptensor::Shape(pdim, pdim));
    for (int j = 0; j < pdim; ++j) {
      for (int k = 0; k < pdim; ++k) {
        id.set_value(mptensor::Index(j, k), (j == k ? 1.0 : 0.0));
      }
    }
    op_identity.push_back(id);
  }

  std::mt19937 gen(peps_parameters.seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  int nr;

  std::string const &load_dir = peps_parameters.tensor_load_dir;
  if (load_dir.empty()) {
    Index index;
    for (int i = 0; i < lattice.N_UNIT; ++i) {
      const auto pdim = lattice.physical_dims[i];
      const auto vdim = lattice.virtual_dims[i];

      const size_t ndim = vdim[0] * vdim[1] * vdim[2] * vdim[3] * pdim;
      std::vector<double> ran(ndim);

      for (int j = 0; j < ndim; j++) {
        ran[j] = dist(gen);
      }
      auto &dir = lattice.initial_dirs[i];
      if (std::all_of(dir.begin(), dir.end(),
                      [=](double x) { return x == 0.0; })) {
        // random
        dir.resize(pdim);
        for (int j = 0; j < pdim; ++j) {
          dir[j] = dist(gen);
        }
      }
      for (int n = 0; n < Tn[i].local_size(); ++n) {
        index = Tn[i].global_index(n);
        if (index[0] == 0 && index[1] == 0 && index[2] == 0 && index[3] == 0) {
          Tn[i].set_value(index, dir[index[4]]);
        } else {
          nr = index[0] + index[1] * vdim[0] + index[2] * vdim[0] * vdim[1] +
               index[3] * vdim[0] * vdim[1] * vdim[2] +
               index[4] * vdim[0] * vdim[1] * vdim[2] * vdim[3];
          Tn[i].set_value(index, lattice.noises[i] * ran[nr]);
        }
      }
    }
  } else {
    // load from the checkpoint
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
    std::vector<double> ls;
    if (mpirank == 0) {
      for (int i = 0; i < N_UNIT; ++i) {
        const auto vdim = lattice.virtual_dims[i];
        std::ifstream ifs(load_dir + "/lambda_" + std::to_string(i) + ".dat");
        for (int j = 0; j < nleg; ++j) {
          for (int k = 0; k < vdim[j]; ++k) {
            double temp = 0.0;
            ifs >> temp;
            ls.push_back(temp);
          }
        }
      }
    }
    int lsize = ls.size();
    MPI_Bcast(&lsize, 1, MPI_INT, 0, comm);
    if (mpirank != 0) {
      ls.resize(lsize);
    }
    MPI_Bcast(&(ls[0]), lsize, MPI_DOUBLE, 0, comm);
    int index = 0;
    for (int i = 0; i < N_UNIT; ++i) {
      const auto vdim = lattice.virtual_dims[i];
      for (int j = 0; j < nleg; ++j) {
        for (int k = 0; k < vdim[j]; ++k) {
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
std::vector<std::vector<double>> TeNeS<ptensor>::measure_onesite(bool save) {
  Timer<> timer;
  const int nlops = num_onesite_operators;
  std::vector<std::vector<double>> local_obs(
      nlops,
      std::vector<double>(N_UNIT, std::numeric_limits<double>::quiet_NaN()));

  std::vector<double> norm(N_UNIT);
  for (int i = 0; i < N_UNIT; ++i) {
    norm[i] = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                eTb[i], eTl[i], Tn[i], op_identity[i]);
  }
  for (auto const &op : onesite_operators) {
    const int i = op.source_site;
    double val = Contract_one_site(C1[i], C2[i], C3[i], C4[i], eTt[i], eTr[i],
                                   eTb[i], eTl[i], Tn[i], op.op);
    local_obs[op.group][i] = val / norm[i];
  }

  time_observable += timer.elapsed();

  if (save && mpirank == 0) {
    std::string filename = outdir + "/onesite_obs.dat";
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save onesite observables to " << filename << std::endl;
    }
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << "# $1: op_group\n";
    ofs << "# $2: site_index (-1: mean, -2: sum) \n";
    ofs << "# $3: real\n";
    ofs << "# $4: imag\n";
    ofs << std::endl;

    for (int ilops = 0; ilops < nlops; ++ilops) {
      int num = 0;
      double sum = 0.0;
      for (int i = 0; i < N_UNIT; ++i) {
        if (std::isnan(local_obs[ilops][i])) {
          continue;
        }
        num += 1;
        sum += local_obs[ilops][i];
        ofs << ilops << " " << i << " " << local_obs[ilops][i] << " " << 0.0
            << std::endl;
      }
      ofs << ilops << " " << -1 << " " << sum / num << " " << 0.0 << std::endl;
      ofs << ilops << " " << -2 << " " << sum << " " << 0.0 << std::endl;
    }
  }
  return local_obs;
}

template <class ptensor>
std::vector<std::map<Bond, double>> TeNeS<ptensor>::measure_twosite(bool save) {
  Timer<> timer;

  const int nlops = num_twosite_operators;
  std::vector<std::map<Bond, double>> ret(nlops);

  constexpr int nmax = 4;

  // norm[{site, dx, dy}]
  std::map<std::tuple<int, int, int>, double> norms;
  // std::vector<std::vector<std::vector<double>>> norms(
  //    N_UNIT, std::vector<double>(d_max, std::vectorstd::numeric_limits<double>::quiet_NaN()));

  for (const auto &op : twosite_operators) {
    const int source = op.source_site;
    const int x_source = lattice.x(source);
    const int y_source = lattice.y(source);
    const int target = op.target_site;
    const int x_target =
        lattice.x(target) + op.offset_x * lattice.LX + op.offset_y * lattice.skew;
    const int y_target = lattice.y(target) + op.offset_y * lattice.LY;
    const int dx = x_target - x_source;
    const int dy = y_target - y_source;

    const int ncol = std::abs(dx)+1;
    const int nrow = std::abs(dy)+1;
    if (ncol > nmax || nrow > nmax) {
      std::cerr << "Warning: now version of TeNeS does not support too long "
                   "operator"
                << std::endl;
      std::cerr << "group = " << op.group << " (dx = " << dx << ", dy = " << dy
                << ")" << std::endl;
      continue;
    }

    std::vector<const ptensor*> C_(4, nullptr);
    std::vector<const ptensor*> eTt_(ncol, nullptr);
    std::vector<const ptensor*> eTr_(nrow, nullptr);
    std::vector<const ptensor*> eTb_(ncol, nullptr);
    std::vector<const ptensor*> eTl_(nrow, nullptr);

    /*
     * Caution: orders of tensors in unitcell and Contract_* function are different
     *
     * Lattice:
     *
     *    y
     *    ^
     *    |
     *    0--> x
     *
     * Contract_*:
     *    0-->col
     *    |
     *    v
     *    row
     *
     */
    std::vector<std::vector<const ptensor*>> Tn_(nrow, std::vector<const ptensor*>(ncol, nullptr));
    std::vector<std::vector<const ptensor*>> op_(nrow, std::vector<const ptensor*>(ncol, nullptr));

    std::vector<std::vector<int>> indices(nrow, std::vector<int>(ncol));

    int source_col, source_row, target_col, target_row;

    if(dx >= 0){
      source_col = 0;
      target_col = ncol-1;
    }else{
      source_col = ncol-1;
      target_col = 0;
    }
    if(dy >= 0){
      source_row = nrow-1;
      target_row = 0;
    }else{
      source_row = 0;
      target_row = nrow-1;
    }

    for(int row=0; row<nrow; ++row){
      for(int col=0; col<ncol; ++col){
        const int index = lattice.other(source, col-source_col, source_row-row);
        indices[row][col] = index;
        op_[row][col] = &(op_identity[index]);
        Tn_[row][col] = &(Tn[index]);
      }
      eTl_[row] = &(eTl[indices[row][0]]);
      eTr_[row] = &(eTr[indices[row][ncol-1]]);
    }
    for(int col=0; col<ncol; ++col){
      eTt_[col] = &(eTt[indices[0][col]]);
      eTb_[col] = &(eTb[indices[nrow-1][col]]);
    }
    C_[0] = &(C1[indices[0][0]]);
    C_[1] = &(C2[indices[0][ncol-1]]);
    C_[2] = &(C3[indices[nrow-1][ncol-1]]);
    C_[3] = &(C4[indices[nrow-1][0]]);

    const auto norm_key = std::make_tuple(indices[0][0], nrow, ncol);
    auto norm = (norms.count(norm_key) ? norms[norm_key] : std::numeric_limits<double>::quiet_NaN() );
    if(std::isnan(norm)){
      norm = Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      norms[norm_key] = norm;
    }

    double value = 0.0;
    if(op.ops_indices.empty()){
      ptensor U, VT;
      std::vector<double> s;
      mptensor::svd(op.op, {0, 2}, {1, 3}, U, s, VT);
      const int ns = s.size();
      for(int is=0; is<ns; ++is){
        ptensor source_op = reshape(slice(U, 2, is, is + 1), {U.shape()[0], U.shape()[0]});
        op_[source_row][source_col] = &source_op;
        ptensor target_op = reshape(slice(VT, 0, is, is + 1), {VT.shape()[1], VT.shape()[1]});
        op_[target_row][target_col] = &target_op;
        auto localvalue = Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
        value += localvalue * s[is];
      }
    }else{
      op_[source_row][source_col] = &(onesite_operators[siteoperator_index(op.source_site, op.ops_indices[0])].op);
      op_[target_row][target_col] = &(onesite_operators[siteoperator_index(op.target_site, op.ops_indices[1])].op);
      auto localvalue = Contract(C_, eTt_, eTr_, eTb_, eTl_, Tn_, op_);
      value += localvalue;
    }
    ret[op.group][{op.source_site, op.target_site, op.offset_x,
                   op.offset_y}] = value / norm;
  }

  time_observable += timer.elapsed();

  if (save && mpirank == 0) {
    std::string filename = outdir + "/twosite_obs.dat";
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "    Save twosite observables to " << filename << std::endl;
    }
    std::ofstream ofs(filename.c_str());
    ofs << std::scientific
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << "# $1: op_group\n";
    ofs << "# $2: source_site (-1: mean, -2: sum)\n";
    ofs << "# $3: target_site\n";
    ofs << "# $4: offset_x\n";
    ofs << "# $5: offset_y\n";
    ofs << "# $6: real\n";
    ofs << "# $7: imag\n";
    ofs << std::endl;
    for (int ilops = 0; ilops < nlops; ++ilops) {
      double sum = 0.0;
      int num = 0;
      for (const auto &r : ret[ilops]) {
        auto bond = r.first;
        auto value = r.second;
        sum += value;
        num += 1;
        ofs << ilops << " " << bond.source_site << " " << bond.target_site
            << " " << bond.offset_x << " " << bond.offset_y << " " << value
            << " " << 0.0 << std::endl;
      }
      ofs << ilops << " " << -1 << " 0 0 0 " << sum / num << " " << 0.0
          << std::endl;
      ofs << ilops << " " << -2 << " 0 0 0 " << sum << " " << 0.0 << std::endl;
    }
  }

  return ret;
}

/*
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

  std::vector<Correlation> correlations;
  for (int left_index = 0; left_index < N_UNIT; ++left_index) {
    const auto vdim = lattice.virtual_dims[left_index];
    ptensor correlation_T(Shape(CHI, CHI, vdim[0], vdim[0]));
    ptensor correlation_norm(Shape(CHI, CHI, vdim[0], vdim[0]));
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
                         Tn[left_index], op_identity[left_index]);

        for (int r = 0; r < r_max; ++r) {
          const int right_x = (left_x + r + 1) % LX;
          const int right_y = left_y;
          const int offset_x = (left_x + r + 1) / LX;
          const int offset_y = 0;
          const int right_index = lattice.index(right_x, right_y);
          double norm = FinishCorrelation(
              correlation_norm, C2[right_index], C3[right_index],
              eTt[right_index], eTr[right_index], eTb[right_index],
              Tn[right_index], op_identity[right_index]);
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
                         op_identity[left_index]);

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
                                          tn, op_identity[right_index]);
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
*/

template <class ptensor> void TeNeS<ptensor>::measure() {
  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "Start calculating observables" << std::endl;
    std::clog << "  Start updating environment" << std::endl;
  }
  update_CTM();

  // const int nlops = lops.size();

  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "  Start calculating local operators" << std::endl;
  }
  auto local_obs = measure_onesite(true);

  if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
    std::clog << "  Start calculating NN correlation" << std::endl;
  }
  auto NN_obs = measure_twosite(true);
  auto energy = 0.0;
  for (const auto &nn : NN_obs[0]) {
    energy += nn.second;
  }

  /*
  if (corparam.r_max > 0) {
    if (peps_parameters.print_level >= PEPS_Parameters::PrintLevel::info) {
      std::clog << "  Start calculating long range correlation" << std::endl;
    }
    auto correlations = measure_correlation(true);
  }
  */

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

      for (int ilops = 0; ilops < num_onesite_operators; ++ilops) {
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
        for (int k = 0; k < lattice.virtual_dims[i][j]; ++k) {
          ofs << lambda_tensor[i][j][k] << std::endl;
        }
      }
    }
  }
}

template <class tensor>
int tenes(MPI_Comm comm, PEPS_Parameters peps_parameters, Lattice lattice,
          NNOperators<tensor> simple_updates, NNOperators<tensor> full_updates,
          Operators<tensor> onesite_operators,
          Operators<tensor> twosite_operators, CorrelationParameter corparam) {
  TeNeS<tensor> tns(comm, peps_parameters, lattice, simple_updates,
                    full_updates, onesite_operators, twosite_operators,
                    corparam);
  tns.optimize();
  tns.save_tensors();
  tns.measure();
  return 0;
}

// template specialization
using d_tensor = mptensor::Tensor<mptensor_matrix_type, double>;
template int tenes<d_tensor>(MPI_Comm comm, PEPS_Parameters peps_parameters,
                             Lattice lattice,
                             NNOperators<d_tensor> simple_updates,
                             NNOperators<d_tensor> full_updates,
                             Operators<d_tensor> onesite_operators,
                             Operators<d_tensor> twosite_operators,
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
