#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <vector>
#include <fstream>

#include <util/string.cpp>
#include <Lattice.cpp>
#include <edge.cpp>
#include <PEPS_Parameters.cpp>
#include <load_toml.cpp>

TEST_CASE("input"){
  using ptensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;

  auto input_toml = cpptoml::parse_file("data/check_input.toml");

  SUBCASE("parameter"){
    PEPS_Parameters peps_parameters;
    peps_parameters.set(input_toml->get_table("parameter"));

    CHECK(peps_parameters.D == 4);
    CHECK(peps_parameters.CHI == 16);
    
    CHECK(peps_parameters.num_simple_step == 1000);
    CHECK(peps_parameters.Inverse_lambda_cut == 1e-10);

    CHECK(peps_parameters.num_full_step == 1);
    CHECK(peps_parameters.Inverse_Env_cut == 1e-10);
    CHECK(peps_parameters.Full_Inverse_precision == 1e-10);
    CHECK(peps_parameters.Full_Convergence_Epsilon == 1e-10);
    CHECK(peps_parameters.Full_max_iteration == 100);
    CHECK(peps_parameters.Full_Gauge_Fix == false);
    CHECK(peps_parameters.Full_Use_FastFullUpdate == false);

    CHECK(peps_parameters.Inverse_projector_cut == 1e-10);
    CHECK(peps_parameters.CTM_Convergence_Epsilon == 1e-8);
    CHECK(peps_parameters.Max_CTM_Iteration == 10);
    CHECK(peps_parameters.CTM_Projector_corner == true);
    CHECK(peps_parameters.Use_RSVD == true);
    CHECK(peps_parameters.RSVD_Oversampling_factor == 10);

    CHECK(peps_parameters.Lcor == 5);
  }

  SUBCASE("lattice"){
    auto toml_lattice = input_toml->get_table("lattice");
    Lattice lattice = gen_lattice(toml_lattice);
    CHECK(lattice.LX == 3);
    CHECK(lattice.LY == 2);
  }

  SUBCASE("evolution"){
    auto toml_evolution= input_toml->get_table("evolution");

    {
      INFO("simple_update");
      const auto simple_edges = gen_edges(toml_evolution, "simple_update", "evolution");
      CHECK(simple_edges[0].dir == Edge::horizontal);
      CHECK(simple_edges[0].source_site == 0);
      CHECK(simple_edges[0].target_site == 1);
      CHECK(simple_edges[0].op_id == 0);

      CHECK(simple_edges[1].dir == Edge::vertical);
      CHECK(simple_edges[1].source_site == 0);
      CHECK(simple_edges[1].target_site == 1);
      CHECK(simple_edges[1].op_id == 1);
    }
    {
      INFO("full_update");
      const auto full_edges = gen_edges(toml_evolution, "full_update", "evolution");
      CHECK(full_edges[0].dir == Edge::horizontal);
      CHECK(full_edges[0].source_site == 0);
      CHECK(full_edges[0].target_site == 2);
      CHECK(full_edges[0].op_id == 0);

      CHECK(full_edges[1].dir == Edge::vertical);
      CHECK(full_edges[1].source_site == 0);
      CHECK(full_edges[1].target_site == 2);
      CHECK(full_edges[1].op_id == 2);
    }

    {
      INFO("matrix");
      const auto evolutions = gen_matrices<ptensor>(toml_evolution, "matrix", "evolution");
      CHECK(evolutions[0].shape() == mptensor::Shape(4,4));
      for(int i=0; i<4; ++i){
        for(int j=0; j<4; ++j){
          double val;
          evolutions[0].get_value(mptensor::Index(i,j), val);
          if(i==j){
            CHECK(val == 1.0);
          }else{
            CHECK(val == 0.0);
          }
        }
      }
    }
  }

  SUBCASE("observable"){
    auto toml_observable = input_toml->get_table("observable");

    {
      INFO("local_operator");
      const auto lops = gen_matrices<ptensor>(toml_observable, "local_operator", "observable");
      CHECK(lops.size() == 2);
      CHECK(lops[0].shape() == mptensor::Shape(2,2));

      double val;
      lops[0].get_value(mptensor::Index(0,0), val);
      CHECK(val == 0.5);
      lops[0].get_value(mptensor::Index(0,1), val);
      CHECK(val == 0.0);
      lops[0].get_value(mptensor::Index(1,0), val);
      CHECK(val == 0.0);
      lops[0].get_value(mptensor::Index(1,1), val);
      CHECK(val == -0.5);

      lops[1].get_value(mptensor::Index(0,0), val);
      CHECK(val == 0.0);
      lops[1].get_value(mptensor::Index(0,1), val);
      CHECK(val == 0.5);
      lops[1].get_value(mptensor::Index(1,0), val);
      CHECK(val == 0.5);
      lops[1].get_value(mptensor::Index(1,1), val);
      CHECK(val == 0.0);
    }

    const auto hams = gen_matrices<ptensor>(toml_observable, "hamiltonian", "observable");
    {
      INFO("hamiltonian");
      const auto hams = gen_matrices<ptensor>(toml_observable, "hamiltonian", "observable");
      CHECK(hams.size() == 1);
      CHECK(hams[0].shape() == mptensor::Shape(4,4));

      double val;
      for(int i=0; i<4; ++i)
      for(int j=0; j<4; ++j){
        double val;
        hams[0].get_value(mptensor::Index(i,j), val);
        if(i==0 && j==0){
          CHECK(val == 0.25);
        }else if(i==1 && j==1){
          CHECK(val == -0.25);
        }else if(i==1 && j==2){
          CHECK(val == 0.5);
        }else if(i==2 && j==1){
          CHECK(val == 0.5);
        }else if(i==2 && j==2){
          CHECK(val == -0.25);
        }else if(i==3 && j==3){
          CHECK(val == 0.25);
        }else{
          CHECK(val == 0.0);
        }
      }
    }

    {
      INFO("hamiltonian_bonds");
      const auto ham_edges = gen_edges(toml_observable, "hamiltonian_bonds", "observable");
      CHECK(ham_edges[0].dir == Edge::horizontal);
      CHECK(ham_edges[0].source_site == 0);
      CHECK(ham_edges[0].target_site == 1);
      CHECK(ham_edges[0].op_id == 0);

      CHECK(ham_edges[1].dir == Edge::vertical);
      CHECK(ham_edges[1].source_site == 0);
      CHECK(ham_edges[1].target_site == 1);
      CHECK(ham_edges[1].op_id == 2);
    }

  }

}

