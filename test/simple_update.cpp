#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <vector>
#include <fstream>

#include <PEPS_Basics.hpp>
#include <PEPS_Parameters.cpp>

TEST_CASE("testing simple update"){
  using tensor = mptensor::Tensor<mptensor::scalapack::Matrix, double>;

  const int ldof = 2;
  const int D = 2;
  const int nleg = 4;

  const double tol = 1.0e-8;

  // make inputs

  tensor Tn1(Shape(D, D, D, D, ldof));
  tensor Tn2(Shape(D, D, D, D, ldof));
  std::vector<std::vector<double>> lambda_1(nleg, std::vector<double>(D, 1.0));
  std::vector<std::vector<double>> lambda_2(nleg, std::vector<double>(D, 1.0));

  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j)
  for(int k=0; k<D; ++k)
  for(int l=0; l<D; ++l)
  for(int m=0; m<ldof; ++m)
  {
    Tn1.set_value(Index(i,j,k,l,m), 1.0);
    Tn2.set_value(Index(i,j,k,l,m), 1.0);
  }

  tensor op(Shape(ldof, ldof, ldof, ldof));

  for(int i=0; i<ldof; ++i)
  for(int j=0; j<ldof; ++j)
  for(int k=0; k<ldof; ++k)
  for(int l=0; l<ldof; ++l)
  {
    op.set_value(Index(i,j,k,l), 0.0);
  }
  for(int i=0; i<ldof; ++i)
  {
    op.set_value(Index(i,i,i,i), 1.0);
  }

  // load answer

  std::ifstream ifs("data/simple_update.dat");

  tensor ans_Tn1(Shape(D, D, D, D, ldof));
  tensor ans_Tn2(Shape(D, D, D, D, ldof));
  std::vector<double> ans_lambda(D, 0.0);

  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j)
  for(int k=0; k<D; ++k)
  for(int l=0; l<D; ++l)
  for(int m=0; m<ldof; ++m)
  {
    double val;
    ifs >> val;
    ans_Tn1.set_value(Index(i,j,k,l,m), val);
  }
  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j)
  for(int k=0; k<D; ++k)
  for(int l=0; l<D; ++l)
  for(int m=0; m<ldof; ++m)
  {
    double val;
    ifs >> val;
    ans_Tn2.set_value(Index(i,j,k,l,m), val);
  }
  for(int i=0; i<D; ++i){
    double val;
    ifs >> val;
    ans_lambda[i] = val;
  }

  // calculation

  PEPS_Parameters peps_parameters;
  int connect = 2;
  tensor new_Tn1, new_Tn2;
  std::vector<double> new_lambda;

  Simple_update_bond(Tn1, Tn2,
      lambda_1, lambda_2,
      op, connect, peps_parameters,
      new_Tn1, new_Tn2, new_lambda
      );


  // check results

  std::ofstream ofs("res_simple_update.dat");
  ofs << std::setprecision(std::numeric_limits<double>::digits10);

  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j)
  for(int k=0; k<D; ++k)
  for(int l=0; l<D; ++l)
  for(int m=0; m<ldof; ++m)
  {
    double result, answer;
    new_Tn1.get_value(Index(i,j,k,l,m), result);
    ans_Tn1.get_value(Index(i,j,k,l,m), answer);
    CHECK(result == doctest::Approx(answer).epsilon(tol));
    ofs << result << " ";
  }
  ofs << std::endl;

  for(int i=0; i<D; ++i)
  for(int j=0; j<D; ++j)
  for(int k=0; k<D; ++k)
  for(int l=0; l<D; ++l)
  for(int m=0; m<ldof; ++m)
  {
    double result, answer;
    new_Tn2.get_value(Index(i,j,k,l,m), result);
    ans_Tn2.get_value(Index(i,j,k,l,m), answer);
    CHECK(result == doctest::Approx(answer).epsilon(tol));
    ofs << result << " ";
  }
  ofs << std::endl;

  for(int i=0; i<D; ++i){
    double result = new_lambda[i];
    double answer = ans_lambda[i];
    CHECK(result == doctest::Approx(answer).epsilon(tol));
    ofs << new_lambda[i] << " ";
  }
  ofs << std::endl;
}

