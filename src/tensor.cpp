#include "tensor.hpp"

#include <sstream>
#include <numeric>
#include "util/abs.hpp"

using dcomplex = std::complex<double>;

// LAPACK routine
extern "C" {
void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* wr,
            double* wi, double* vl, int* ldvl, double* vr, int* ldvr,
            double* work, int* lwork, int* info);
void zgeev_(char* jobvl, char* jobvr, int* n, dcomplex* a, int* lda,
            dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr,
            dcomplex* work, int* lwork, double* rwork, int* info);
}

namespace tenes {

void eigen_impl(small_tensor<double> const& A, std::vector<dcomplex>& eigvals,
                std::vector<dcomplex>& eigvecs_last) {
  mptensor::lapack::Matrix<double> a = A.get_matrix();
  char jobl = 'N';
  char jobr = 'V';
  int ld = A.shape()[0];
  int dim = ld;
  int N = ld * dim;
  std::vector<double> wr(dim), wi(dim);
  std::vector<double> vl(1), vr(N);
  int lwork = -1;
  std::vector<double> work(1);
  int info = 0;

  dgeev_(&jobl, &jobr, &dim, a.head(), &ld, wr.data(), wi.data(), vl.data(),
         &ld, vr.data(), &ld, work.data(), &lwork, &info);
  lwork = static_cast<int>(work[0]);
  work.resize(lwork);
  dgeev_(&jobl, &jobr, &dim, a.head(), &ld, wr.data(), wi.data(), vl.data(),
         &ld, vr.data(), &ld, work.data(), &lwork, &info);

  eigvals.resize(dim);
  eigvecs_last.resize(dim);
  for (int i = 0; i < dim; ++i) {
    eigvals[i] = std::complex<double>(wr[i], wi[i]);
  }
  for (int i = 0; i < dim; ++i) {
    if (wi[i] == 0.0) {
      eigvecs_last[i] = vr[ld - 1 + i * ld];
    } else {
      double re = vr[ld - 1 + i * ld];
      double im = vr[ld - 1 + (i + 1) * ld];
      eigvecs_last[i] = std::complex<double>(re, im);
      eigvecs_last[i + 1] = std::complex<double>(re, -im);
      ++i;
    }
  }
}

void eigen_impl(small_tensor<dcomplex> const& A, std::vector<dcomplex>& eigvals,
                std::vector<dcomplex>& eigvecs_last) {
  mptensor::lapack::Matrix<dcomplex> a = A.get_matrix();
  char jobl = 'N';
  char jobr = 'V';
  int ld = A.shape()[0];
  int dim = ld;
  int N = ld * dim;
  std::vector<dcomplex> vl(1), vr(N);
  int lwork = -1;
  std::vector<dcomplex> work(1);
  std::vector<double> rwork(2 * dim);
  int info = 0;

  eigvals.resize(dim);

  zgeev_(&jobl, &jobr, &dim, a.head(), &ld, eigvals.data(), vl.data(), &ld,
         vr.data(), &ld, work.data(), &lwork, rwork.data(), &info);
  lwork = static_cast<int>(work[0].real());
  work.resize(lwork);
  zgeev_(&jobl, &jobr, &dim, a.head(), &ld, eigvals.data(), vl.data(), &ld,
         vr.data(), &ld, work.data(), &lwork, rwork.data(), &info);

  eigvecs_last.resize(dim);
  for (int i = 0; i < dim; ++i) {
    eigvecs_last[i] = vr[i * ld + ld - 1];
  }
}

template <class value_type>
void eigen(small_tensor<value_type> const& A, std::vector<dcomplex>& eigvals,
           std::vector<dcomplex>& eigvecs_last, int nev) {
  int k = A.shape()[0];
  eigen_impl(A, eigvals, eigvecs_last);

  std::vector<int> sorted_indices(k);
  std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
  std::partial_sort(sorted_indices.begin(), sorted_indices.begin() + nev,
                    sorted_indices.end(), [&eigvals](int a, int b) {
                      return util::abs2(eigvals[a]) > util::abs2(eigvals[b]);
                    });

  std::vector<dcomplex> evals_ret(nev);
  std::vector<dcomplex> evecs_ret(nev);
  for (int i = 0; i < nev; ++i) {
    evals_ret[i] = eigvals[sorted_indices[i]];
    evecs_ret[i] = eigvecs_last[sorted_indices[i]];
  }
  eigvals.swap(evals_ret);
  eigvecs_last.swap(evecs_ret);
}

template <class tensor>
tensor resize_tensor(tensor const& src, mptensor::Shape target_shape) {
  mptensor::Shape shape = src.shape();
  const size_t ndim = shape.size();
  if (target_shape.size() != ndim) {
    std::stringstream ss;
    ss << "dimension mismatch in resize_tensor: source = " << ndim
       << ", target = " << target_shape.size();
    tenes::logic_error(ss.str());
  }
  mptensor::Shape zero = shape;

  bool to_extend = false;
  bool to_shrink = false;

  for (size_t i = 0; i < ndim; ++i) {
    if (shape[i] < target_shape[i]) {
      to_extend = true;
      shape[i] = target_shape[i];
    } else if (shape[i] > target_shape[i]) {
      to_shrink = true;
    }
    zero[i] = 0;
  }
  tensor A;
  if (to_extend) {
    A = mptensor::extend(src, shape);
  } else {
    A = src;
  }
  if (to_shrink) {
    return mptensor::slice(A, zero, target_shape);
  } else {
    return A;
  }
}

template real_tensor resize_tensor(
    real_tensor const& src, mptensor::Shape target_shape);

template complex_tensor resize_tensor(
    complex_tensor const& src, mptensor::Shape target_shape);

template void eigen(small_tensor<double> const& A,
                    std::vector<std::complex<double>>& eigvals,
                    std::vector<std::complex<double>>& eigvecs_last, int nev);

template void eigen(small_tensor<dcomplex> const& A,
                    std::vector<std::complex<double>>& eigvals,
                    std::vector<std::complex<double>>& eigvecs_last, int nev);

}  // namespace tenes
