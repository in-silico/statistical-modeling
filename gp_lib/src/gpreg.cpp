
#include "gplib.h"

using namespace arma;
using namespace std;

namespace gplib {

  struct GPReg::Implementation {
    shared_ptr<Kernel> kernel;
    mat X; //Matrix of inputs
    vec y; //vector of outputs

    vec evalMean(const arma::mat& data) {
      // For the moment just use the zero mean
      return zeros<mat>(data.n_rows, data.n_rows);
    }

    MVGauss predict(const arma::mat& newData) {
      mat M = join_vert(X, newData);
      int N = X.n_rows, Nval = newData.n_rows;
      mat cov = kernel->eval(M,M);
      vec mean = evalMean(M);
      MVGauss gd(mean, cov);
      vector<bool> observed(N+Nval, false);
      for (int i=0; i<N; i++) observed[i] = true;
      return gd.conditional(y, observed);
    }
  };

  GPReg::GPReg() {
    pimpl = new Implementation();
  }

  GPReg::~GPReg() {
    delete pimpl;
  }

  void GPReg::setKernel(const std::shared_ptr<Kernel>& k) {
    pimpl->kernel = k;
  }

  shared_ptr<Kernel> GPReg::getKernel() const {
    return pimpl->kernel;
  }

  void GPReg::setTrainingSet(const arma::mat &X, const arma::vec& y) {
    pimpl->X = X;
    pimpl->y = y;
  }

  //void train();

  arma::vec GPReg::predict(const arma::mat& newData) const {
    MVGauss g = pimpl->predict(newData);
    return g.getMean();
  }

};
