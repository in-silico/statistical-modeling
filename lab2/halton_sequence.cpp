
#include <bits/stdc++.h>
#include <kmeans.h>

using namespace std;
using namespace LibAnn;

template<class T>
ostream& operator<<(ostream& a, const vector<T>& v) {
  a << "{";
  if (v.size()>0) a << v[0];
  for (int i=1; i<v.size(); i++) a << ", " << v[i];
  a << "}";
  return a;
}

template<typename T>
struct Fraction {
	T p, q;

	Fraction() : p(0), q(1) {}
	Fraction(T P) : p(P), q(1) {}
	Fraction(T P, T Q) : p(P), q(Q) { simplify(); }
	void simplify() {
		T g = __gcd(p, q);
		p /= g;
		q /= g;
	}
	Fraction operator+(const Fraction &f) const {
		return Fraction(p * f.q + q * f.p, q * f.q);
	}
	Fraction operator-(const Fraction &f) const {
		return Fraction(p * f.q - q * f.p, q * f.q);
	}
	Fraction operator*(const Fraction &f) const {
		return Fraction(p * f.p, q * f.q);
	}
	Fraction operator/(const Fraction &f) const {
		return Fraction(p * f.q, q * f.p);
	}
	Fraction operator%(int m) const {
		return Fraction(p % (m*q), q);
	}
	Fraction operator-() const { return Fraction(-p, q); }
	bool operator<(const Fraction &f) const { return p*f.q < q*f.p; }
	bool operator>(const Fraction &f) const { return p*f.q > q*f.p; }
	bool operator<=(const Fraction &f) const { return p*f.q <= q*f.p; }
	bool operator>=(const Fraction &f) const { return p*f.q >= q*f.p; }
	bool operator==(const Fraction &f) const { return p == f.p && q == f.q; }
};

template <class T>
ostream& operator<<(ostream &a, const Fraction<T>& f) {
  a << f.p << "/" << f.q;
  return a;
}

/**** Code of Halton *****/

typedef Fraction<long long int> frac;

vector<int> genPrimes(int up_to) {
  vector<bool> isPrime((up_to-3)/2 + 1, true);
  vector<int> ans;
  ans.push_back(2);
  for (int p=3; p<=up_to; p+=2) {
    if (isPrime[(p-3)/2]) {
      ans.push_back(p);
      for (int x=p*p; x<=up_to; x += 2*p) {
        isPrime[(x-3)/2] = false;
      }
    }
  }
  return ans;
}

class Halton {
  int ndims;
  vector<int> primes;
public:
  Halton(int ndims) {
    primes = genPrimes(ndims*10);
    this->ndims = ndims;
  }

  vector< vector<frac> > genSequence(int N) {
    vector< vector<frac> > ans(N);
    for (int d = 0; d < ndims; d++) {
      int prime = primes[d];
      for (int index=1; index <= N; index++) {
        frac f(1, prime);
        frac result(0,1);
        for (int i=index; i>0; i /= prime) {
          result = result + f*(i % prime);
          f = f / prime;
        }
        ans[index-1].push_back(result);
      }
    }
    return ans;
  }
};

double toDouble(const frac& f) {
  return ((double)f.p)/f.q;
}

void generateVoronoi(int nclusters, int nsamples, int ndims) {
  Halton h(ndims);
  Mat* pts = new Mat(ndims*nsamples);
  pts->setSize(nsamples, ndims);
  vector< vector<frac> > samples = h.genSequence(nsamples);
  /*vector< vector<double> > corrected;
  for (int i=0; i<nsamples; i++) {
    double wl = toDouble(samples[i][0])*2 + 4;
    double ww = toDouble(samples[i][1])*2 + 3.5;
    double tl = toDouble(samples[i][2])*3 + 6;
    double al = toDouble(samples[i][3])*2 + 10;
    if ((
  }*/
  for (int i=0; i<nsamples; i++) {
    for (int j = 0; j < ndims; ++j){
      double xsample = ((double)samples[i][j].p) / samples[i][j].q;
      pts->get(i,j) = xsample;
    }
  }

  KMeansConf kf;
  kf.maxIter = 1000;
  kf.nthreads = 1;
  Mat* c = new Mat(nclusters*ndims);
  c->setSize(nclusters, ndims);
  kmeansInit(c, pts, nclusters);
  kmeans(c, pts, c, &kf);

  //print all the centers
  for (int i=0; i<nclusters; i++) {
    double wl = c->get(i,0)*2 + 4;
    double ww = c->get(i,1)*2 + 3.5;
    double tl = c->get(i,2)*3 + 6;
    double al = c->get(i,3)*2 + 10;
    /*for (int j=0; j<ndims; j++) {
      if (j!=0) cout << " ";
      cout << c->get(i,j);
    }
    cout << endl;*/
    cout << wl << " " << ww << " " << tl << " " << al << endl;
  }
  delete pts;
}

int main() {
  int n, d, nclusters;
  cin >> n >> d >> nclusters;
  /*Halton h(d);
  vector<vector<frac> > f = h.genSequence(n);
  for(int i=0; i<n; i++) {
    for (int j=0; j<d; j++) {
      if (j!=0) cout << " ";
      double x = ((double)f[i][j].p) / f[i][j].q;
      cout << x;
    }
    cout << endl;
  }*/
  generateVoronoi(nclusters, n, d);
}
