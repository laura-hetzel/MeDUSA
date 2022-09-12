#include <RcppArmadillo.h>
#include <armadillo>
using namespace Rcpp;
using namespace arma;

typedef std::vector<double> stdvec;
typedef std::vector< std::vector<double> > stdvecvec;

// [[Rcpp::depends(RcppArmadillo)]]

NumericVector diffC(NumericVector intensities){
    int n = intensities.size() - 1;
    NumericVector res(n);
    for (int i = 0; i < n; i++) {
        res[i] = intensities[i + 1] - intensities[i];
    }
    return res;
}

NumericVector signC(NumericVector differences){
    int n =  differences.size();
    NumericVector res(n);
    for (int i = 0; i < n; i++) {
        res[i] = (differences[i] > 0) - (differences[i] < 0);
    }
    return res;
}



template<class Object>
arma::mat elementwise_pow(const Object& base, const Object& p) {
    assert(base.n_elem == p.n_elem);
    Object result;
    result.copy_size(base);
    for (std::size_t i = 0; i < result.n_elem; ++i) {
        result[i] = std::pow(base[i], p[i]);
    }
    return (mat) result;
}

arma::mat savgolCoefC(int m = 2, int nk = 3){
    nk = nk + 1;
    int nm = 2 * m + 1;

    mat K(nm, nk);
    for (int i = 0; i < nk; i++){
        for (int j = 0; j < nm; j++){
            K.col(i)[j] = i;
        }
    }

    mat F(nm, nm);

    for (int i = 0; i < m + 1; i++){

        mat M(nm, nk);
        for (int j = 0; j < nk; j++){
            for (int k = 0; k < nm; k++){
                M.col(j)[k] = k - i;
            }
        }
        mat X = elementwise_pow(M, K);
        mat Tr = X.t();
        mat T = solve(Tr * X, Tr);
        F.row(i) = T.row(0);
    }

    int rows = nm - (m + 1);
    F.tail_rows(rows) = reverse(reverse(F.head_rows(rows), 1));
    return F;
}

arma::vec savgolFilterC(NumericVector vec_rcpp, int hws = 2){
    mat coef = savgolCoefC(hws);

    vec x = as<vec>(wrap(vec_rcpp));
    vec y = conv(x, coef.row(hws), "same");


    int w = 2 * hws + 1;
    y.head(hws) = coef.head_rows(hws) * x.head(w);
    y.tail(hws) = coef.tail_rows(hws) * x.tail(w);
    return y;
}

// [[Rcpp::export]]
std::vector<int> getLocalMaximaC(NumericVector intensities){
    intensities = savgolFilterC(intensities);
    int n = intensities.size();
    for (int i = 0; i < n; i++){
        if (intensities[i] < 0){
            intensities[i] = 0;
        }
    }

    NumericVector maximas = diffC(signC(diffC(intensities)));
    std::vector<int> res;
    int count = 0;
    for (int i = 0; i < n; i++){
        if (maximas[i] == -2) {
            res.push_back(i + 2); // + 2 since C starts with 0, while R starts with 1
            count += 1;
        }
    }
    return res;
}
