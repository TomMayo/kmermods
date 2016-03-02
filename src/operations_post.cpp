#include <Rcpp.h>
using namespace Rcpp;
//' Creates the update for the parameter vector according to the error term 
//' and the kmers involved
//'
//' Returns a vector of length equal to the parameter vector, which is update 
//' vector to be added to the parameter vector
//' 
//' @param kmers A vector of integers representing the kmers that were 
//' present leading to the error term
//' @param update_vec The parameter update vector so far
//' @param err_term The error term used for the update
//' @param err_term The error term used for the update
//' @return add_to_update TRUE/FALSE, whether the calculated update vector for 
//' the run is added automatically to the input update vector, so that the 
//' function returns a running total, which is the faster option.
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector update_paras(IntegerVector kmers, NumericVector update_vec,  double err_term,
                           bool add_to_update = true){
    NumericVector update;
    if (add_to_update==true) {
        update = clone(update_vec);
    } else {
        update = rep(0.0, update_vec.size());    
    }
    int n = kmers.size();
    for(int i = 0; i< n; i++) {
        int ind = kmers[i];
        update[ind] += err_term;
    }
    return update;
}

#include <Rcpp.h>
using namespace Rcpp;
typedef Rcpp::Nullable<Rcpp::NumericVector> nullable_t;


//' Computes the dot product of the parameters with the kmer counts, with or 
//' without it being warped
//'
//' Computes the dot product of the parameters with the kmer counts, with or 
//' without it being warped
//' 
//' @param kmers is a vector of integers of any length representing kmers over
//'  a window
//' @param paras is a vector of length equal to the total number of kmers
//' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry 
//' @return A double, representing the dot product of the parameters with what 
//' would usually be the vector of kmer abundances (warped or not)
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @examples 
//' kmers <- c(0,7,89,45,75,65,22,12)
//' paras <- rep(1,100)
//' kmer_dot_prod(kmers,paras)
//' @export
// [[Rcpp::export]]
double kmer_dot_prod_c(Rcpp::IntegerVector indices, Rcpp::NumericVector params, nullable_t warp_ = R_NilValue) {
    if (Rcpp::na_omit(indices).size() != indices.size()) {
        return NA_REAL;
    }
    if (warp_.isNull()) {
        Rcpp::NumericVector tmp = params[indices];
        return Rcpp::sum(tmp);    
    } else {
        Rcpp::NumericVector warp(warp_), tmp = params[indices];
        return Rcpp::sum(tmp * warp); 
    }
}

#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' Takes a longer kmer and converts it into a vector of its shorter kmers
//'
//' Takes in a single long kmer, for example a 26mer, and returns the composite 
//' shorter kmers for some k, for example, all the 10mers within the 26mer. 
//' It is designed to be used to unwrap a dense representation of the kmers: 
//' for example, given the 26mers every 19 base pairs, we can quickly extract 
//' all 8mers.
//' 
//' @param kmer A number representing a single large kmer, with k = old_len
//' @param old_len The length of the kmer represented by kmer (ie k)
//' @param new_len The length of the kmers we want to extract, eg, if we want 
//' 8mers, new_len = 8
//' @param base The length of the alphabet. For normal DNA sequence this is 4.
//' @param num_kmers The number of the new kmers to return. The functions 
//' returns the first num_kmers in the vector. Default value is 0, which returns 
//' all kmers.
//' @return A vector of doubles (actually integers represented as doubles), 
//' representing all the (new_len)-mers contained within the single (old_len)-mer.
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector unwrap_kmers(double kmer, double old_len = 26, double new_len = 8, 
                             double base = 4, int num_kmers = 0){
    if (num_kmers == 0){
        num_kmers = old_len - new_len + 1;    
    }
    NumericVector new_kmers(num_kmers); 
    for (int i = 0;i < num_kmers; i++){
        double top = fmod(kmer,  pow(base, old_len - i));
        double bottom = pow(base, old_len - new_len - i);
        new_kmers[i] = floor(top / bottom);
    }
    return new_kmers;
}

//' Takes in a single long kmer, for example a 26mer, and returns the composite 
//' shorter kmers for some k, for example, all the 10mers within the 26mer. 
//' It is designed to be used to unwrap a dense representation of the kmers: 
//' for example, given the 26mers every 19 base pairs, we can quickly extract 
//' all 8mers.
//' 
//' @param kmer A number representing a single large kmer, with k = old_len
//' @param old_len The length of the kmer represented by kmer (ie k)
//' @param new_len The length of the kmers we want to extract, eg, if we want 
//' 8mers, new_len = 8
//' @param base The length of the alphabet. For normal DNA sequence this is 4.
//' @param num_kmers The number of the new kmers to return. The functions 
//' returns the first num_kmers in the vector. Default value is 0, which returns 
//' all kmers.
//' @return A vector of doubles (actually integers represented as doubles), 
//' representing all the (new_len)-mers contained within the single (old_len)-mer.
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector unwrap_kmers_vect(NumericVector kmers, double old_len = 26, double new_len = 8, 
                           double base = 4, int num_kmers = 0){
    if (num_kmers == 0){
        num_kmers = old_len - new_len + 1;    
    }
    int n = kmers.size();
    NumericVector new_kmers(num_kmers * n);
    for (int j = 0; j < n; j++){
        double kmer = kmers[j];
        for (int i = 0; i < num_kmers; i++){
            double top = fmod(kmer,  pow(base, old_len - i));
            double bottom = pow(base, old_len - new_len - i);
            new_kmers[(j * num_kmers) + i] = floor(top / bottom);
        }
    }
    return new_kmers;
}