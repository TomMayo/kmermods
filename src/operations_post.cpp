#include <Rcpp.h>
#include <algorithm>
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

//' Computes the vector with which to update the parameters in a logistic regression
//' onto the peaks
//'
//' This calculates the prediction for each region, in a non-sliding scheme, 
//' calculates the error and returns the update vector for all of the regions.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' //' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A vector, representing the amount to update the parameter vector
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector params_peaks_noslide(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params, 
                                   NumericMatrix peaks, int win_size, int chrom_loc,
                                   nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    NumericVector update;
    int num_params = params.size();
    update = rep(0.0, num_params);
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    int ind;
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            err_term = peak - pred;
            for(int j = 0; j < win_size; j++) {
                ind = kmers[j];
                update[ind] += err_term;
            }
            update[num_params - 1] += err_term;
        }
    }
    return update;
}


//' Computes the total error of the predictions.
//'
//' This calculates the total of the absolute values of the errors, given a set
//' of parameters and outputs, over a given region.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return The total of the absolute errors, with a breakdown
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector total_error(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params, 
                                           NumericMatrix peaks, int win_size, int chrom_loc,
                                           nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    int num_params = params.size();
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    int ind;
    NumericVector err_sum;
    err_sum = rep(0.0, 5);    
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            err_term = peak - pred;
            err_sum[0] += fabs(err_term);
            if(pred >= 0.5){
                if (peak==1.0){
                    err_sum[1] += 1.0; // true pos
                } else {
                    err_sum[2] += 1.0; // false pos
                }
            } else {                
                if (peak==1.0){
                    err_sum[3] += 1.0; // false neg 
                } else {
                    err_sum[4] += 1.0; // true neg 
                }
            }
        }
    }
    return err_sum;
}

//' Computes the vector with which to update the parameters in a logistic regression
//' onto the peaks and also returns error calculations
//'
//' This calculates the prediction for each region, in a non-sliding scheme, 
//' calculates the error and returns the update vector for all of the regions.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' //' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A list, the first element is a vector representing the amount to 
//' update the parameter vector, the second is the sum of errors of different 
//' types
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
Rcpp::List params_peaks_noslide_w_error(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params, 
                                   NumericMatrix peaks, int win_size, int chrom_loc,
                                   nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    NumericVector update;
    int num_params = params.size();
    update = rep(0.0, num_params);
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    int ind;
    NumericVector err_sum;
    err_sum = rep(0.0, 5);
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            err_term = peak - pred;
            for(int j = 0; j < win_size; j++) {
                ind = kmers[j];
                update[ind] += err_term;
            }
            update[num_params - 1] += err_term;
            err_term = peak - pred;
            // now calculate error
            err_sum[0] += fabs(err_term);
            if(pred >= 0.5){
                if (peak==1.0){
                    err_sum[1] += 1.0; // true pos
                } else {
                    err_sum[2] += 1.0; // false pos
                }
            } else {                
                if (peak==1.0){
                    err_sum[3] += 1.0; // false neg 
                } else {
                    err_sum[4] += 1.0; // true neg 
                }
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("update") = update,
                              Rcpp::Named("error") = err_sum);
}


//' Computes the predictions for the logistic regression
//'
//' This calculates the predictions at each window, returning a probability of 
//' seeing a peak there
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A vector of probabilities
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericMatrix predict_peaks(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params,
                            NumericMatrix peaks, int win_size, int chrom_loc,
                          nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    if (num_res < 0){
        NumericMatrix ret(2,0);
        return ret;
    }
    int num_params = params.size();
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    int num_cols = reg_len / win_size;
    NumericMatrix ret(2, num_cols);
    int current_reg = 0;
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            ret(0, current_reg) = pred;
            ret(1, current_reg) = peak;
            current_reg += 1;
        }
    }
    return ret;
}


//' L1 regulatisation proximal operator
//'
//' This function computes the proximal operator for L1- regularised regression
//' (lasso) and returns the new vector.
//' 
//' @param params is a vector of length equal to the total number of kmers, 
//' representing the parameters in the model
//' @param thresh is the threshold for the proximal operator for l1 regularised
//' regression
//' @return A vector of parameters
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector l1_prox_op(Rcpp::NumericVector params, double thresh) {
    int num_params = params.size();
    NumericVector new_params = rep(0.0,num_params);
    for (int i = 0; i < num_params; i ++){
        if (fabs(params[i]) > thresh){
            new_params[i] = params[i] - thresh;
            if (params[i] < 0){
                new_params[i] += (2 * thresh);
            }
        }
    }
    return(new_params);
}


//' Calculates the log likelihood of the data, given some parameters for a 
//' logistic regression
//'
//' Given the kmers, set of peaks and parameters, this returns the log likelihood.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' //' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A double, representing the log likelihood
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
double loglik_logreg(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params, 
                     NumericMatrix peaks, int win_size, int chrom_loc,
                     nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    int num_params = params.size();
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    double loglik = 0;
    double pseudo_prob = 0.00001;
    int ind;
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            pred = std::min(pred, 1.0 - pseudo_prob);
            pred = std::max(pseudo_prob, pred);
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            loglik += peak * log(pred) + (1.0 - peak) * log(1 - pred);
        }
    }
    return loglik;
}

//' Computes the gradient of the log likelihood of a logistic regression
//'
//' This calculates the prediction for each region, in a non-sliding scheme, 
//' calculates the error and returns the update vector for all of the regions.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' //' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A vector, representing the gradient (direction to update the 
//' parameter vector)
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector loglik_log_reg_grad(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params, 
                                  NumericMatrix peaks, int win_size, int chrom_loc,
                                  nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    int num_params = params.size();
    NumericVector update;
    update = rep(0.0, num_params);
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    int ind;
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            err_term = peak - pred;
            for(int j = 0; j < win_size; j++) {
                ind = kmers[j];
                update[ind] += err_term;
            }
            update[num_params - 1] += err_term;
        }
    }
    return update;
}

//' Computes the gradient of the log likelihood of a logistic regression
//'
//' This calculates the prediction for each region, in a non-sliding scheme, 
//' calculates the error and returns the update vector for all of the regions.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' //' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A vector, representing the gradient (direction to update the 
//' parameter vector)
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
Rcpp::List log_reg_wrapper(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params, 
                           NumericMatrix peaks, int win_size, int chrom_loc,
                           nullable_t warp_ = R_NilValue) {
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    NumericVector update;
    int num_params = params.size();
    update = rep(0.0, num_params);
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    double loglik = 0;
    int ind;
    NumericVector err_sum;
    err_sum = rep(0.0, 5);    
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        double lin_prod;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_prod = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            lin_prod =Rcpp::sum(tmp) + params[num_params - 1];    
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            lin_prod = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_prod); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            pred = 1.0 / (1.0 + exp(-lin_prod));
            // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
            // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
            // peak_stop << ", peak count " << peak_count;
            err_term = peak - pred;
            for(int j = 0; j < win_size; j++) {
                ind = kmers[j];
                update[ind] += err_term;
            }
            update[num_params - 1] += err_term;
            loglik += peak * log(pred) + (1.0 - peak) * log(1 - pred);
            err_sum[0] += fabs(err_term);
            if(pred >= 0.5){
                if (peak==1.0){
                    err_sum[1] += 1.0; // true pos
                } else {
                    err_sum[2] += 1.0; // false pos
                }
            } else {                
                if (peak==1.0){
                    err_sum[3] += 1.0; // false neg 
                } else {
                    err_sum[4] += 1.0; // true neg 
                }
            }
        }
    }
    return Rcpp::List::create(Rcpp::Named("gradient") = update,
                              Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("error") = err_sum);
}

//' Calculates the log likelihood of the data, given some parameters for a 
//' logistic regression
//'
//' Given the kmers, set of peaks and parameters, this returns the log likelihood.
//' 
//' @param kmers_win is a vector of integers of any length representing kmers in
//' a region
//' @param paras is a vector of length equal to the total number of kmers
//' @param grad is a vector of length equal to the total number of kmers, 
//' representing the gradient of the log likelihood at that point
//' @param alphas is a vector of alphas, giving us the values for the line 
//' search, this will test loglik at params + alpha*grad for each alpha
//' @param peaks is a matrix giving the locations of the peaks on the chromosome,
//' the first column is starts, second is ends, inclusive, indexed from 1
//' @param win_size is the length of the sliding window we are using
//' @param chrom_loc is the position of the first kmer along the chromosome - 
//' this avoids indexing errors when splitting up the data
//' //' @param warp is a vector of length as long as the kmer vector, with the 
//' multiplicative weights for how much to warp the entry
//' @return A double, representing the log likelihood
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector loglik_linesearch(Rcpp::IntegerVector kmers_win, Rcpp::NumericVector params,
                         Rcpp::NumericVector grad, Rcpp::NumericVector alphas,
                         NumericMatrix peaks, int win_size, int chrom_loc,
                         nullable_t warp_ = R_NilValue) {
    int num_alphas = alphas.size();
    NumericVector logliks = rep(0.0, num_alphas);
    int reg_len = kmers_win.size();
    int num_res = reg_len - win_size + 1; // number of outcomes
    int num_params = params.size();
    int half_win = floor(win_size / 2);
    int peak_count = 0;
    int num_peaks = peaks.nrow();
    int peak_start = peaks(0,0);
    int peak_stop = peaks(0,1);
    bool more_peaks = true;
    double pseudo_prob = 0.00001;
    int ind;
    double lin_prod;
    double lin_sum;
    double grad_sum;
    for (int i = 0; i < num_res; i = i + win_size){
        // define the kmers for the window
        IntegerVector kmers(win_size);
        for(int j = 0; j < win_size; j++){
            kmers[j] = kmers_win[i+j];
        };
        // run the dot product
        NumericVector lin_prods;
        if (Rcpp::na_omit(kmers).size() != kmers.size()) {
            lin_sum = NA_REAL;
        } else if (warp_.isNull()) {
            Rcpp::NumericVector tmp = params[kmers];
            Rcpp::NumericVector tmp_grad = grad[kmers];
            lin_sum = Rcpp::sum(tmp) + params[num_params - 1];    
            grad_sum = Rcpp::sum(tmp_grad) + grad[num_params - 1];
        } else {
            Rcpp::NumericVector warp(warp_), tmp = params[kmers];
            Rcpp::NumericVector tmp_grad = grad[kmers];
            lin_sum = Rcpp::sum(tmp * warp) + params[num_params - 1]; 
            grad_sum = Rcpp::sum(tmp_grad * warp) + grad[num_params - 1];
        }
        //apply sigmoid
        double pred;
        double err_term;
        bool test = NumericVector::is_na(lin_sum); 
        if(!test){
            int peak_loc = i + half_win + chrom_loc; 
            double peak = 0.0;
            if (more_peaks){
                if(peak_loc < peak_start){
                    peak = 0.0;
                } else if (peak_loc < peak_stop){
                    peak = 1.0;
                } else {
                    while ((peak_count < (num_peaks-1)) && (peak_loc > peak_stop)){
                        //  Rcout << "\nGoing up the peak list at i = " << i;
                        peak_count += 1;
                        peak_start = peaks(peak_count, 0);
                        peak_stop = peaks(peak_count, 1);
                    }
                    if(peak_loc < peak_start){
                        peak = 0.0;
                    } else if (peak_loc < peak_stop){
                        peak = 1.0;
                    } else if (peak_count==(num_peaks-1) && peak_loc > peak_stop){
                        more_peaks = false;
                    }
                }
            }
            for (int l = 0; l < num_alphas; l++){   
                lin_prod = lin_sum + (alphas[l] * grad_sum);
                pred = 1.0 / (1.0 + exp(-lin_prod));
                pred = std::min(pred, 1.0 - pseudo_prob);
                pred = std::max(pseudo_prob, pred);
                // Rcout << "\npred " << pred << ", peak " << peak << ", peak_loc "<<
                // peak_loc << ", peak_start" << peak_start << ", peak_stop " <<
                // peak_stop << ", peak count " << peak_count;
                logliks[l] += peak * log(pred) + (1.0 - peak) * log(1.0 - pred);
            }
        }
    }
    return logliks;
}


