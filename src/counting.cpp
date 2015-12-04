#include <Rcpp.h>
using namespace Rcpp;

//' Finds the integer that represents the nucleotide
//'
//' Returns an integer for {G,A,T,G} and NA for 'N'
//' 
//' @param letter A single nucleotide, as a character.
//' @param alph_vect The alphabet we are using. A dataframe created using 
//' build_alphabet()
//' @return An integer, or NA if the input is 'N'
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
// [[Rcpp::export]]
int let2base_c(String letter, CharacterVector alph_vect){
    CharacterVector temp = alph_vect;
    temp.push_back("N");
    int len = temp.size();
    int i = 0;
    bool test = true;
    while (test){
        if (temp[i] == letter){
            test = false;
        }
        if(test){
            i++;
        }
        if(i == len){
            stop("Invalid letter, must be ACGTM or N");
        }
    }
    if(i == len - 1){
        double ret = NA_REAL;
        return ret;
    }
    return i;
}

//' Converts a base representation (base 4 default) to normal base 10
//'
//' Returns an integer in base 10, calculated from the base 4 (default) 
//' representation
//' 
//' @inheritParams base2kmer
//' @inheritParams kmer2base
//' @return An integer, corresponding to the number in base 4(default)
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
// [[Rcpp::export]]
int base5to10_c(NumericVector number,int k,int base = 5){
    if(!all(number < base).is_true()){
        stop("Some integers are greater than the base");
    }
    if(k != number.size()){
        stop("k does not equal length of kmer");
    }
    int counter = 0;
    for (int i = 0;i < k; i++){
        counter += number[i] * pow(base, k - i - 1);
    }
    return counter;
}

// // [[Rcpp::export]]
// int kmer_counter_c(std::string dna, NumericVector alph_vect, int k){
//     int base = alph_vect.size();
//     int len = dna.size();
//     NumericVector base_vect(len);
//     for (int i = 0; i < len; i++){
//         int base_vect[i] = let2base_c(dna[i], alph_vect);
//     }
//     NumericVector ret(len - k);
//     for (int i = 0; i < len - k; i++){
//         NumericVector number = base_vect[i:i + k - 1];
//         if(any(is.na(number[i:(i + k - 1)]))){
//             double ret[i] = NA_REAL;
//         } else {
//             int ret[i] = base5to10_c(number, k, base);
//         }
//     }
//     return ret;
// }