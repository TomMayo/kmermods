#include <Rcpp.h>
#include <math.h>
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

//' Converts a normal (base 10) integer to base representation (base 4 default) 
//'
//' Returns a base 4 (default) number as a vector calculated from the  base 10 
//' equivalent
//' 
//' @param number_b10 An integer in base 10, representing a kmer
//' @inheritParams base2kmer
//' @inheritParams kmer2base
//' @return A vector of integers representing a number in base 4 (default)
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
// [[Rcpp::export]]
NumericVector convert10to5_c(int number_b10, int k, int base = 4){
    if(number_b10 >= pow(base, k)){
        stop("k is not large enough to contain the number in this base");
    }
    NumericVector ret(k);
    if(k == 1){
        return number_b10;
    }
    for (int i = 1;i < k; i++){
        int entry = 0;
        int upper = (entry + 1) * pow(base ,(k - i));
        while (number_b10 >= upper){
            entry = entry + 1;
            upper = (entry + 1) * pow(base, (k - i));
        }
        ret[i-1] = entry;
        number_b10 = number_b10 - entry * pow(base, (k - i));
    }
    ret[k-1] = number_b10;
    return ret;
}

//' Converts a normal (base 10) integer to base representation (base 4 default) 
//'
//' Returns a base 4 (default) number as a vector calculated from the  base 10 
//' equivalent. Same as convert10to5_c but slightly faster
//' 
//' @param number_b10 An integer in base 10, representing a kmer
//' @inheritParams base2kmer
//' @inheritParams kmer2base
//' @return A vector of integers representing a number in base 4 (default)
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
NumericVector convert10tobase_c(double number_b10, int k, int base = 4){
    if(number_b10 >= pow(base, k)){
        stop("k is not large enough to contain the number in this base");
    }
    NumericVector ret(k);
    if(k == 1){
        return number_b10;
    }
    for (int i = 1 ;i < k; i++){
        double entry = floor(number_b10 / pow(base, (k - i)));
        ret[i-1] = entry;
        number_b10 = number_b10 - entry * pow(base, (k - i));
    }
    ret[k-1] = number_b10;
    return ret;
}

//' Counts the number of mismatches between two kmers, when represented in 
//' integer format 
//'
//' Returns a count of the number of mismatches between two kmers (of equal 
//' length), which are represented in integer format. For instance the two kmers
//' AAA and ATA have one mismatch. This function bypasses the conversion to 
//' strings.
//' 
//' @param kmer_1 An integer in base 10, representing a kmer
//' @param kmer_2 An integer in base 10, representing a kmer
//' @inheritParams base2kmer
//' @inheritParams kmer2base
//' @return A count of the mismatches
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
//' @export
// [[Rcpp::export]]
int mismatch_kmers(double kmer_1, double kmer_2, int k, int base = 4){
    NumericVector vec_1 = convert10tobase_c(kmer_1, k, base);
    NumericVector vec_2 = convert10tobase_c(kmer_2, k, base);
    int n = vec_1.size();
    double ret = 0;
    for (int i = 0; i < n; i++){
        if(vec_1[i] != vec_2[i]){ret += 1;}
    }
    return ret;
}
