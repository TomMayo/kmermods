#' Finds the reverse complementary kmer, in integer form, without methylation
#'
#' Taking in an integer in base 10, it gives the integer in base 10 that 
#' represents the reverse complement. NOT FULLY TESTED
#' 
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @inheritParams convert10to5
#' @return An integer, representing the reverse complementary kmer
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' alph <- kmermods::build_alphabet()
#' rev_comp_meth(100,alph,3)
#' @export
rev_comp <- function(number_b10, alph, k){
    if(k == 1){
        stop("k = 1, this function is not for trivial case of 1-mers")
    }
    
    comp <- 4^k - 1 - number_b10
    base_num <- convert10to5(comp, k, base = 4)
    
    revcomp <- rev(base_num)    # reverse
    return(base5to10(revcomp, k, base = 4))    # convert back to integer
}


#' Finds the reverse complementary kmer, in integer form, using methylation
#'
#' Taking in an integer in base 10, it gives the integer in base 10 that 
#' represents the reverse complement. NOT FULLY TESTED
#' 
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @inheritParams convert10to5
#' @return An integer, representing the reverse complementary kmer
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' alph <- kmermods::build_alphabet()
#' rev_comp_meth(100,alph,3)
#' @export
rev_comp_meth <- function(number_b10, alph, k){
    if(k == 1){
        stop("k = 1, this funvtion is not for trivial case of 1-mers")
    }
    base5 <- convert10to5(number_b10, k)
    # store M locations
    Ms <- which(base5 == kmer2base("M", alph, 1))
    # change all Ms to Cs
    base5[Ms] <- kmer2base("C", alph, 1)
    # comp
    comp <- rep(3, k) - base5
    # add back Ms in CpG methylation
    for (i in Ms) {
        if(i == 1){
            if(comp[i + 1] == kmer2base("C", alph, 1)){
                comp[i + 1] <- kmer2base("M", alph, 1)
            }
        } else if (i == k){
            if(comp[k - 1] == kmer2base("C",alph,1)){
                comp[k - 1] <- kmer2base("M",alph,1)
            }
        } else {
            if(comp[i - 1] == kmer2base("C",alph,1)){
                comp[i - 1] <- kmer2base("M",alph,1)
            }
            if(comp[i + 1] == kmer2base("C",alph,1)){
                comp[i + 1] <- kmer2base("M",alph,1)
            }
        }}
    revcomp <- rev(comp)    # reverse
    return(base5to10(revcomp,k))    # convert back to integer
}

#' Finds the (k-1)mer integer representation ffrom the kmer integer
#'
#' Taking in an integer representing the kmer at that site, it returns the 
#' integer for the (k-1)mer starting at the same place
#' 
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @inheritParams convert10to5
#' @return An integer, representing the (reverse complementary kmer)k-1)mer
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' lower_kmer(100, 3, base = 4)
#' @export
lower_kmer <- function(number_b10, k, base = 5){
    # put into base5, delete the last entry, put it back
    # same as round down to nearest 5 and then divide by 5
    return(floor(number_b10/base))
}

#' Computes the dot product of the parameters with the kmer counts, with or 
#' without it being warped
#'
#' Taking in an integer representing the kmer at that site, it returns the 
#' integer for the (k-1)mer starting at the same place
#' 
#' @param kmers is a vector of integers of any length representing kmers over
#'  a window
#' @param paras is a vector of length equal to the total number of kmers
#' @param warp is a vector of length as long as the kmer vector, with the 
#' multiplicative weights for how much to warp the entry 
#' @return An double, representing the dot product of the parameters with what 
#' would usually be the vector of kmer abundances (warped or not)
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples 
#' kmers <- c(0,7,89,45,75,65,22,12)
#' paras <- rep(1,100)
#' kmer_dot_prod(kmers,paras)
#' @export
kmer_dot_prod <- function(kmers, paras, warp = NA){
    if(is.na(warp[1])){
        return(sum(sapply(kmers, function(ind) paras[ind + 1])))
    } else {
        return(sum(sapply(1:length(kmers), function(i){
            ind <- kmers[i]
            paras[ind + 1] * warp[i]
            })))
    }
}
