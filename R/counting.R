#' Sets up the alphabet for indexing nucleotides
#'
#' Returns a dataframe with each row being a nucleotide and an integer
#' 
#' @param int A vector of the integers to be used, in order. Default
#' settings are recommended.
#' @param let A vector of the letters corresponding to those integers. Default
#' settings are recommended.
#' @return Returns  a dataframe with each row being a nucleotide and an integer
#' rows denote equivalence betwwen the nucleotide and the integer
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' alph <- kmermods::build_alphabet()
#' int <- alph$int[alph$let=="G"]
#' @export
build_alphabet <- function(int = c(0,1,2,3,4), let = c("G","A","T","C","M")){
    if(length(int) != length(let)){
        stop("Integers and Letters are not the same length")
    }
    alph <- data.frame(let, int)
    return(alph)
}

#' Finds the integer that represents the nucleotide
#'
#' Returns an integer for {G,A,T,G} and NA for 'N'
#' 
#' @param letter A single nucleotide, as a character.
#' @param alph The alphabet we are using. A dataframe created using 
#' build_alphabet()
#' @return An integer, or NA if the input is 'N'
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
let2base <- function(letter, alph){
    if(!any(c("N" == letter, alph$let == letter))){
        stop("Invalid letter, must be ACGTM or N")
    }
    if(letter == "N"){
        return(NA)
    }
    return(alph$int[alph$let == letter])
}

#' Converts a kmer to a number in base 4 (default)
#'
#' Returns a vector of integers, which correspond to a single number in base 4
#' (default) or whatever is defined in build_alphabet()
#' 
#' @param kmer A character vector of nucleotides to be converted the 
#' representation
#' @param k The length of kmers we are assessing. Used as a check.
#' @inheritParams let2base
#' @return A vector of integers, representing a number in base 4(default)
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
kmer2base <- function(kmer, alph, k){
    if(k != nchar(kmer)){
        stop("k does not equal length of kmer")
    }
    sapply(1:k, function(i) {
        alph$int[alph$let == substr(kmer,i, i)]
    })
}


#' Converts a base representation (base 4 default) to the kmer
#'
#' Returns the kmer denoted by the integer
#' 
#' @param number A vector of integers representing a number in base 4 (default)
#' which correspond to a kmer.
#' @param base The base of our representation, default is 4. This corresponds to
#' the number of distinct nucleotides we are processing.
#' @inheritParams let2base
#' @inheritParams kmer2base
#' @return A kmer, corresponding to the number in base 4(default)
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
base2kmer <- function(number, alph, k, base = 5){
    if(length(number) != k){
        stop("The number should be entered as a vector of length k")
    }
    letters <- sapply(1:k, function(i) {
        alph$let[alph$int == number[i]]
        })
    kmer <- paste(letters, collapse = "",sep = "")
    return(kmer)
}

#' Converts a base representation (base 4 default) to normal base 10
#'
#' Returns an integer in base 10, calculated from the base 4 (default) 
#' representation
#' 
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @return An integer, corresponding to the number in base 4(default)
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
base5to10 <- function(number, k, base = 5){
    if(!all(number < base)){
        stop("Some integers are greater than the base")
    }
    if(k != length(number)){
        stop("k does not equal length of kmer")
    }
    sum(sapply(1:k, function(i) number[i] * base ^ (k - i)))
}


#' Converts a normal (base 10) integer to base representation (base 4 default) 
#'
#' Returns a base 4 (default) number as a vector calculated from the  base 10 
#' equivalent
#' 
#' @param number_b10 An integer in base 10, representing a kmer
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @return A vector of integers representing a number in base 4 (default)
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
convert10to5 <- function(number_b10, k, base = 5){
    if(number_b10 >= base ^ k){
        stop("k is not large enough to contain the number in this base")
    }
    ret <- rep(0, k)
    if(k == 1){
        return(number_b10)
    }
    for (i in 1:(k - 1)){
        entry <- 0
        upper <- (entry + 1) * (base ^ (k - i))
        while (number_b10 >= upper){
            entry <- entry + 1
            upper <- (entry + 1) * (base ^ (k - i))
        }
        ret[i] <- entry
        number_b10 <- number_b10 - entry * (base ^ (k - i))
    }
    ret[k] <- number_b10
    return(ret)
}

#' Converts a normal (base 10) integer to the equivalent kmer
#'
#' Returns the kmer that corresponds to the integer input
#' 
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @inheritParams convert10to5
#' @return A kmer corresponding to the integer representation
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' alph <- kmermods::build_alphabet()
#' convert10tokmer(112, alph, 4, base = 4)
#' @export
convert10tokmer <- function(number_b10, alph, k, base = 5){
    base5 <- convert10to5(number_b10, k, base)
    kmer <- base2kmer(base5, alph, k, base)
    return(kmer)
}

#' Converts a kmer to it's normal (base 10) integer representation
#'
#' Returns the integer in base 10 that corresponds to the kmer
#' 
#' @inheritParams base2kmer
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @return An integer corresponding to the kmer
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' alph <- kmermods::build_alphabet()
#' kmerto10('GATC', alph, base = 4)
#' @export
kmerto10 <- function(kmer, alph, base = 5){
    k <- nchar(kmer)
    base5 <- kmer2base(kmer, alph, k)
    return(base5to10(base5, k, base))
}




#' Represents the kmers present in the input with their integer values, 
#' as a vector
#'
#' This runs over the input calculating the integer value for every kmer that 
#' begins at that location. It returns the resutls as a vector.
#' 
#' @param dna The DNA sequence, as a character
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @return An integer vector corresponding to the kmers at each starting point
#' in the DNA
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' alph <- kmermods::build_alphabet()
#' kmer_counter('GATCATCAT', alph, 4)
#' @export
kmer_counter <- function(dna,alph,k){
    dna <- strsplit(dna, "")[[1]]
    num_b5 <- sapply(dna, function(i) let2base(i, alph))
    ret <- sapply(1:(length(num_b5) - k + 1), function(i){
        if(any(is.na(num_b5[i:(i + k - 1)]))){
            NA
        } else {
            base5to10(num_b5[i:(i + k - 1)], k)
        }
    })
    return(ret)
}

#' Represents the kmers present in the input with their integer values, 
#' as a vector
#'
#' This runs over the input calculating the integer value for every kmer that 
#' begins at that location. It returns the resutls as a vector. It acts as a
#' wrapper to kmer_counter, and takes in DNAStringSet objects from the Biostrings
#' package (available from Bioconductor)
#' 
#' @param dna_input The DNA sequence, as a DNAStringSet object 
#' (see Biostrings package)
#' @param chunk_size The length of DNA to send in each chunk for parallel 
#' processing by kmer_counter
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @return An integer vector corresponding to the kmers at each starting point
#' in the DNA
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' \donttest{
#' data(dna_input)
#' alph <- kmermods::build_alphabet()
#' kmer_counter_wrapper(dna_input,100,alph,4)}
#' @export
kmer_counter_wrapper <- function(dna_input, chunk_size, alph, k){
    len <- Biostrings::width(dna_input)
    num_chunks <- ceiling(len / chunk_size)
    kmer_vector <- unlist(pbapply::pblapply(1:num_chunks, function(i){
        start <- (i - 1) * chunk_size + 1
        end <- min(len, i * chunk_size + k - 1)
        dna <- as.character(dna_input[[1]][start:end])
        kmer_counter(dna, alph, k)
    }))
    return(kmer_vector)
}