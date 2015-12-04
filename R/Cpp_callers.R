#' Represents the kmers present in the input with their integer values, 
#' as a vector (with C++ implementations)
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
#' alph <- c('G','A','T','C')
#' kmer_counter_c('GATCATCAT', alph, 4)
#' @export
kmer_counter_c <- function(dna,alph,k){
    dna <- strsplit(dna, "")[[1]]
    num_b5 <- sapply(dna, function(i) let2base_c(i, alph))
    ret <- sapply(1:(length(num_b5) - k + 1), function(i){
        if(any(is.na(num_b5[i:(i + k - 1)]))){
            NA
        } else {
            base5to10_c(num_b5[i:(i + k - 1)], k)
        }
    })
    return(ret)
}

#' Represents the kmers present in the input with their integer values, 
#' as a vector (with C++ implementations)
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
#' alph <- c('G','A','T','C')
#' kmer_counter_wrapper(dna_input,100,alph,4)}
#' @export
kmer_counter_wrapper_c <- function(dna_input, chunk_size, alph, k){
    len <- Biostrings::width(dna_input)
    num_chunks <- floor(len / chunk_size)
    kmer_vector <- unlist(pbapply::pblapply(1:num_chunks, function(i){
        start <- (i - 1) * chunk_size + 1
        end <- min(len, i * chunk_size + k - 1)
        dna <- as.character(dna_input[[1]][start:end])
        kmer_counter_c(dna, alph, k)
    }))
    return(kmer_vector)
}

#' Represents the kmers present in the input with their integer values, 
#' as a vector, done with a parallelised version
#'
#' This runs over the input calculating the integer value for every kmer that 
#' begins at that location. It returns the resutls as a vector. It acts as a
#' wrapper to kmer_counter, and takes in DNAStringSet objects from the Biostrings
#' package (available from Bioconductor)
#' 
#' @param dna_input The DNA sequence, normall a whole chromosom as a 
#' DNAStringSet object (see Biostrings package)
#' @param num_cores The number of cores to use. The default setting is the 
#' maximum available
#' @inheritParams kmer2base
#' @inheritParams let2base
#' @return An integer vector corresponding to the kmers at each starting point
#' in the DNA
#' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
#' @examples
#' \donttest{
#' data(dna_input)
#' alph <- c("G","A","T","C")
#' kmer_counter_para(dna_input,4,alph,4)}
#' @export
kmer_counter_para_c <- function(dna_input, num_cores=NaN, alph, k){
    if (is.na(num_cores)){
        num_cores <- detectCores()
    }
    alph_list <- lapply(1:num_cores, function(i) alph)
    k_list <- lapply(1:num_cores, function(i) k)
    len <- Biostrings::width(dna_input)
    
    chunk_size <- floor(len / num_cores)
    dna_list <- lapply(1:num_cores, function(i){
        start <- (i - 1) * chunk_size + 1
        end <- min(len, i * chunk_size + k - 1)
        if (i==num_cores) {end <- len}
        as.character(dna_input[[1]][start:end])
    })
    kmer_vector <- mcmapply(function(x, y, z) kmer_counter_c(x,y,z), dna_list, alph_list, 
             k_list, mc.cores=num_cores, mc.preschedule=FALSE)
    return(kmer_vector)
}