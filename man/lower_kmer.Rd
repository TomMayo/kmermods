% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operations_post.R
\name{lower_kmer}
\alias{lower_kmer}
\title{Finds the (k-1)mer integer representation ffrom the kmer integer}
\usage{
lower_kmer(number_b10, k, base = 5)
}
\arguments{
\item{number_b10}{An integer in base 10, representing a kmer}

\item{k}{The length of kmers we are assessing. Used as a check.}

\item{base}{The base of our representation, default is 4. This corresponds to
the number of distinct nucleotides we are processing.}
}
\value{
An integer, representing the (reverse complementary kmer)k-1)mer
}
\description{
Taking in an integer representing the kmer at that site, it returns the 
integer for the (k-1)mer starting at the same place
}
\examples{
lower_kmer(100, 3, base = 4)
}
\author{
Tom Mayo \email{t.mayo@ed.ac.uk}
}

