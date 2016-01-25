library(kmermods)
context("Test Counting Functions")


test_that("kmer2base converts kmers to lower base representation", {
    alph <- c('G','A','T','C')
    expect_equal(kmer2base('G',alph,1),0)
    expect_equal(kmer2base('GAAAT',alph,5),c(0,1,1,1,2))
    expect_equal(kmer2base('MCTACCCCCC',alph,10),c(4,3,2,1,3,3,3,3,3,3))
    expect_error(kmer2base('MCTACCCCCC',alph,1))
})

test_that("base2kmer converts  lower base representation to kmers", {
    alph <- c('G','A','T','C')
    expect_equal(base2kmer(0,alph,1),'G')
    expect_equal(base2kmer(c(0,1,1,1,2),alph,5),'GAAAT')
    expect_equal(base2kmer(c(4,3,2,1,3,3,3,3,3,3),alph,10),'MCTACCCCCC')
    expect_error(base2kmer(150,alph,2))
})


test_that("base5to10 converts  lower base representation to integers", {
    expect_equal(base5to10(4,1),4)
    expect_equal(base5to10(c(1,0),2),5)
    expect_equal(base5to10(c(2,3),2),13)
    expect_equal(base5to10(c(1,1,1),3),25+5+1)
    expect_error(base5to10(c(1,1,1),2))
    expect_error(base5to10(c(5,4,1,0),4))
})

test_that("convert10to5 converts  integers to lower base representation", {
    expect_equal(convert10to5(4,1),4)
    expect_equal(convert10to5(5,2),c(1,0))
    expect_equal(convert10to5(25,3),c(1,0,0))
    expect_equal(convert10to5(27,3),c(1,0,2))
    expect_equal(convert10to5(30,3),c(1,1,0))
    expect_equal(convert10to5(39,3),c(1,2,4))
    expect_error(convert10to5(30,2))
    expect_error(convert10to5(7,1))
})

test_that("convert10tokmer converts  integers back to kmers", {
    alph <- c('G','A','T','C')
    expect_equal(convert10tokmer(0,alph,1),'G')
    expect_equal(convert10tokmer(0,alph,3),'GGG')
    expect_equal(convert10tokmer(14,alph,3),'GTM')
    expect_equal(convert10tokmer(39,alph,5),'GGATM')
    expect_equal(convert10tokmer(39,alph,8),'GGGGGATM')
    expect_equal(convert10tokmer(5^8-1,alph,8),'MMMMMMMM')
    expect_equal(convert10tokmer(5^8-5,alph,8),'MMMMMMMG')    
})

test_that("kmerto10 converts kmers to integers", {
    alph <- build_alphabet()
    expect_equal(kmerto10('G',alph),0)
    expect_equal(kmerto10('GGG',alph),0)
    expect_equal(kmerto10('GTM',alph),14)
    expect_equal(kmerto10('GGATM',alph),39)
    expect_equal(kmerto10('GGGGGATM',alph),39)
    expect_equal(kmerto10('MMMMMMMM',alph),5^8-1)
    expect_equal(kmerto10('MMMMMMMG',alph),5^8-5)      
})

test_that("convert10tokmer and kmerto10 are inverse functions", {
    alph <- build_alphabet()
    expect_equal(convert10tokmer(kmerto10('MMMMMMMG',alph),alph,8),'MMMMMMMG')
    expect_equal(convert10tokmer(kmerto10('MMMMMMMM',alph),alph,8),'MMMMMMMM')
    expect_equal(convert10tokmer(kmerto10('MCTCATMG',alph),alph,8),'MCTCATMG')
    expect_equal(convert10tokmer(kmerto10('ATATATAT',alph),alph,8),'ATATATAT')
})

test_that("kmer_counter returns the correct vector of integers", {
    alph <- build_alphabet()
    expect_equal(kmer_counter('GGG',alph,3),0)
    expect_equal(kmer_counter('GGGGGGGGGG',alph,3),rep(0,8))
    expect_equal(kmer_counter('GGGGGGGGGG',alph,3),rep(0,8))    
    expect_equal(kmer_counter('AAATT',alph,3),c(31,32,37))
    expect_equal(kmer_counter('GAATT',alph,4),c(32,125+37)) 
    expect_equal(kmer_counter('NNN',alph,2),c(NA,NA))
    expect_equal(kmer_counter('GGGNGG',alph,3),c(0,rep(NA,3)))    
})

test_that("kmer_dot_prod calculates the dot product as expected (note indexing
          is from 1 in a vector!)", {
    expect_equal(kmer_dot_prod(c(0,rep(1,5)),c(1,0,0,4,5)),1)
    expect_equal(kmer_dot_prod(c(0,rep(3,5)),c(1,0,0,4,5)),21)
    expect_equal(kmer_dot_prod(c(0,rep(4,5)),c(1,0,0,4,5)),26)
    expect_equal(kmer_dot_prod(c(0,rep(2,5),rep(3,2)),c(1,0,0,4,5)),9)
    expect_equal(kmer_dot_prod(c(4,4,4,3,3,2,1,2,1,2,0,0,0,4),c(1,0,0,4,5),warp=rep(1,14)),31)    
})