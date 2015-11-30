library(kmermods)
context("Test Operations Post Counting Functions")

test_that("rev_comp converts does the reverse complement properly", {
    alph <- build_alphabet()
    expect_error(rev_comp(0,alph,1))
    expect_equal(rev_comp(5,alph,2),17)
    expect_equal(rev_comp(15,alph,2),15) # CG to CG
    expect_equal(rev_comp(20,alph,2),20) #MG to MG
    expect_equal(rev_comp(5^8-1,alph,8),0) #Ms to Gs
})

test_that("lower_kmer calculates the (k-1)mer properly", {
    alph <- build_alphabet()
    expect_equal(lower_kmer(20,2),4) # MG to M
    expect_equal(lower_kmer(17,5),3) #GGGCT to GGGC
    expect_equal(lower_kmer(32,3),6) # AAT to AA
    expect_equal(lower_kmer(6,2),1) # AA to A
    
})

test_that("kmer_dot_prod calculates the dot product as expected (note indexing
          is from 1 in a vector!)", {
              expect_equal(kmer_dot_prod(c(0,rep(1,5)),c(1,0,0,4,5)),1)
              expect_equal(kmer_dot_prod(c(0,rep(3,5)),c(1,0,0,4,5)),21)
              expect_equal(kmer_dot_prod(c(0,rep(4,5)),c(1,0,0,4,5)),26)
              expect_equal(kmer_dot_prod(c(0,rep(2,5),rep(3,2)),c(1,0,0,4,5)),9)
              expect_equal(kmer_dot_prod(c(4,4,4,3,3,2,1,2,1,2,0,0,0,4),c(1,0,0,4,5),warp=rep(1,14)),31)    
          })