context("test-phylo2L")

test_that("phylogeny to L table and back results in the original phylogeny", {

  phylogeny <- ape::read.tree(text = "((t1:1, t3:1):1, t2:2);")
  l_table <- phylo2L(phylogeny)
  phylogeny_again <- DDD::L2phylo(l_table)
  expect_equal(phylogeny, phylogeny_again)
})
