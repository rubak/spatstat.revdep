library(testthat)
context("get_nodes_list")

test_that("length of returning object from get_nodes_list() should be the 
          same as number of split nodes in a tree if user newick's string 
          or object return from tree_info() is true", {
  # newick format string of a decision tree for breast cancer
  breast_cancer<-"(((BMI=25.745,BMI=29.722)Resistin=13.248)Age>44.5,(((((Age=70)Adiponectin<9.3482)BMI<32.275)Glucose<111)Leptin>7.93315)Age>48.5)Glucose=91.5;"
  result = newickToTree(breast_cancer)

  # correctly returning a nested list of nodes number along each path from root to such node.
  # here we want to check if the length of this nested list matches the total number of nodes of a tree for breast cancer data
  expect_equal(length(get_nodes_list(result$node_index)),11)

})

test_that("testing for empty node list", {
  # newick format string of a decision tree for breast cancer
  emptyNodeList <- NULL
  # correctly returning a nested list of nodes number along each path from root to such node.
  # here we want to check if the length of this nested list matches the total number of nodes of a tree for breast cancer data
  expect_error(get_nodes_list(emptyNodeList),"there's no splitting node")
  
})