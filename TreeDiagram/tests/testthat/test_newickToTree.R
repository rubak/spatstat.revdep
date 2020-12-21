library(testthat)
context("newickToTree")

# test 1) correctly turning a newick's string to a list of node number
test_that("correctly turning a newick's string to a list of node number", {
  # newick format string of a decision tree for breast cancer
  breast_cancer<-"(((BMI=25.745,BMI=29.722)Resistin=13.248)Age>44.5,(((((Age=70)Adiponectin<9.3482)BMI<32.275)Glucose<111)Leptin>7.93315)Age>48.5)Glucose=91.5;"
  result = newickToTree(breast_cancer)
  
  expect_equal(result$node_index, c(1,2,3,5,7,10,11,15,30,60,120)) 
})


# test 2) return error message if user typed tree string with wrong Newick's format
test_that("wrong Newick's format string", {
  # newick format string of a decision tree for breast cancer
  breast_cancer<-"((Resistin=13.248(BMI=25.745,BMI=29.722))Age>44.5,(((((Age=70)Adiponectin<9.3482)BMI<32.275)Glucose<111)Leptin>7.93315)Age>48.5)Glucose=91.5;"

  # correctly turning a newick's string to a list of node number
  expect_error(newickToTree(breast_cancer), "Error: Incorrect Newick's format is detected. Please check if your parent node is at the end of the brackets of child nodes") 
})



# test 3) return error if user forget to input variable name or value for split node
test_that("missing variable name(s) or value(s) for split node(s)", {
  
  # newick format string with missing values 
  # note that BMI has missing value: "(((BMI=,BMI=29.722)...."
  missing_value_breast_cancer<-"(((BMI=,BMI=29.722)Resistin=13.248)Age>44.5,(((((Age=70)Adiponectin<9.3482)BMI<32.275)Glucose<111)Leptin>7.93315)Age>48.5)Glucose=91.5;"
  # expecting an warming message pointing which row contains missing value(s)
  expect_warning(newickToTree(missing_value_breast_cancer))
  
  # newick format string with missing variable name
  # note that string below has missing variable names"(((BMI=25.745,=29.722)=13.248)..." 
  missing_variable_breast_cancer<-"(((BMI=25.745,=29.722)=13.248)Age>44.5,(((((Age=70)Adiponectin<9.3482)BMI<32.275)Glucose<111)Leptin>7.93315)Age>48.5)Glucose=91.5;"
  # expecting an warming message pointing which row contains missing variable names(s)
  expect_warning(newickToTree(missing_variable_breast_cancer))
  
})