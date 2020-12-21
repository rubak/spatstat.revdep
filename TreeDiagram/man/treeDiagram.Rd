\name{treeDiagram}

\alias{treeDiagram}

\title{Draw a tree diagram for herarchical tree-based method}

\usage{treeDiagram(data,treedat,cat_var,filename,pic_height=10,pic_width=10)}

\arguments{
  \item{data}{An non-empty data frame.Variable names should not contain speical characters such as ">", "<", or "="}
  \item{treedat}{A character string (in extended Newick's format as described below) or the first object returned by \code{tree()} from \pkg{tree} package.}
  \item{cat_var}{A character string indicaing the predictor(classification) variable name in tree method.}
  \item{filename}{A character string. The name of output plot}
  \item{pic_height}{A real positive numeric value. This argument is optional, and the default argument is set as 10}
  \item{pic_width}{A real positive numeric value. This argument is optional, and the default argument is set as 10}
}

\author{Wen Tian (Wendy) Wang (\email{wtwang@sfu.ca}), Lloyd T. Elliott}

\description{\code{treeDiagram} generates a tree diagram of any tree-based method and save automatically into user's current working directory. The classifying tree acts on the points provided in 'data' variable. The tree is specified by 'treedat' (for example, this could be a decision tree, or a character string generated from extended Newick's format, as described below). Both these variables must be provided. Output tree diagram will also provide tree information (e.g. node number, variable and value used at a split, minimum and maximum values of the split variable) for the first three nodes.
}

\details{\code{treeDiagram} visualizes hierarchical clustering by projecting each tree split into an one dimensional density plot of the partitioned data. These projections are then arranged through rotation and translation to indicate the topology of the tree. As these projections are rotated and translated, a single plot with a fractal-like organisation is formed. The tree diagram allows users to access the quality of cut in terms of linearly separability, depth, balance of tree and distribution and classification of partitioned data in each cut.
}

\note{
  \itemize{
    \item Predictor variable must be a binary factor variable.
    \item The extended Newick's format is revised such that it follows the following rules:
      \enumerate{
          \item If a parent node has two children nodes, use the form of "variable = value"
          \item If a parent node has only one child node, use the form of "var < value" if child node is on the left; otherwise use "var > value" for child node on the right

          \item parent node is always at the end of brackets.
          \item tip node is written in the form of "var = value"
          \item tree is read from left to right
      }
  \item Six required packages, "ggplot2","ape","cowplot","tree","stringr","spatstat", will be automatically installed if user has not yet installed.
  }
}

\examples{
\dontrun{
# read breast cancer data from UCI database website
cancer <- read.csv(
  url("https://archive.ics.uci.edu/ml/machine-learning-databases/00451/dataR2.csv"))
# make sure predictor variable is factored
cancer$Classification <- factor(cancer$Classification)
summary(cancer) #optional step giving an overview of the dataset

# set working directory (this step is optional)
# setwd("~/user location")

# e.g.1 : draw tree diagram according to the first object returned by 'tree()'
# create decision tree
library(tree)
t_cancer <- tree(Classification ~ ., data=cancer)

# plot tree diagram and save to your working directory
treeDiagram(cancer,t_cancer[[1]],"Classification","tree diagram for tree()")

# e.g.2 : draw tree diagram giving a newick format file of a tree
# newick format string of a decision tree for breast cancer
breast_cancer <- paste0(
  "(((BMI=25.745,BMI=29.722)Resistin=13.248)Age>44.5,",
  "(((((Age=70)Adiponectin<9.3482)BMI<32.275)Glucose<",
  "111)Leptin>7.93315)Age>48.5)Glucose=91.5;")


# plot tree diagram and save to your working directory
treeDiagram(cancer,breast_cancer,"Classification","tree diagram for newicks format file")
}}

\keyword{tree}
