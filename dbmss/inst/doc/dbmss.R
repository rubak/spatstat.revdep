## ----global_options, include=FALSE--------------------------------------------
set.seed(2018)

## ----wmppp, warning=FALSE, message=FALSE--------------------------------------
library("dbmss")
# Draw the coordinates of 10 points
X <- runif(10)
Y <- runif(10)
# Draw the point types.
PointType   <- sample(c("A", "B"), 10, replace=TRUE)
# Plot the point pattern. Weights are set to 1 ant the window is adjusted
par(mar=c(1,1,1,1))
plot(wmppp(data.frame(X, Y, PointType)), which.marks=2, main="")

## ----paracou, fig.width=4-----------------------------------------------------
data(paracou16)
# Plot (second column of marks is Point Types) 
par(mar=c(1,1,1,1))
plot(paracou16, which.marks=2, leg.side="right", main="")

## ----m, fig.width=6, fig.height=5---------------------------------------------
plot(Mhat(paracou16, , "V. Americana", "Q. Rosea"), main="")

## ---- fig.width=6, fig.height=5-----------------------------------------------
plot(KdEnvelope(paracou16, , ReferenceType="Q. Rosea", Global=TRUE), main="")

