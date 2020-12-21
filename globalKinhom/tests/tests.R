require(globalKinhom)

set.seed(49)

lambda <- funxy(function(x,y) 200*(.5*exp(-((x-.5)^2 + (y-.5)^2)*2) + .5), owin())

X <- rpoispp(lambda)
Y <- rpoispp(lambda)
 
g.truef <- pcfglobal(X,lambda=lambda)
c.truef <- pcfcross.global(X, Y,lambdaX=lambda, lambdaY=lambda)
K.truef <- Kglobal(X,lambda=lambda)
Kcross.truef <- Kcross.global(X, Y,lambdaX=lambda, lambdaY=lambda)
print("done with estimates using true gamma")

g <- pcfglobal(X, interpolate=TRUE)
c <- pcfcross.global(X, Y, interpolate=TRUE)
K <- Kglobal(X, isotropic=TRUE, interpolate=TRUE)
Kcross <- Kcross.global(X, Y, isotropic=TRUE, interpolate=TRUE)
print("done with estimates using estimated gamma")

print(g.truef$global)
print(g$global)
print(c.truef$global)
print(c$global)
print(K.truef$global)
print(K$global)
print(Kcross.truef$global)
print(Kcross$global)
