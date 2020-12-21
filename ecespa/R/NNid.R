# changed in ecespa 1.1.10. to remove dependency of splancs
# using nnwhich() instead of nndistG()$neighs
# It was
# `NNid` <-function (xy, splancs = TRUE)

`NNid` <-
function (xy) 
{
    nnwhich(xy)
    # removred in ecespa 1.1.10 from here
 #   if (splancs) {
 #       nndistG(xy)$neighs
 #   }
 #    to here

}

