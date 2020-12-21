#'
#'   undocumented functions in spatstat.gui
#'

require(spatstat.gui)
local({

  ## infrastructure of iplot, istat
  fakepanel <- list(x=redwood,
                    xname="redwood",
                    envel="none",
                    stat="data",
                    sigma=0.08,
                    pcfbw=0.01,
                    simx=rpoispp(ex=redwood, nsim=39))

  a <- do.istat(fakepanel)
  for(envel in c("none", "pointwise", "simultaneous")) {
    fakepanel$envel <- envel
    for(stat in c("density",
                  "Kest", "Lest", "pcf",
                  "Kinhom", "Linhom", "Fest", "Gest", "Jest")) {
      fakepanel$stat <- stat
      a <- do.istat(fakepanel)
    }
  }

})


