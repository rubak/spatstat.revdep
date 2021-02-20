Deployment status for the new spatstat 2.0 structure
================
Ege Rubak and Adrian Baddeley

# Repository contents

This repository contains the result of running a script to update
reverse dependencies of spatstat to the new spatstat 2.0 structure (see
<https://github.com/rubak/spatstat.revdep.helpers> for details).

# Current status and recommendations

The rest of this README is intended to give a continuous update of the
current recommendations and state of affairs for spatstat 2.0, with
planned submission to CRAN on **Monday 22 Feb. 2021**.

To adapt your package to spatstat please follow the advice sent by email
to authors of reverse dependencies and further detailed below. Once you
have followed the steps below to update your package to spatstat 2.0 you
may experience problems checking your package if it depends on other
packages that have not yet been adapted to statstat 2.0. In particular
you may need updated versions of `sf` and `maptools` as explained in
[this issue](https://github.com/rubak/spatstat.revdep/issues/2), where
you can also help us by announcing a developement version of your own
package for others to install. If you need an updated dependency you can
search for it in the issue, and if it is missing you may try to reach
out to Ege via the issue or by direct email.

## Current status of packages

We have registred the following packages as updated and ready for CRAN
submission once spatstat 2.0 is accepted on CRAN:

``` r
csvfile <- file.path(tempdir(), "spatstat.csv")
download.file("https://docs.google.com/spreadsheets/d/1RTmP2iwZGlxSUToYaht-H1SlS4CyFgbQ0Umv5i3x-1Q/export?format=csv", 
              destfile = csvfile)
done <- trimws(unlist(strsplit(read.csv(csvfile)[,2], ",")))
unlink(csvfile)
sort(done)
```

    ##  [1] "AROC"           "CalSim"         "dixon"          "ecespa"        
    ##  [5] "envi"           "ETAS"           "gateR"          "GmAMisc"       
    ##  [9] "graph4lg"       "idar"           "inlabru"        "maptools"      
    ## [13] "overlapptest"   "rcarbon"        "replicatedpp2w" "RImageJROI"    
    ## [17] "ROCnReg"        "selectspm"      "sf"             "shar"          
    ## [21] "spagmix"        "sparr"          "sparrpowR"      "stars"         
    ## [25] "tabularaster"   "trajectories"   "trip"

If your package should be included on this list of packages please
submit its name using [this
webform](https://forms.gle/W2VzUshfM3zWaQow9).

## Recommended procedure for updating

1.  Install all the sub-packages from CRAN by just typing

    ``` r
    install.packages('spatstat.linnet', dependencies=TRUE)
    ```

2.  Install spatstat 2.0-0 from <https://spatstat.r-universe.dev>

    ``` r
    install.packages("spatstat", repos = "https://spatstat.r-universe.dev")
    ```

3.  Modify your package so that it works with spatstat 2.0-0 and the
    relevant spatstat subpackages (please see details further below).
    You can find scripts to help you with this at
    <https://github.com/rubak/spatstat.revdep.helpers> and a result of
    running these scripts for your package at the repository you are
    currently in: <https://github.com/rubak/spatstat.revdep>

4.  Run the CRAN package checker on your modified package (possibly with
    development versions of other packages installed as detailed above)

5.  Let us know that your package is updated and ready for submission by
    typing its name in [the webform mention
    previously](https://forms.gle/W2VzUshfM3zWaQow9). This is to help us
    tell CRAN that maintainers are ready with updated packages once
    spatstat 2.0-0 is on CRAN and makes old packages break.

We strongly recommend that you

**DO NOT submit your updated package to CRAN yet**

**DO help us by adding `spatstat (>= 2.0-0)` to
`Imports`/`Depends`/`Suggests` as detailed below even if it isn’t
strictly necessary.**

(Explanation: many other packages depend on `spatstat`. There are still
many copies of the old spatstat lying around. When your package depends
on `spatstat.geom` (for example) then the old and new code may coexist
on some systems. and this will cause problems, so we really want to get
spatstat 2.0-0 to any system that installs spatstat.xxxx. See ‘Collision
Avoidance’ below)

How to use spatstat 2.0-0 in your package depends on how you currently
use spatstat:

#### Depends:

If your package currently ‘Depends’ on spatstat, the general advice is
to have the following in `DESCRIPTION` (in addition to other packages
you use):

    Depends: spatstat (>= 2.0-0)
    Imports: spatstat.geom, spatstat.core, spatstat.linnet

It may seem excessive to have all the sub-packages here, but they are
strictly needed if you currently use the double-colon idiom `spatstat::`
in your code. If not you could omit those that you don’t use.

In `NAMESPACE` you need something like:

    import(spatstat.geom, spatstat.core, spatstat.linnet, spatstat)

#### Suggests:

If your package currently ‘Suggests’ spatstat, the general advice which
should work out of the box for most cases (together with the automated
code changes by the script):

    Suggests: spatstat.geom, spatstat.core, spatstat.linnet, spatstat (>= 2.0-0)

You can omit any `spatstat.xxxx` that you don’t use, and in the long run
you can omit `spatstat (>= 2.0-0)`, but we would really appreciate it if
you can keep it in for now.

See ‘Collision Avoidance’ below, for some code which you can use to
automatically check for possible conflicts.

#### Imports:

For DESCRIPTION:

    Imports: spatstat.geom, spatstat.core, spatstat.linnet, spatstat (>= 2.0-0)

Where you *must* delete any sub-packages you don’t actually use. For
`NAMESPACE` you have to add mandatory `imports(spatstat.xxxx)` or
preferably `importFrom(spatstat.xxxx, yyyy)` for all the
subpackages/functions you use, and if you kindly put
`spatstat (>= 2.0-0)` in `DESCRIPTION` then you should add a line in
`NAMESPACE`:

    import("spatstat")

This line can be deleted in a few months together with the dependence on
umbrella spatstat when the dust has settled, and then you should only
depend on relevant sub-packages.

Much of this code update may be achieved by using the automated script,
but many will also need manual changes to `DESCRIPTION` and `NAMESPACE`
as detailed above. Please reach out if you are in doubt or need help. In
particular if your package is on GitHub Ege is happy to run the script
for you and make a pull request with the relevant changes if you send
him an email with the repository information.

## Collision Avoidance

Maintainers of packages with `spatstat` in `Suggests` will typically
change their package so that it `Suggests` the relevant sub-packages
`spatstat.xxxx` as well. Any code which requires `spatstat` would
typically have been wrapped in

``` r
if(require(spatstat)) {
    # Run spatstat code
}
```

This should now change to `require(spatstat.geom)`, etc. We offer the
following replacement which can be used to check that the relevant
subpackage is available (and that an old spatstat is no longer
installed). Note that this will only pass the R package checker if
`spatstat` is in the list of Suggested packages.

Example:

``` r
if(check_spatstat("spatstat.geom")){
      # Run spatstat code
}
```

Function definition:

``` r
check_spatstat <- function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE)){
    stop(paste("Package", pkg, "required; please install it (or the full spatstat package) first."))
  } else{
    spst_ver <- try(packageVersion("spatstat"), silent = TRUE)
    if(!inherits(spst_ver, "try-error") && spst_ver < 2.0-0){
      stop(paste0("You have an old version of spatstat installed which is incompatible with ", 
                 pkg, 
                 ". Please update spatstat (or uninstall it)."))
    }
  }
}
```
