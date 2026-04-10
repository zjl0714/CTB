pkgname <- "CTTB"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "CTTB-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CTTB')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CTTB")
### * CTTB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CTTB
### Title: Covariate-tightened trimming bounds
### Aliases: CTTB

### ** Examples

## Not run: 
##D data(simData)
##D res <- CTTB(
##D   data = simData,
##D   Y = "Y", D = "D", S = "S",
##D   X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
##D   Pscore = "Ps",
##D   aggBounds = TRUE,
##D   LeeBounds = TRUE
##D )
##D print(res)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CTTB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.CTTB")
### * plot.CTTB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.CTTB
### Title: Plot CTTB results
### Aliases: plot.CTTB

### ** Examples

## Not run: 
##D data(simData)
##D res <- CTTB(simData, Y = "Y", D = "D", S = "S",
##D             X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
##D             Pscore = "Ps", aggBounds = TRUE, LeeBounds = TRUE)
##D plot(res)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.CTTB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.CTTB")
### * print.CTTB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.CTTB
### Title: Print a CTTB object
### Aliases: print.CTTB

### ** Examples

## Not run: 
##D data(simData)
##D res <- CTTB(simData, Y = "Y", D = "D", S = "S",
##D             X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
##D             Pscore = "Ps", aggBounds = TRUE, LeeBounds = TRUE)
##D print(res)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.CTTB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.CTTB")
### * summary.CTTB

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.CTTB
### Title: Summarize a CTTB object
### Aliases: summary.CTTB

### ** Examples

## Not run: 
##D data(simData)
##D res <- CTTB(simData, Y = "Y", D = "D", S = "S",
##D             X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
##D             Pscore = "Ps", aggBounds = TRUE, LeeBounds = TRUE)
##D summary(res)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.CTTB", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
