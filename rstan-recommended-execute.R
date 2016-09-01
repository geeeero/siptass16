# create configuration file
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function", 
    file = M, sep = "\n", append = TRUE)

#install.packages('rstan', repos = 'https://cloud.r-project.org/', dependencies=TRUE)
#install.packages("rstan", type = "source")
install.packages("lattice")
install.packages("coda")
install.packages("shinystan")

# verify your toolchain
fx <- inline::cxxfunction( signature(x = "integer", y = "numeric" ) , '
    return ScalarReal( INTEGER(x)[0] * REAL(y)[0] ) ;
' )
fx( 2L, 5 ) # should be 10

# recommended in startup message
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#