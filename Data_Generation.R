# This code generates the data used for POD_Tree_Schrodinger.R

library(ReacTran)

# ---------------
# ---------------
# schrodinger equation
# ---------------
# ---------------
Schrodinger <- function(t, u, parms) {
  du <- 1i * tran.1D (C = u, D = 1, dx = xgrid)$dC +
    1i * gam * abs(u)^2 * u
  list(du)
}

alf.range = seq(0.05,0.5,0.01)
u.list = list()
ii = 1

for (i in 1:length(alf.range)){
  alf <- alf.range[i]
  gam <- 1
  N <- 300
  xgrid <- setup.grid.1D(-20, 80, N = N)
  x <- xgrid$x.mid
  c1 <- 1
  c2 <- 0.1
  sech <- function(x) 2/(exp(x) + exp(-x))
  soliton <- function (x, c1)
    sqrt(2*alf/gam) * exp(0.5*1i*c1*x) * sech(sqrt(alf)*x)
  yini <- soliton(x, c1) + soliton(x-25, c2)
  times <- seq(0, 40, by = 0.1)
  
  out <- ode.1D(y = yini, parms = NULL, func = Schrodinger,
                  times = times, dimens = 300, method = "adams")
  
  a = image(abs(out), grid = x, mfrow = NULL, ylab = "Distance, x",
        main = paste("amplitude","=", alf,sep="")  )
  print(a)
  
  u.list[[ii]] = t( abs(out) )[-1,]
  ii = ii + 1
}

save(u.list, file="u.list.schrodinger.RData") # save the data






