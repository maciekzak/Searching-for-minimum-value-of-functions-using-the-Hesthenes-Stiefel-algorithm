#Optymalizacja nieliniowa projekt2


install.packages('numDeriv')
library(numDeriv)

install.packages('rgl')
library(rgl)

#funkcja ternary uzywana jest jako krok w funkcji hestenes-stiefel
ternary <- function(f, lower, upper, tol) {
  f.lower <- f(lower)
  f.upper <- f(upper)
  while (abs(upper - lower) > 2 * tol) {
    x1 <- (2 * lower + upper) / 3
    f.x1 <- f(x1)
    x2 <- (lower + 2 * upper) / 3
    f.x2 <- f(x2)
    if (f.x1 < f.x2) {
      upper <- x2
      f.upper <- f.x2
    } else {
      lower <- x1
      f.lower <- f.x1
    }
  }
  return((upper + lower) / 2)
}

#implementacja funkcji znajdujacej minimum
hestenes.stiefel <- function(f, x, tol) {
  beta <- 1
  d <- -grad(f,x)
  iter <- 0
  repeat {
    g <- function(a) {f(x + a * d)}
    step <- ternary(g,0,5,tol)
    new.x <- x + step * d
    if (dist(rbind(new.x,x)) < tol) {
    cat("Wektor Kierunku: ",d, "\nWielkość kroku: ", step, 
        "\nBeta: ",beta,"\nIteracja: ", iter, "\nPunkt: ",new.x,"\n\n")
      return(new.x)
    }
    beta <- (grad(f,new.x) %*% (grad(f,new.x) - grad(f,x))) / (d %*% (grad(f,new.x) - grad(f,x)))
    d <- -grad(f,new.x) + as.vector(beta) * d
    x <- new.x
    iter <- iter + 1
    cat("Wektor Kierunku: ",d, "\nWielkość kroku: ", step,
        "\nBeta: ",beta,"\nIteracja: ", iter, "\nPunkt: ",new.x,"\n\n")
  }
}



funkcja<- function(x) {
  return((1.5 - x[1] + x[1]*x[2])^2 + (2.25 - x[1] + x[1]*x[2]^2)^2 + (2.625 - x[1] + x[1]*x[2]^3)^2)
}

tolerancja <- 1e-14

punk_poczatkowy <- c(-1.5,-1)

#wywolanie funkcji
hestenes.stiefel(funkcja, punk_poczatkowy, tolerancja)



#WIZUALIZACJA
#############
# Podaję Funkcję
Beale <- function(x, y) {
  (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
}


#Podaję zakres
x <- seq(-10, 10, by = 0.01)
y <- seq(-10, 10, by = 0.01)


# Obliczam wartosci
z <- outer(x, y,Beale )

# W zaleznosci od wartosci z dziele na przedzialy by uzyskac ladniejszy kolor
intervals <- cut(z, breaks = 100)

colors <- rainbow(length(levels(intervals)))[as.numeric(intervals)]


# Tworze wykres
persp3d(x, y, z, col = colors)
#dodaje zrodlo swiatla zeby wykres byl jasniejszy
rgl.light(x=2, y=2, z=2, ambient="white", diffuse="white")


minimum.point <- c(3,0.5,0)
rgl.spheres(minimum.point, radius=0.01,color="black")


# #Podaję zakres
x <- seq(2.8, 3.2, by = 0.01)
y <- seq(0.4, 0.6, by = 0.01)


