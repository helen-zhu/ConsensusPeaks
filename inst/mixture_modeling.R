library(ConsensusPeaks)
library(rstan)

# Question: Can we model it directly, if given the "correct" classes?

g1 <- rgamma(100, shape = 2, rate = 1/2)
n1 <- rnorm(100, mean = 15)
x <- obs.to.int.hist(c(g1, n1))
plot(x)
s = unlist(ftc.helen(x))
plot.segments(x, s)


set.seed(13)
u1 <- runif(300, max = 5)
x2 <- c(u1, n1)
x2.hist <- obs.to.int.hist(x2)
plot(x2.hist)
s2 = unlist(ftc.helen(x2.hist))
plot.segments(x2.hist, s2)

N <- length(x2)
y <- x2
input_data = list(N = N, y = y)
mixture_fit <- stan(file='models/mixture_model.stan', data=input_data, chains=4, iter = 1e3)
mixture_fit
