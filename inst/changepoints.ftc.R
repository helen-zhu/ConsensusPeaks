
library(ConsensusPeaks)

# Preamble ----------------------------------------------------------------
# FTC function seems to depend heavily on which changepoints are identified
# This was implemented mainly to speed up the computations

# Preliminary -------------------------------------------------------------

# Simulating Histograms of Bed Region Pile-ups ----------------------------

gtf.path = system.file("extdata", package = "ConsensusPeaks")
gtf = paste0(gtf.path, "/test.gtf")

peaks = simulate.gaussian.peaks(
  mu = c(20, 50, 150),
  sd = c(5, 10, 30),
  extend.width = c(10, 30, 20),
  nsamples = c(15, 20, 40),
  gene = "ENSGXX",
  gtf = gtf,
  seed = 123
)

# Convert Peaks to Histogram

# Calculating Changepoints

# Running FTC

# Plotting Results

# Running FTC -------------------------------------------------------------



# Testing this
filename = paste0("~/figures/local.minmax.", gene, ".pdf")
pdf(filename)
opar = par(mfrow = c(3,1), mar = c(2,2,2,2))

# p.init - all the endpoints of the bed files
plot(hist, type = "h")
points(p.init, hist[p.init], col = "orange")
p.init.ftc = ftc.helen(hist, p.init, eps)
points(p.init.ftc, hist[p.init.ftc], col = "orange", pch = 16)

# find.changepoints - all the endpoints + the points right before
plot(hist, type = "h")
points(chg.pts, hist[chg.pts], col = "red")
chg.pts.ftc = ftc.helen(hist, chg.pts, eps)
points(chg.pts.ftc, hist[chg.pts.ftc], col = "red", pch = 16)

# min max
plot(hist, type = "h")
min.max = c(1, sort(unlist(local.minmax(hist))))
points(min.max, hist[min.max], col = "blue")
min.max.ftc = ftc.helen(hist, min.max, eps)
points(min.max.ftc, hist[min.max.ftc], col = "blue", pch = 16)
dev.off()


# New ---------------------------------------------------------------------
# Using a combination of the changpoints and min-max to speed things up
