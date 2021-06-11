library(ConsensusPeaks)
library(isotone)
library(ggplot2)
library(BoutrosLab.plotting.general)
library(reshape2)

# Preamble ----------------------------------------------------------------
# Trying to implement a gap detection thing based on the Delon papers


# Functions ---------------------------------------------------------------

# Take a vector of values and get the histogram for integer breaks
obs.to.int.hist = function(x) {
  table(cut(x, breaks = (floor(min(x)) - 1):(ceiling(max(x)) + 1)))
}

# Relative entropy
rel.entropy = function(h, p, a, b) {
  interval = a:b
  # Round to prevent floating point issues
  hab = round(sum(h[interval]), digits = 14)
  pab = round(sum(p[interval]), digits = 14)
  if(pab == 0 || pab == 1 || hab == 0 || hab == 1) {
    return(0)
  } else {
    hab * log(hab / pab) + (1 - hab) * log((1 - hab) / (1 - pab))
  }
}

# Plots the vector x of counts (or table) and the optional segment points s
plot.segments = function(x, s = NULL, threshold = 0, ...) {
  index = seq_along(x)
  plot(x, type = "h", ...)

  if(!is.null(s)) {
    points(s, x[s], col = "orange")
  }
}

# Uniform Generation
generate.unif = function(x){
  rep(1/length(x), length(x))
}

calc.prob.diff = function(h, p, a, b){
  interval = a:b
  # Round to prevent floating point issues
  hab = round(sum(h[interval]), digits = 14)
  pab = round(sum(p[interval]), digits = 14)
  hab > pab
}

# Meaningful interval
meaningful.interval = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  prob.diff = calc.prob.diff(h, p, a, b)
  c(mint = relative.entropy >= (1/N)*log(L*(L+1)/2) && prob.diff, entropy = relative.entropy)
}

# Meaningful gap
meaningful.gap = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  prob.diff = calc.prob.diff(h, p, a, b)
  mgap = (relative.entropy >= (1/N)*log(L*(L+1)/2) && !prob.diff) || (all(h == 0))
  c(mgap = mgap, entropy = relative.entropy)
}

# Test Cases --------------------------------------------------------------

set.seed(5)
n1 = rnorm(100, mean = 5, sd = 2)
n2 = rnorm(50, mean = 10, sd = 3)
n3 = rnorm(200, mean = 30, sd = 1)
n4 = rnorm(150, mean = 50, sd = 1)
n.combined = c(n1, n2, n3, n4)
x = obs.to.int.hist(n.combined)
plot.segments(x)

# Running FTC
s = ftc.helen(x, s = sort(unlist(local.minmax(x))))
s2 = ftc(x)
plot.segments(x, unlist(s))
plot.segments(x, s2)

# Trying Maximal Meaningful
span = seq_along(x)
todo = expand.grid(span, span)
todo = todo[todo$Var2 > todo$Var1,]
mint = do.call(rbind, lapply(1:nrow(todo), function(i) {
  meaningful.interval(
    h = x/sum(x),
    p = generate.unif(x),
    a = todo$Var1[i],
    b = todo$Var2[i],
    N = sum(x),
    L = length(x)
  )
}))
mgap = do.call(rbind, lapply(1:nrow(todo), function(i) {
  meaningful.gap(
    h = x/sum(x),
    p = generate.unif(x),
    a = todo$Var1[i],
    b = todo$Var2[i],
    N = sum(x),
    L = length(x)
  )
}))
df = cbind(todo, mint, mgap)
# df = df[order(df$Var2, df$Var1),]
df = df[order(df$entropy), ]
df$mint = as.numeric(df$mint)
df$mgap = as.numeric(df$mgap)
df$scaled_entropy = (df$entropy - min(df$entropy)) / (max(df$entropy) - min(df$entropy))

# mint.cond = df$mint > 0
# plotting.mint = matrix(NA, nrow = sum(mint.cond), ncol = 35)
# plotting.mint

# Create A Heatmap Looking For Meaningful Segments
# plotting.mint = acast(df, Var1 ~ Var2, value.var = "mint")
# plotting.mgap = acast(df, Var1 ~ Var2, value.var = "mgap")

bp = create.barplot(
  Freq ~ Var1,
  data.frame(x),
  xaxis.cex = 0,
  xlab.cex = 0,
  ylab.label = "Counts",
  xaxis.tck = 0,
  yaxis.tck = 0,
  yaxis.cex = 0.8,
  ylab.cex = 1
)

seg.int.data = df[df$mint > 0, ]
seg.int.data$index = 1:nrow(seg.int.data)

seg.gap.data = df[df$mgap > 0, ]
seg.gap.data$index = 1:nrow(seg.gap.data)

# Quick and dirty implementation of finding maximal meaningful intervals (or gaps)
maximal.meaningful = function(x) {
  curr.df = x
  max.intervals = data.frame()
  while(nrow(curr.df) > 0) {
    max.entropy.index = which.max(curr.df$entropy)
    max.entropy = curr.df[max.entropy.index, ]
    max.entropy.seq = seq(max.entropy$Var1, max.entropy$Var2)
    max.intervals = rbind(max.intervals, max.entropy)

    # Find all of the segments that overlap.
    # These will all be less than the maximum
    overlap.max = mapply(function(from, to) {
      s = seq(from, to)
      length(intersect(max.entropy.seq, s)) > 0
    }, from = curr.df$Var1, to = curr.df$Var2)

    curr.df = curr.df[!overlap.max, ]
  }
  # Preserve the old index
  max.intervals$index = rownames(max.intervals)
  rownames(max.intervals) <- NULL
  max.intervals
}

max.intervals = maximal.meaningful(seg.int.data)
max.gaps = maximal.meaningful(seg.gap.data)

# Construct data for heatmap
maximal.data <- rep("none", length(x))
for(i in rownames(max.intervals)) {
  r = max.intervals[i,]
  int.seq = seq(r$Var1, r$Var2)
  maximal.data[int.seq] <- "interval"
}

for(i in rownames(max.gaps)) {
  r = max.gaps[i,]
  gap.seq = seq(r$Var1, r$Var2)
  maximal.data[gap.seq] <- "gap"
}
maximal.data <- factor(maximal.data, levels = c("none", "gap", "interval"))

colour.scheme = c('white', 'blue', 'orange')
maximal.heatmap <- create.heatmap(
  x = t(as.matrix(as.numeric(maximal.data))),
  clustering.method = 'none',
  scale.data = FALSE,
  colour.scheme = colour.scheme,
  total.colours = 1 + length(colour.scheme),
  grid.col = TRUE,
  print.colour.key = FALSE,
  force.grid.col = TRUE,
  axes.lwd = 0
)

filename = "plots/Meaningful.Segments.png"
#pdf(filename, height = 12, width = 8)
create.multipanelplot(
  plot.objects = list(bp, maximal.heatmap),
  plot.objects.heights = c(1, 0.1),
  y.spacing = -1,
  ylab.label = "Segments",
  xlab.label = "Index",
  ylab.cex = 2,
  xlab.cex = 2,
  height = 6,
  width = 8,
  filename = filename
)
#dev.off()
