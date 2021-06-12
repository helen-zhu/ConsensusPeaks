library(ConsensusPeaks)
library(BoutrosLab.plotting.general)

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

find.all.meaningful.gap = function(x) {
  span = seq_along(x)
  todo = expand.grid(span, span)
  todo = todo[todo$Var2 > todo$Var1,]
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
  df = cbind(todo, mgap)
  # df = df[order(df$Var2, df$Var1),]
  df = df[order(df$entropy), ]
  df$mgap = as.numeric(df$mgap)
  # df$scaled_entropy = (df$entropy - min(df$entropy, na.rm = T)) / (max(df$entropy, na.rm = T) - min(df$entropy, na.rm = T))

  seg.gap.data = df[df$mgap > 0 & !is.na(df$mgap), ]
  # seg.gap.data$index = 1:nrow(seg.gap.data)

  maximal.meaningful(seg.gap.data)
}

# Finds the meaningful gaps between the points in s
meaningful.gaps.local = function(x, seg.points) {

  max.gaps.list <- lapply(seq(2, length(seg.points)), function(i) {
    x.sub = x[seg.points[i-1]:seg.points[i]]

    max.gaps = find.all.meaningful.gap(x.sub)
    if(nrow(max.gaps) > 0) {
      max.gaps[, c('Var1','Var2')] <- max.gaps[, c('Var1','Var2')] + seg.points[i-1]
      max.gaps$seg.start = seg.points[i-1]
      max.gaps$seg.end = seg.points[i]
      max.gaps
    }
  })

  max.gaps <- do.call(rbind.data.frame, max.gaps.list)
  max.gaps
}

find.new.segments = function(gaps.df) {
  new.segments = c()
  for(i in rownames(gaps.df)) {
    r = gaps.df[i, ]
    # Have a large gap
    if(r$Var2 - r$Var1 > 1) {
      left.diff = abs(r$seg.start - r$Var1)
      right.diff = abs(r$seg.end-r$Var2)
      if(left.diff < right.diff && right.diff > 2) {
        new.segments = c(new.segments, r$Var2 - 1)
      } else if(left.diff > right.diff && left.diff > 2){
        new.segments = c(new.segments, r$Var1 + 1)
      }
    }
  }
  new.segments
}

maximal.gaps.multiplot = function(x, seg.points, max.gaps, new.segments = NULL, ...) {
  # Add the gap data to a vector for heatmap
  maximal.data <- rep("none", length(x))
  for(i in rownames(max.gaps)) {
    r = max.gaps[i,]
    gap.seq = seq(r$Var1, r$Var2) - 1
    maximal.data[gap.seq] <- "gap"
  }
  maximal.data <- factor(maximal.data, levels = c("none", "gap"))

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

  colour.scheme = c('white', 'orange')
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

  segment.data <- rep("none", length(x))
  segment.data[seg.points] <- "segment"
  if(is.null(new.segments)) {
    colour.scheme.segments <- c('white', 'red')
    segment.data <- factor(segment.data, levels = c("none", "segment"))
  } else {
    colour.scheme.segments <- c('white', 'red', 'blue')
    segment.data[new.segments] <- "extra_segment"
    segment.data <- factor(segment.data, levels = c("none", "segment", "extra_segment"))
  }

  segments.heatmap <- create.heatmap(
    x = t(as.matrix(as.numeric(segment.data))),
    clustering.method = 'none',
    scale.data = FALSE,
    colour.scheme = colour.scheme.segments,
    total.colours = 1 + length(colour.scheme.segments),
    grid.col = TRUE,
    print.colour.key = FALSE,
    force.grid.col = TRUE,
    axes.lwd = 0
  )

  create.multipanelplot(
    plot.objects = list(bp, segments.heatmap, maximal.heatmap),
    plot.objects.heights = c(1, 0.15, 0.15),
    y.spacing = -0.75,
    ylab.label = "Segments",
    xlab.label = "Index",
    ylab.cex = 2,
    xlab.cex = 2,
    height = 6,
    width = 8,
    ...
  )
}

set.seed(5)
n1 = rnorm(100, mean = 5, sd = 2)
n2 = rnorm(50, mean = 10, sd = 3)
n3 = rnorm(200, mean = 30, sd = 1)
n4 = rnorm(150, mean = 50, sd = 1)
n.combined = c(n1, n2, n3, n4)
x = obs.to.int.hist(n.combined)
plot.segments(x)

# Running FTC
s = unlist(ftc.helen(x, s = sort(unlist(local.minmax(x)))))
plot.segments(x, s)

mgaps = meaningful.gaps.local(x, s)


filename = "plots/Meaningful.Segments.png"
maximal.gaps.multiplot(x, s, mgaps, filename = filename)

set.seed(13)
g1 <- rgamma(100, shape = 2, rate = 1/2)
n1 <- rnorm(100, mean = 15)
u1 <- runif(300, max = 5)

x2 <- c(u1, n1)
tab.x <- obs.to.int.hist(x2)
s2 = unlist(ftc.helen(tab.x))
plot.segments(tab.x, s2)

mgaps2 = meaningful.gaps.local(tab.x, s2)

new.segments = find.new.segments(mgaps2)
filename = "plots/Meaningful.Segments.unif.png"
maximal.gaps.multiplot(tab.x, s2, mgaps2, new.segments, filename = filename)

mod = list()
fits = list()
mod.optim = list()
jc.optim.multi.distr = list()
for(i in seq(1, length(s2) - 1)) {
  seg.start = s2[i]
  seg.end = s2[i + 1]
  x.sub <- x2[x2 >= seg.start & x2 <= seg.end]
  x.adjusted <- (x.sub - seg.start) + 1e-10
  x.range = seg.start:seg.end

  mod[[i]] = fit.continuous.distributions(
    x = as.numeric(x.adjusted[x.range]),
    seg.start = seg.start,
    seg.end = seg.end,
    fit.mixtures = c('unif', 'tnorm', 'tgamma', 'tgamma_flip'),
    max.iterations = 500)

  fits[[i]] = extract.distribution.parameters(
    mod = mod[[i]],
    x = x.adjusted)

  mod.optim[[i]] = which(fits[[i]]$aic == min(fits[[i]]$aic))
  jc.optim.multi.distr[[i]] = fits[[i]]$jc[mod.optim[[i]]]
}
