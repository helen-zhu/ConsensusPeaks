
remove.max.gaps = function(seg.gr, max.gaps){
  mgap.gr = GenomicRanges::GRanges(seqnames = geneinfo$chr, IRanges::IRanges(start = max.gaps$Var1, end = max.gaps$Var2), strand = geneinfo$strand)
  seg.gr = GenomicRanges::setdiff(seg.gr, mgap.gr)
  seg.gr
}
