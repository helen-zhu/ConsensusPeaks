
remove.max.gaps = function(geneinfo, seg.gr, max.gaps){
  if(nrow(max.gaps) > 0){

    # Creating a GenomicRanges Object with replicates of max gap
    mgap.gr = GenomicRanges::GRanges(seqnames = geneinfo$chr, IRanges::IRanges(start = max.gaps$Var1, end = max.gaps$Var2), strand = geneinfo$strand)
    mgap.id = sort(rep(1:length(seg.gr), length(mgap.gr)))
    mgap.gr = rep(mgap.gr, length(seg.gr))
    mgap.gr = split(mgap.gr, mgap.id)

    # Identifying the difference in set size
    seg.gr = split(seg.gr, 1:length(seg.gr))
    seg.gr = GenomicRanges::setdiff(seg.gr, mgap.gr)
    seg.gr = unlist(seg.gr)
  }
  seg.gr
}
