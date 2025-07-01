matchProbePair <- function(Fprobe, Rprobe, subject, algorithm="auto",
                           logfile=NULL, verbose=FALSE, fixed = FALSE)
{
  ## This won't copy the data if Fprobe and Rprobe are already DNAString objects
  F <- DNAString(Fprobe)
  R <- DNAString(Rprobe)
  
  ## F and R hits on the + strand
  Fp_hits <- start(matchPattern(F, subject, algorithm=algorithm, fixed = fixed))
  Rp_hits <- start(matchPattern(R, subject, algorithm=algorithm, fixed = fixed))
  
  ## F and R hits on the - strand
  Fm_hits <- end(matchPattern(reverseComplement(F), subject, algorithm=algorithm, fixed = fixed))
  Rm_hits <- end(matchPattern(reverseComplement(R), subject, algorithm=algorithm, fixed = fixed))
  
  if (verbose) {
    cat("Fp_hits:", Fp_hits, "  Rp_hits:", Rp_hits,
        "  Fm_hits:", Fm_hits, "  Rm_hits:", Rm_hits, "\n")
  }
  
  matches0 <- Biostrings:::reduceProbePairMatches(c(Fp_hits, Rp_hits), c(Fm_hits, Rm_hits))
  ans <- Views(subject, start=matches0$start, end=matches0$end)
  
  if (!is.null(logfile)) {
    nFp <- length(Fp_hits)
    nRp <- length(Rp_hits)
    nFm <- length(Fm_hits)
    nRm <- length(Rm_hits)
    nmatches0 <- length(ans)
    ## cat("", ..., sep="\t") is a trick to get an extra tab
    cat("", nFp, nRp, nFm, nRm, nmatches0, file=logfile, sep="\t")
  }
  ans
}