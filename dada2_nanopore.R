list.of.packages <- c("dada2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dada2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "otu_sequences.fa"
}

path <- args[1]
#fnFs <- list.files(path, full.names = TRUE)
fnFs <- args[1]
fnFs_filtered <- args[3]
filF <- filterAndTrim(
  fnFs,
  filt = fnFs_filtered,
  maxEE = args[4],
  truncLen = 0,
  maxN = 0,
  truncQ = args[5],
  trimRight = 0,
  minQ = 0,
  rm.lowcomplex = 0,
  orient.fwd = NULL,
  matchIDs = FALSE,
  id.sep = "\\s",
  id.field = NULL,
  n = 1e+05,
  OMP = TRUE,
  qualityType = "Auto",
  trimLeft = 0,
  minLen = 50,
  maxLen = Inf,
  rm.phix = FALSE,
  compress = FALSE
)
errF <- learnErrors(
  fnFs_filtered,
  errorEstimationFunction = PacBioErrfun,
  nbases = 1e8,
  nreads = NULL,
  randomize = TRUE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  multithread = TRUE
)
dadaFs <- dada(
  fnFs_filtered,
  err = errF,
  selfConsist = FALSE,
  priors = character(0),
  DETECT_SINGLETONS = FALSE,
  GAPLESS = TRUE,
  GAP_PENALTY = -8,
  GREEDY = TRUE,
  KDIST_CUTOFF = 0.42,
  MATCH = 5,
  MAX_CLUST = 0,
  MAX_CONSIST = 10,
  MIN_ABUNDANCE = 1,
  MIN_FOLD = 1,
  MIN_HAMMING = 1,
  MISMATCH = -4,
  OMEGA_A = 1e-40,
  OMEGA_C = 1e-40,
  OMEGA_P = 1e-4,
  PSEUDO_ABUNDANCE = Inf,
  PSEUDO_PREVALENCE = 2,
  SSE = 2,
  USE_KMERS = TRUE,
  USE_QUALS = TRUE,
  VECTORIZED_ALIGNMENT = TRUE,
  BAND_SIZE = 16,
  HOMOPOLYMER_GAP_PENALTY = NULL,
  pool = FALSE,
  multithread = TRUE
)
i<-0
for(seq in dadaFs$sequence){
i<-i+1
write(c(paste(">Otu",toString(i),sep="")),file=args[2],append=TRUE)
write(c(paste(seq)),file=args[2],append=TRUE)
}

