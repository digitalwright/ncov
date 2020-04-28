library(DECIPHER) # MUST INSTALL DECIPHER FIRST

# need to set the working directory
setwd("<<path to ncov directory>>/") # needs trailing slash

###
# Load files
###

source('./RobustRegression_v1.R')
source('./map_v1.R')
source('./movavg_v1.R')
source('./coordinates_v1.R')
ref <- readDNAStringSet("./NC_045512.2.fas")
dna <- readDNAStringSet("./gisaid_cov2020_sequences-Apr14.fasta.gz") # expect warning for ignored spaces

cols <- c(A="#C53932", C="#56BBCC", G="#529D3F", T="#D57FBF", U="#F08536")

# prepare genomes and record date of collection
dna <- RemoveGaps(dna)
ns <- strsplit(names(dna), "|", fixed=TRUE)
date <- sapply(ns, tail, n=1)
date <- as.Date(date, "%Y-%m-%d")
min_date <- min(date, na.rm=TRUE)
day <- date - min_date
sum(!is.na(day)) # number with collection date

###
# Create Figure 1
###

muts <- rep(NA_real_, length(dna))
pBar <- txtProgressBar(style=3)
for (i in seq_along(dna)) {
	DNA <- AlignProfiles(dna[i], ref)
	d <- DistanceMatrix(DNA,
		type="dist",
		verbose=FALSE,
		penalizeGapLetterMatches=FALSE)
	t <- TerminalChar(DNA)
	muts[i] <- d[1]*(width(DNA)[1] - max(t[, 1]) - max(t[, 2]))
	setTxtProgressBar(pBar, i/length(dna))
}

dev.new(height=3.5, width=4)
p <- par(mar=c(4, 4, 2, 1))
plot(day,
	muts,
	ylab="Differences from reference (Δ)",
	xlab=paste0("Days since first genome (",
		min_date,
		")"),
	pch=NA)
rect(day - 0.5,
	muts - 0.5,
	day + 0.5,
	muts + 0.5,
	col="#22222222",
	border=NA)
ab <- lm(muts~unclass(day))
abline(ab, lwd=2)
text(0,
	(par("usr")[4] - par("usr")[3])*0.9 + par("usr")[3],
	paste0("Δ ≈ ",
		round(coef(ab)[2], 3),
		"*days + ",
		round(coef(ab)[1], 1)),
	pos=4)
par(p)

###
# Map independent substitutions (Fig. 2b)
###

DNA <- AlignSeqs(dna,
	processors=NULL)
DNA2 <- AlignProfiles(ref, DNA)

d <- DistanceMatrix(DNA,
	type="dist",
	penalizeGapLetterMatches=FALSE,
	processors=NULL)
c <- IdClusters(d,
	method="NJ",
	processors=NULL,
	type="dendrogram",
	collapse=-1)
c <- reorder(c, rev(unlist(c)))
o <- unlist(c)

dev.new(width=8.5, height=7)
p <- par(mar=c(0.1, 0.1, 10, 30))
plot(c, leaflab="none", yaxt="n", horiz=TRUE)

dx <- par("xaxp")
dx <- dx[1] - dx[2]
dy <- par("yaxp")
dy <- dy[2] - dy[1]
x <- par("pin")
segments(0.001,
	0.98*dy,
	0.001,
	0.98*dy + 2e-4/dx*dy*x[1]/x[2])
text(0.001,
	1.093*dy,
	"0.0002",
	pos=2,
	srt=90)

delta <- -6e-5
space <- -2e-5
min_width <- length(o)*0.001 # set min_width = 1 to max all rectanges exact below leaves, but might be invisible because too small to render
sep <- 1 # minimum horizontal separation on the tree
count <- 0L # number of parallel changes
ref_pos <- 0L # position in reference sequence
l <- width(DNA2)[1]

results <- data.frame(ref_pos=integer(),
	reps=integer(),
	sub=character(),
	nsyn=integer(),
	syn=integer(),
	gene=character(),
	aa=character(),
	stringsAsFactors=FALSE)

subs <- matrix(0L,
	nrow=4, # from base
	ncol=4, # to base
	dimnames=list(DNA_BASES,
		DNA_BASES))

pBar <- txtProgressBar(style=3)
for (i in seq_len(l)) {
	setTxtProgressBar(pBar, i/l)
	s <- subseq(DNA2, i, i)
	if (s[1]=="-")
		next # gap in ancestor
	ref_pos <- ref_pos + 1L
	
	# number of cases of parallel mutation
	w1 <- which(s[-1][o] %in% DNA_BASES)
	w2 <- which(s[-1][o][w1] != s[1])
	if (length(w2)==0) {
		reps <- 0
	} else {
		# need to resolve substitutions on the tree
		reps <- map(tree=c,
			bases=as.character(s[-1]),
			ancestor=as.character(s[1]))
	}
	
	if (reps > 0) {
		w <- which(s[-1][o] != s[1] &
			s[-1][o] %in% DNA_BASES)
		start <- tail(cds[cds <= ref_pos], n=1)
		
		sub <- paste0(ifelse(s[1]=="T", "U", s[1]),
					ref_pos,
					names(IUPAC_CODE_MAP)[IUPAC_CODE_MAP==paste(sort(unique(s[-1][o][w1][w2])), collapse="")])
		
		if (length(start) > 0 &&
			ref_pos - start < 3) { # coding region
			start <- start + i - ref_pos
			codons <- subseq(DNA2[c(1,
					seq_along(s)[-1][o][w])],
				start,
				start + 2)
			codons <- as.character(codons)
			aas <- GENETIC_CODE[codons]
			aas <- aas[!is.na(aas)]
			ref_aa <- aas[1]
			aas <- table(aas[-1])
			syn <- aas[ref_aa]
			if (is.na(syn)) {
				syn <- FALSE
			} else {
				syn <- TRUE
			}
			nsyn <- names(aas[names(aas) != ref_aa])
			
			gene <- tail(protein_s[ref_pos - protein_s >= 0], n=1)
			if (ref_pos >= 13468 && # after frameshift
				ref_pos <= 16236) # within nsp12 still
				gene <- gene - 1L
			
			results <- rbind(results,
				data.frame(ref_pos=ref_pos,
					reps=as.integer(reps),
					sub=sub,
					nsyn=as.integer(length(nsyn)),
					syn=as.integer(syn),
					gene=names(gene),
					aa=paste0(ref_aa,
						ceiling((ref_pos - gene + 1)/3),
						paste(nsyn,
							collapse="/")),
					stringsAsFactors=FALSE))
			nsyn <- length(nsyn)
		} else { # non-coding region
			codons <- NULL
			results <- rbind(results,
				data.frame(ref_pos=ref_pos,
					reps=as.integer(reps),
					sub=sub,
					nsyn=NA_integer_,
					syn=NA_integer_,
					gene=NA_character_,
					aa=NA_character_,
					stringsAsFactors=FALSE))
		}
		
		fromto <- cbind(as.character(s[1]),
			as.character(unique(s[-1][o][w1][w2])))
		subs[fromto] <- subs[fromto] + 1L
		
		if (reps >= 7 && # evidence of parallelism is >= N
			ref_pos > 33 && # 5'UTR > 265
			ref_pos < 29675) { # 3'UTR < 29675
			count <- count + 1L
			rect(delta*(count - 1) + space*count,
				seq_along(DNA)[w] - min_width/2,
				delta*count + space*count,
				seq_along(DNA)[w] + min_width/2,
				col=cols[ifelse(as.character(s)[-1][o]=="T",
					"U",
					as.character(s)[-1][o])][w],
				border=NA,
				xpd=TRUE)
			vals <- sort(unique(s[-1][o][w1][w2]))
			lab <- paste0(ifelse(s[1]=="T", "U", s[1]),
					ref_pos,
					paste(ifelse(vals=="T", "U", vals),
						collapse="/"),
					"(",
					reps,
					")")
			if (!is.null(codons)) {
				if (syn) {
					if (nsyn==0) {
						lab <- paste0(lab,
							" S; ")
					} else {
						lab <- paste0(lab,
							" S/N; ")
					}
				} else {
					lab <- paste0(lab,
						" N; ")
				}
				lab <- paste0(lab,
					tail(results$gene, n=1),
					" (",
					tail(results$aa, n=1),
					")")
			}
			text(x=delta*(count - 1.2) + space*count,
				y=length(o)*1.01,
				labels=lab,
				pos=4,
				cex=0.8,
				srt=90,
				xpd=TRUE)
		}
	}
}
par(p)

#write.csv(results,
#	file="./results_v3.csv")

###
# Compute statistics including dN/dS
###

# total number of mutations
nrow(results)

# number of non-coding mutations
sum(!is.na(results$syn))

# number of independent mutations
sum(results$reps > 1)

# number independent in non-coding regions
sum(!is.na(results$syn[results$reps > 1]))

# N and S overall
sum(results[, "nsyn"], na.rm=TRUE) # N
sum(results[, "syn"], na.rm=TRUE) # S

# N and S on parallel
sum(results[results$reps > 1, "nsyn"], na.rm=TRUE) # N
sum(results[results$reps > 1, "syn"], na.rm=TRUE) # S

###
# Determine expected frequency of N and S
###

a <- alphabetFrequency(ref, as.prob=TRUE)[1:4] # base frequency
subs # substitutions [from, to]

t <- subseq(rep(ref, length(protein_s)),
	protein_s,
	protein_e)
t <- trinucleotideFrequency(t)
t <- colSums(t[, -12]) # drop nsp12 because frame-shifted
t <- t/sum(t)

p <- subs/sum(subs)
N <- 0
S <- 0
for (i in seq_along(t)) {
	s1 <- names(t)[i]
	s <- strsplit(s1, "")[[1]]
	for (j in seq_len(4)) {
		for (k in 1:3) {
			s2 <- s
			s2[k] <- colnames(p)[j]
			s2 <- paste(s2, collapse="")
			if (GENETIC_CODE[s2]==GENETIC_CODE[s1]) { # synonmyous
				S <- S + t[i]*p[s[k], colnames(p)[j]]
			} else { # non-synonymous
				N <- N + t[i]*p[s[k], colnames(p)[j]]
			}
		}
	}
}
S
N
N/(S + N) # fraction N expected

###
# Draw Figure 2a
###

arrow <- 300 # maximum arrow length
h <- 0.1 # height of gene bars

dev.new(width=7.3, height=2)
layout(matrix(1:2), heights=c(1.5, 0.5))
p <- par(mar=c(0, 4, 1, 4))

plot(NA,
	xlim=c(0, width(ref)),
	ylim=c(0, max(results$reps)),
	ylab="",
	xlab="",
	xaxt="n",
	bty="n")

# add conserved regions
zeros <- rep(0, width(ref))
zeros[results$ref_pos] <- 1 # mutations
zeros <- movavg(zeros)
zeros <- rle(zeros < 0.06) # absence of mutations
ranges <- cumsum(zeros$lengths)
zeros <- which(zeros$values &
	zeros$lengths >= 100)
rect(ranges[zeros - 1],
	0,
	ranges[zeros],
	par("usr")[2],
	col="#F9D8B4",
	border=NA)
for (i in seq_along(zeros)) {
	begin <- ranges[zeros[i] - 1]
	end <- ranges[zeros[i]]
	w <- which(protein_s < begin &
		protein_e > begin)
	cat("\n",
		i,
		"\t",
		begin,
		"\t",
		end,
		"\t",
		names(w),
		"\t",
		ceiling((begin - protein_s[w] + 1)/3),
		"\t",
		ceiling((end - protein_s[w] + 1)/3),
		sep="")
}
# NOTE: Can go out of bounds in the protein sequence

# add dn/ds
dnds <- avg_dnds <- rep(NA_real_, width(ref))
dnds[results$ref_pos] <- ifelse(is.na(results$nsyn) | results$nsyn==0, 0, 1) - results$syn
for (i in seq_along(protein_s)) {
	index <- protein_s[i]:protein_e[i]
	avg_dnds[index] <- movavg(dnds[index], 0.01)
}
ylims <- axTicks(2)
ylims <- c(ylims[1], ylims[length(ylims)])
mymax <- ceiling(max(abs(avg_dnds), na.rm=TRUE)*50)/50
scaled_dnds <- avg_dnds/mymax # normalize avg_dnds to mymax
scaled_dnds <- scaled_dnds*diff(range(ylims))/2
scaled_dnds <- scaled_dnds + mean(ylims) # shift mean to center-line
lines(scaled_dnds,
	col="#AAAAAA")
ticks <- seq(-1, 1, 0.5)
ticks <- ticks*diff(range(ylims))/2
ticks <- ticks + mean(ylims)
axis(4,
	at=ticks,
	labels=c(-mymax, "", "0", "", mymax),
	col="#AAAAAA",
	col.ticks="#AAAAAA",
	col.axis="#AAAAAA")
mtext("Moving average\n(number of N - S)",
	4,
	line=3,
	las=3,
	col="#AAAAAA")
segments(0,
	ticks[3],
	width(ref),
	ticks[3],
	lty=2,
	col="#AAAAAA")

# add mutations
title(ylab="Independent\nsubstitutions",
	line=2)
segments(results$ref_pos,
	0,
	results$ref_pos,
	results$reps,
	col=ifelse(is.na(results$nsyn),
		"#D57FBF", # pink = non-coding
		ifelse(results$nsyn > 0,
			ifelse(results$syn > 0,
				"#BDBC35", # yellowish = both
				"#3976AF"), # blue = non-synonymous
			"#85584E"))) # brown = synonymous

par(mar=c(0, 4, 0, 4))
plot(NA,
	xlim=c(0, width(ref)),
	ylim=c(-2.3*h, h),
	xlab="",
	ylab="",
	yaxt="n",
	xaxt="n",
	bty="n")
segments(1,
	h/2,
	width(ref),
	h/2,
	lwd=4)
for (i in seq_along(bound_s)) {
	a <- min(bound_e[i] - bound_s[i],
		arrow)
	rect(bound_s[i],
		0,
		bound_e[i] - a + 1,
		h,
		col="#8D6BB8",
		border=NA)
	polygon(c(bound_e[i],
		bound_e[i] - a,
		bound_e[i] - a),
		c(h/2, 0, h),
		col="#8D6BB8",
		border=NA)
}
w <- which(!(protein_s %in% bound_s))
segments(protein_s[w],
	0,
	protein_s[w],
	h)
w <- which(!(protein_e %in% bound_e))
w <- w[-length(w)] # drop final end because only a few positions from the protein terminus
segments(protein_e[w],
	0,
	protein_e[w],
	h)
xs <- (protein_s + protein_e)/2
xs <- xs + 540 # center the labels
w <- which(names(xs) %in% c("nsp7", "nsp9", "nsp11", "ORF6", "ORF7b"))
text(x=xs[-w],
	y=-h/10,
	labels=names(protein_s)[-w],
	pos=2,
	cex=0.6,
	srt=90)
par(p)
