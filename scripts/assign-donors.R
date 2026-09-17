#!/usr/bin/env Rscript

# Confident triplets from three-pool-donors
# Confident doublets from two-pool donors
# Confident doublets from three-pool-donors (third sublibrary is uncertain)
# Any singlets remaining that overlap with one-pool donors
# Assign zero-degrees-of-freedom identity and check remaining matches

library(data.table)
library(foreach)
best_lod_min <- 50
second_best_lod_min <- 50
options(datatable.fwrite.sep = "\t")

# Q: Does any pair of libraries ever have more than one sample in common?
# Stated differently, can a pair of libraries uniquely be assigned to a donor?

# Generates symmetric LOD comparison table and saves it only if it doesn't exist


fingerprints_LODs_fn <- 'DATA/FREEMUXLET/FINGERPRINTS/LOD-scores.tsv'
if(! file.exists(fingerprints_LODs_fn)) {
    file_list <- list.files('DATA/FREEMUXLET/FINGERPRINTS/', pattern='*.txt', full.names=T)

    dat <- foreach(file=file_list, .combine='rbind') %do% {
        fread(file)
    }
    
    
    # Duplicate table with inverse matches
    dat2 <- copy(dat)
    setcolorder(dat2, c('RIGHT_GROUP_VALUE','LEFT_GROUP_VALUE','LOD_SCORE'))
    setnames(dat2, c('LEFT_GROUP_VALUE','RIGHT_GROUP_VALUE','LOD_SCORE'))
    dat <- rbindlist(list(dat, dat2))
    setnames(dat, 'LEFT_GROUP_VALUE', 'SUBLIB1')
    setnames(dat, 'RIGHT_GROUP_VALUE', 'SUBLIB2')
    dat[, subset_pair := paste(sort(c(SUBLIB1,SUBLIB2)), collapse='_'), by=.I]


    dat[, 'LIB1' := tstrsplit(SUBLIB1, split='\\.')[1]]
    dat[, 'LIB2' := tstrsplit(SUBLIB2, split='\\.')[1]]
    dat <- dat[!is.na(LOD_SCORE)][LIB1 != LIB2]

    fwrite(dat, file=fingerprints_LODs_fn, quote=FALSE)
    rm(dat)
    rm(dat2)
}
dat <- fread(fingerprints_LODs_fn)

highLOD <- unique(copy(dat[LOD_SCORE > best_lod_min]))

# Assign reciprocal best match with CheckFingerprints LOD score > 100

reciprocal_matches_fn <- 'DATA/FREEMUXLET/FINGERPRINTS/reciprocal-match-scores.tsv'

if(! file.exists(reciprocal_matches_fn)) {
    lib1s <- copy(highLOD[, .SD, .SDcols=c('SUBLIB1','SUBLIB2','LOD_SCORE')])
    lib2s <- copy(highLOD[, .SD, .SDcols=c('SUBLIB1','SUBLIB2','LOD_SCORE')])

    # Select highest LOD score for each SUBLIB1 or SUBLIB2
    lib1s <- lib1s[lib1s[, .I[which.max(LOD_SCORE)], by = SUBLIB1]$V1]
    lib2s <- lib2s[lib2s[, .I[which.max(LOD_SCORE)], by = SUBLIB2]$V1]

    recip_hits <- merge(lib1s, lib2s, by=c('SUBLIB1','SUBLIB2'), all=F)[SUBLIB1 < SUBLIB2][LOD_SCORE.x == LOD_SCORE.y]
    recip_hits[, 'LOD_SCORE.y' := NULL]
    setnames(recip_hits, 'LOD_SCORE.x','LOD_SCORE')

    recip_hits[, subset_pair := paste(c(SUBLIB1,SUBLIB2), collapse='_'), by=.I]
    recip_hits[, 'LIB1' := tstrsplit(SUBLIB1, split='\\.')[1]]
    recip_hits[, 'LIB2' := tstrsplit(SUBLIB2, split='\\.')[1]]
    recip_hits[, pool_pair := paste0(LIB1, '_', LIB2)]

    fwrite(recip_hits, file=reciprocal_matches_fn, quote=F)
    rm(recip_hits)
}
reciprocal_hits <- fread(reciprocal_matches_fn)
setkey(reciprocal_hits, pool_pair)


combos_fn <- 'DATA/FREEMUXLET/FINGERPRINTS/overlap_donor_combos.tsv'
if(! file.exists(combos_fn)) {
    # Build table of pair-donor combos
    pools <- fread('DATA/syn51753257/AMP-AD_DiverseCohorts_assay_multiome_metadata.csv')
    pools <- pools[, .SD, .SDcols=c('specimenID','libraryBatch')]
    pools[, individualID := tstrsplit(specimenID, split='_')[1]]

    onepool <- pools[individualID %in% pools[, .N, by=individualID][N==1]$individualID]
    twopools <- pools[individualID %in% pools[, .N, by=individualID][N==2]$individualID]
    threepools <- pools[individualID %in% pools[, .N, by=individualID][N==3]$individualID]
    all_libs <- unique(pools$libraryBatch)

    combos <- unique(CJ(L1=all_libs, L2=all_libs)[L1 != L2])[L1 < L2]
    o <- foreach(i=1:nrow(combos), .combine='rbind') %do% {
        lib_1 <- combos[i,L1]
        lib_2 <- combos[i,L2]
        overlap_donors <- intersect(pools[libraryBatch==lib_1]$individualID,
        pools[libraryBatch==lib_2]$individualID)
        data.table('LIB1'=lib_1, 'LIB2'=lib_2, 'overlap_donors'=paste0(overlap_donors, collapse=','))
    }
    o <- o[overlap_donors != '']
    o[, pool_pair := paste0(LIB1, '_', LIB2)]
    fwrite(o, file=combos_fn, quote=F)
    rm(o)
    rm(all_libs)
}
combos <- fread(combos_fn)


# Get table of trios
trios_fn <- 'DATA/FREEMUXLET/FINGERPRINTS/trios.tsv'
if(! file.exists(trios_fn)) {
    # Generate table of expected trios (donors in exactly 3 pools)
    combos2 <- copy(combos)
    setnames(combos2, 'LIB2','LIB3')
    setnames(combos2, 'LIB1','LIB2')
    setkey(combos2, LIB2)
    setkey(combos2, LIB2)
    trios <- merge(combos, combos2, 'LIB2', allow.cartesian=T)[overlap_donors.x == overlap_donors.y]
    trios[, 'overlap_donors.y' := NULL]
    setnames(trios, 'overlap_donors.x', 'overlap_donors')
    setcolorder(trios, c('LIB1','LIB2','LIB3','overlap_donors','pool_pair.x','pool_pair.y'))
    trios[, c('pool_pair.x','pool_pair.y') := NULL]
    fwrite(trios, file=trios_fn, quote=F)
    rm(combos2)
}
trios <- fread(trios_fn)


#===================================================================================================
# TRIPLETS
#===================================================================================================
# For donors that have samples across three pools, it can be expected that these subpools will
# match each other better than any other. To identify these groups of three coming from one donor,
# check the first- and second-best match for every subpool. The three subpools that come from one
# shared donor will, together, have a unique set of only three subpools: each other.

triplet_score_fn <- 'DATA/FREEMUXLET/FINGERPRINTS/triplet-scores.tsv'
if(! file.exists(triplet_score_fn)) {
    
    getCombos <- function(set_of_three) {
        o <- foreach(L1=paste0(set_of_three[1], c('.0','.1','.2')), .combine='rbind') %do% {
                foreach(L2=paste0(set_of_three[2], c('.0','.1','.2')), .combine='rbind') %do% {
                    foreach(L3=paste0(set_of_three[3], c('.0','.1','.2')), .combine='rbind') %do% {
                    return(data.table(L1,L2,L3))
                }
            }
        }
        return(o)
    }
    
    get_triplet_scores <- function(triplet) {
        combo_trios <- getCombos(triplet)
        dat.sub <- dat[LIB1 %in% triplet & LIB2 %in% triplet]
        dat.sub[, rank := frank(-LOD_SCORE), by=SUBLIB1]
        setkey(dat.sub, SUBLIB1, SUBLIB2)
        o <- foreach(L1=combo_trios$L1, L2=combo_trios$L2, L3=combo_trios$L3, .combine='rbind') %do% {
            L1.1 <- dat.sub[list(L1)][rank==1]$SUBLIB2          # best match for sublibrary 1
            L1.1.LOD <- dat.sub[list(L1)][rank==1]$LOD_SCORE    # LOD score of best match
            L1.2 <- dat.sub[list(L1)][rank==2]$SUBLIB2          # second-best match for sublibrary 1
            L1.2.LOD <- dat.sub[list(L1)][rank==2]$LOD_SCORE    # LOD score of second-best match
            L2.1 <- dat.sub[list(L2)][rank==1]$SUBLIB2          # etc
            L2.1.LOD <- dat.sub[list(L2)][rank==1]$LOD_SCORE
            L2.2 <- dat.sub[list(L2)][rank==2]$SUBLIB2
            L2.2.LOD <- dat.sub[list(L2)][rank==2]$LOD_SCORE
            L3.1 <- dat.sub[list(L3)][rank==1]$SUBLIB2
            L3.1.LOD <- dat.sub[list(L3)][rank==1]$LOD_SCORE
            L3.2 <- dat.sub[list(L3)][rank==2]$SUBLIB2
            L3.2.LOD <- dat.sub[list(L3)][rank==2]$LOD_SCORE
            N <- length(unique(c(L1, L2, L3, L1.1, L1.2, L2.1, L2.2, L3.1, L3.2)))
            data.table('triplet'=paste0(triplet,collapse='_'), L1, L2, L3, L1.1,L1.1.LOD,L1.2,L1.2.LOD,L2.1,L2.1.LOD,L2.2,L2.2.LOD,L3.1,L3.1.LOD,L3.2,L3.2.LOD,N)
        }
        return(o)
    }

    unique_trios <- trios[! ( LIB1 %like% '[ABC]$' | LIB2 %like% '[ABC]$' | LIB3 %like% '[ABC]$') ]
    o3 <- foreach(lib1=unique_trios$LIB1, lib2=unique_trios$LIB2, lib3=unique_trios$LIB3, donor=unique_trios$overlap_donors, .combine='rbind') %do% {
        o2 <- get_triplet_scores(c(lib1, lib2, lib3))
        o2[, 'donor' := donor]
        return(o2[])
    }
    fwrite(o3, file=triplet_score_fn, quote=F)
    rm(o3)
    rm(unique_trios)
}
triplet_scores <- fread(triplet_score_fn)
triplet_scores[, c('LIB1','LIB2','LIB3') := tstrsplit(triplet, split='_')]
setkey(triplet_scores, LIB1, LIB2, LIB3)

# Defining good triplets: only three subpools identified (all 3 are reciprocally best matches)
# and best/second-best match LOD > 50

good_triplets <- triplet_scores[N==3][L1.1.LOD > best_lod_min][L2.1.LOD > best_lod_min][L3.1.LOD > best_lod_min][L1.2.LOD > second_best_lod_min][L2.2.LOD > second_best_lod_min][L3.2.LOD > second_best_lod_min]
good_triplets <- good_triplets[, .SD, .SDcols=c('L1','L2','L3','LIB1','LIB2','LIB3','donor')]
setnames(good_triplets, c('L1','L2','L3'), c('SUBLIB1','SUBLIB2','SUBLIB3'))
trios[, triplet := paste0(LIB1, '_', LIB2, '_', LIB3)]
trios[, 'complete' := 0]
trios[triplet %in% good_triplets$triplet, complete := 3]



# Get reciprocal best-matches from remaining three-pool-donors
incomplete_trios <- trios[complete != 3]
incomplete_trios[, pair1 := paste0(LIB1, '_', LIB2)]
incomplete_trios[, pair2 := paste0(LIB1, '_', LIB3)]
incomplete_trios[, pair3 := paste0(LIB2, '_', LIB3)]
incomplete_trio_pairs <- unique(c(incomplete_trios$pair1, incomplete_trios$pair2, incomplete_trios$pair3))
twothirds_triplets <- reciprocal_hits[.(incomplete_trio_pairs)][!is.na(LOD_SCORE)]
twothirds_triplets[, 'LIB3' := NA]
twothirds_triplets[, 'SUBLIB3' := NA]

# Add in overlap donor ID
twothirds_triplets <- merge(combos, twothirds_triplets, by=c('LIB1','LIB2','pool_pair'))[, .SD, .SDcols=c('SUBLIB1','SUBLIB2','SUBLIB3','LIB1','LIB2','LIB3','overlap_donors')]
setnames(twothirds_triplets, 'overlap_donors','donor')

# Double-pool donors
double_pool_donors <- combos[overlap_donors %in% combos[, .N, by=overlap_donors][N==1, overlap_donors]]
double_pools <- merge(double_pool_donors, reciprocal_hits, by=c('pool_pair','LIB1','LIB2'))
double_pools[, 'LIB3' := NA]
double_pools[, 'SUBLIB3' := NA]
setnames(double_pools, 'overlap_donors','donor')
double_pools <- double_pools[, .SD, .SDcols=c('SUBLIB1','SUBLIB2','SUBLIB3','LIB1','LIB2','LIB3','donor')]
# confident about 133/181

good_triplets[, type := 'complete_triplet']
good_triplets           # 163 / 260 Trios completely solved by circular best-match set

twothirds_triplets      # 46 / 97 _remaining_ trios partly solved by reciprocal-best-match
twothirds_triplets[, type := 'partial_triplet']

double_pools            # 133 /181 two-pool donors solved by reciprocal-best-match
double_pools[,type := 'complete_doublet']

# Stragglers

deconvolved <- rbindlist(list(good_triplets, twothirds_triplets, double_pools))

deconvolved.long <- melt(deconvolved, id='donor', measure.vars=c('SUBLIB1','SUBLIB2','SUBLIB3'), variable.name='SUBLIB',value.name='NAME')[!is.na(NAME)]
deconvolved.long[, LIB := gsub('\\.[012]','', NAME), by=.I]
pools <- fread('DATA/syn51753257/AMP-AD_DiverseCohorts_assay_multiome_metadata.csv')
pools <- pools[, .SD, .SDcols=c('specimenID','libraryBatch')]
pools[, c('donor','region') := tstrsplit(specimenID, split='_')]
setnames(pools, 'libraryBatch', 'LIB')
pools[, 'specimenID' := NULL]

final <- merge(deconvolved.long, pools, by=c('donor','LIB'))
final <- final[, .SD, .SDcols=c('LIB','NAME','donor','region')]
setnames(final, 'NAME','SUBLIB')
setnames(final, 'donor','IndividualID')
setnames(final, 'region','REGION')
fwrite(final, file='DATA/FREEMUXLET/FINGERPRINTS/deconcolved.tsv', quote=F)