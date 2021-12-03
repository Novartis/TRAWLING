args <- commandArgs(trailingOnly=TRUE)
library(Rgb)
library(parallel)
library(dplyr)
options(scipen = 999)

inputD <- args[1]
outputD <- args[2]
typeS <- args[3]
dir.create(file.path(outputD, '/gtf'))
dir.create(file.path(outputD, '/fasta'))

bed_outputExon <- read.table(paste0(inputD,'/fasta_bed/homo_sapiens_gencode_exons2.fa'), header = F, stringsAsFactors = F)
bed_outputExon <- bed_outputExon[!duplicated(bed_outputExon[c('V1','V2')]),]

bed_outputJnc <- read.table(paste0(inputD,'/fasta_bed/homo_sapiens_gencode_jnc2.fa'), header = F, stringsAsFactors = F)
bed_outputJnc <- bed_outputJnc[!duplicated(bed_outputJnc[c('V1','V2')]),]

gtf_fileExon<- read.gtf(paste0(inputD,'/gtf_tmp/homo_sapiens_gencode_exons.gtf'))
gtf_fileJnc<- read.gtf(paste0(inputD,'/gtf_tmp/homo_sapiens_gencode_jnc.gtf'))
gtf_fileCdna<- read.gtf(paste0(inputD,'/gtf_tmp/HS_cdna_gencode.gtf'))

cat('create t2g file...\n')
t2gExon <- gtf_fileExon[!duplicated(gtf_fileExon[c('transcript_id','gene_id')]), c('transcript_id','gene_id')]
t2gExon$transcript_id <- paste0(t2gExon$transcript_id,'.s.3prim')
t2gJnc <- gtf_fileJnc[!duplicated(gtf_fileJnc[c('transcript_id','gene_id')]), c('transcript_id','gene_id')]
t2gJnc$transcript_id <- paste0(t2gJnc$transcript_id,'.s.3prim')
t2gCdna <- gtf_fileCdna[!duplicated(gtf_fileCdna[c('transcript_id','gene_id')]), c('transcript_id','gene_id')]
t2gCdna$transcript_id <- paste0(t2gCdna$transcript_id,'.s.3prim')
t2g <- rbind(t2gExon,t2gJnc,t2gCdna)
#TRAWLING
write.table(t2g, file= paste0(outputD,'/fasta/TRAWLING-cdna',typeS,'.t2g'), 
            col.names = F, row.names = F, quote = F, sep='\t')
#GENCODE
write.table(t2gCdna, file= paste0(outputD,'/fasta/HS_gencode-cdna',typeS,'.t2g'), 
            col.names = F, row.names = F, quote = F, sep='\t')

cat('add donor/acceptor sequences to gtf file...\n')

#####Exon
bed_outputExon <- split(bed_outputExon, bed_outputExon$V2)
tmp2Exon <- merge(bed_outputExon$D, bed_outputExon$A, by = 'V1')
colnames(tmp2Exon) <- c('transcript_id', 'D', 'D_motif', 'A', 'A_motif')
x <- gtf_fileExon
x$ID_rows <-  1:nrow(x)
tmp3Exon <- merge(tmp2Exon, x, by= 'transcript_id')
tmp3Exon <- tmp3Exon[order(tmp3Exon$ID_rows), ]
tmp3Exon <- tmp3Exon[,!(names(tmp3Exon)=='ID_rows')]

#create t2a file
tmp_t2aExon <-tmp3Exon[tmp3Exon$feature=='transcript',]
tmp_t2aExon$gene_id2 <- tmp_t2aExon$gene_id
tmp_t2aExon <- tmp_t2aExon[,c('transcript_id', 'source', 'gene_id','seqname','start','end','strand','D_motif','A_motif','gene_id2')]
colnames(tmp_t2aExon)[4] <- 'chr'

#####jnc
bed_outputJnc <- split(bed_outputJnc, bed_outputJnc$V2)
tmp2Jnc <- merge(bed_outputJnc$D, bed_outputJnc$A, by = 'V1')
colnames(tmp2Jnc) <- c('transcript_id', 'D', 'D_motif', 'A', 'A_motif')
x <- gtf_fileJnc
x$ID_rows <-  1:nrow(x)
tmp3Jnc <- merge(tmp2Jnc, x, by= 'transcript_id')
tmp3Jnc <- tmp3Jnc[order(tmp3Jnc$ID_rows), ]
tmp3Jnc <- tmp3Jnc[,!(names(tmp3Jnc)=='ID_rows')]

#create t2a file
tmp_t2aJnc <-tmp3Jnc[tmp3Jnc$feature=='transcript',]
tmp_t2aJnc$gene_id2 <- tmp_t2aJnc$gene_id
tmp_t2aJnc <- tmp_t2aJnc[,c('transcript_id', 'source', 'gene_id','seqname','start','end','strand','D_motif','A_motif','gene_id2')]
colnames(tmp_t2aJnc)[4] <- 'chr'


#create t2a file
tmp_t2aCdna <-gtf_fileCdna[gtf_fileCdna$feature=='transcript',]
tmp_t2aCdna$gene_id2 <- tmp_t2aCdna$gene_id
tmp_t2aCdna$D_motif <- 'NA'
tmp_t2aCdna$A_motif <- 'NA'
tmp_t2aCdna <- tmp_t2aCdna[,c('transcript_id', 'source', 'gene_id','seqname','start','end','strand','D_motif','A_motif','gene_id2')]
colnames(tmp_t2aCdna)[4] <- 'chr'


tmp_t2a <- rbind(tmp_t2aExon,tmp_t2aJnc,tmp_t2aCdna)
cat('save t2a file...')
options(scipen = 999)
write.table(tmp_t2a, file= paste0(outputD,'/gtf/TRAWLING.t2a'), 
            col.names = T, row.names = F, quote = F, sep='\t')

write.table(tmp_t2aCdna, file= paste0(outputD,'/gtf/HS_gencode.t2a'), 
            col.names = T, row.names = F, quote = F, sep='\t')


tmp3Exon$string1 <- paste0('transcript_id',' \"',paste0(tmp3Exon$transcript_id,'.s.3prim'),'\";',' gene_id',' "',tmp3Exon$gene_id,'";',
                           ' D_motif', ' "', tmp3Exon$D_motif, '";',
                           ' A_motif', ' "',tmp3Exon$A_motif,'";')
new_gtfExon <- tmp3Exon[,c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand','frame','string1')]
new_gtfExon$score <- '.'
new_gtfExon$frame <- '.'

tmp3Jnc$string1 <- paste0('transcript_id',' \"',paste0(tmp3Jnc$transcript_id,'.s.3prim'),'\";',' gene_id',' "',tmp3Jnc$gene_id,'";',
                          ' D_motif', ' "', tmp3Jnc$D_motif, '";',
                          ' A_motif', ' "',tmp3Jnc$A_motif,'";')
new_gtfJnc <- tmp3Jnc[,c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand','frame','string1')]
new_gtfJnc$score <- '.'
new_gtfJnc$frame <- '.'


gtf_fileCdna_2 <- gtf_fileCdna
gtf_fileCdna_2$D_motif <- 'NA'
gtf_fileCdna_2$A_motif <- 'NA'
gtf_fileCdna_2$string1 <- paste0('transcript_id',' \"',paste0(gtf_fileCdna_2$transcript_id,'.s.3prim'),'\";',' gene_id',' "',gtf_fileCdna_2$gene_id,'";',
                                 ' D_motif', ' "', gtf_fileCdna_2$D_motif, '";',
                                 ' A_motif', ' "',gtf_fileCdna_2$A_motif,'";')
new_gtfCdna <- gtf_fileCdna_2[,c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand','frame','string1')]
new_gtfCdna$score <- '.'
new_gtfCdna$frame <- '.'


new_gtf <- rbind(new_gtfExon,new_gtfJnc,new_gtfCdna)

cat('save gtf file...')
write.table(new_gtf, file= paste0(outputD,'/gtf/TRAWLING.gtf'), 
            col.names = F, row.names = F, quote = F, sep='\t')

write.table(new_gtfCdna, file= paste0(outputD,'/gtf/HS_gencode.gtf'), 
            col.names = F, row.names = F, quote = F, sep='\t')
