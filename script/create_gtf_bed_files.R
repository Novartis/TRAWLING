args <- commandArgs(trailingOnly=TRUE)
library(Rgb)
library(parallel)
library(dplyr)
input <- args[1]
outputD <- args[2]
n <- as.numeric(args[3])
typeS <- args[4]
print(args)

canonical_transcript = read.gtf(input)
x <- canonical_transcript[which(canonical_transcript$feature %in% c('transcript','exon')),]
string<- unique(na.omit(x$transcript_id))
x_split <- split(x,x$transcript_id)

dir.create(file.path(outputD, '/gtf_tmp'))
dir.create(file.path(outputD, '/bed_tmp'))

create_bed <- function(transc,strand_ext){
  #ID site = id jnt + D (Donor Site) or + A (Acceptor Site)
  bed_f <- data.frame(chr = as.character(transc[1]), start = 0 , end = 0, id = as.character(transc[17]), 
                      val = 0, strand = strand_ext, stringsAsFactors = F)
  #2 lines (donor site + acceptor site)
  bed_f <- bed_f %>% slice(rep(1:n(), each = 2))
  # donor site: 4 + GT + 8 nucl 
  # acceptor site: 8 nucl + AG + 4
  #Add D and A to the jnt ID
  bed_f[1,2:3] <- c(as.numeric(transc[27])-11, as.numeric(transc[27])+3) 
  bed_f[2,2:3] <- c(as.numeric(transc[28])-4, as.numeric(transc[28])+10)
  
  if(strand_ext=='+'){
    bed_f[1,4] <- paste0(bed_f$id[1],'.A') 
    bed_f[2,4] <- paste0(bed_f$id[2],'.D')
    
  } else if(strand_ext=='-'){
    bed_f[1,4] <- paste0(bed_f$id[1],'.D') 
    bed_f[2,4] <- paste0(bed_f$id[2],'.A')
  }
  return(bed_f)
}


create_bed_jnc <- function(jnc,strand_ext){
  #ID site = id jnt + D (Donor Site) or + A (Acceptor Site)
  bed_f <- data.frame(chr = jnc[1], start = 0 , end = 0, transcript_name = jnc[9], 
                      val = 0, strand = strand_ext, stringsAsFactors = F)
  #2 lines (donor site + acceptor site)
  bed_f <- bed_f %>% slice(rep(1:n(), each = 2))
  # donor site: 4 + GT + 8 nucl 
  # acceptor site: 8 nucl + AG + 4
  bed_f[1,2:3] <- c(as.numeric(jnc[12])-5, as.numeric(jnc[12])+9) 
  bed_f[2,2:3] <- c(as.numeric(jnc[13])-10, as.numeric(jnc[13])+4)
  #Add D and A to the jnt ID
  if(strand_ext=='+'){
    bed_f[1,4] <- paste0(bed_f$transcript_name[1],'.D') 
    bed_f[2,4] <- paste0(bed_f$transcript_name[2],'.A')
  } else if(strand_ext=='-'){
    bed_f[1,4] <- paste0(bed_f$transcript_name[1],'.A') 
    bed_f[2,4] <- paste0(bed_f$transcript_name[2],'.D')
  }
  return(bed_f)
}
options(scipen = 999)
create_new_gtf <- function(y,n, type){
  y_ex <- y[y$feature=='exon',]
  strand <- as.character(unique(y$strand))
  y_ex <- y_ex[order(y_ex$exon_number, decreasing = T),]
  y_ex$size <- (y_ex$end-y_ex$start)+1
  cumV <-cumsum(y_ex$size)
  y_ex$oldS <- y_ex$start
  y_ex$oldE <- y_ex$end
  if(type=='3prim'){
    idx <- which(cumV>=400)[1]
    if(!is.na(idx)){
      y_ex <- y_ex[1:idx,]
      if(strand=='-'){
        y_ex$end[idx] <- y_ex$end[idx]-cumV[idx]+400
      } else{
        y_ex$start[idx] <- y_ex$start[idx]+cumV[idx]-400
      }
    }
  } else if(type=='bulk'){
    y_ex = y_ex
  } 
  y_ex$transcript_name <- gsub('\\.','_',y_ex$transcript_name)
  y_ex$gene_name <- gsub('\\.','_',y_ex$gene_name)
  y_ex$transcript_nameTrue <- y_ex$transcript_name
  y_ex$transcript_name  <-paste0(y_ex$seqname,':',y_ex$transcript_name,'_', seq(1,length(y_ex$transcript_name)))
  y_ex$transcript_name <- gsub("[.]\\d+(?=:)", "", y_ex$transcript_name, perl=TRUE)
  tmp <- y_ex[,c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand','frame','transcript_name','gene_name', 'exon_number','size','transcript_nameTrue')]
  tmp[is.na(tmp)] <- '.'
  y_ex2 <- tmp[,-c(11)]
  tmp$sizeN <- (tmp$end-tmp$start)+1
  ##############################
  #####   Create gtf cdna  #####
  ##############################
  y_cdna <- tmp
  y_cdna <- y_cdna[order(y_cdna$exon_number),]
  y_cdna <- rbind(y_cdna[1,], y_cdna)
  y_cdna$source <-'cdna'
  y_cdna[1,'feature'] <- 'transcript'
  if(strand=='+'){
    y_cdna[1,4] <- y_cdna[2,4] 
    y_cdna[1,5] <- y_cdna[dim(y_cdna)[1],5]
  } else if(strand=='-'){
    y_cdna[1,5] <- y_cdna[2,5] 
    y_cdna[1,4] <- y_cdna[dim(y_cdna)[1],4]
  }
  y_cdna$transcript_nameTrue  <-paste0(y_cdna$seqname,':','f_',y_cdna$transcript_nameTrue)
  y_cdna$transcript_nameTrue  <- gsub("[.]\\d+(?=:)", "", y_cdna$transcript_nameTrue, perl=TRUE)
  string2 <- paste0('transcript_id',' \"',y_cdna$transcript_nameTrue,'\";',' gene_id',' "',y_cdna$gene_name,'";')
  y_cdna <- y_cdna[,-c(9:13)]
  y_cdna$string2 <- string2
  
  ##############################
  ##### extract junctions ######
  ##############################
  #remove junctions if length 
  #exon <n
  
  ## if more than 1 exon
  if(dim(y_ex2)[1]>1){
    list_jnc <- vector('list',dim(tmp)[1]-1 )
    if(strand=='+'){
      tmp <- tmp[order(tmp$exon_number),]
    }
    for(jnD in 1:(dim(tmp)[1]-1)){
      
      if(tmp[jnD,'sizeN']>=n & tmp[jnD+1,'sizeN']>=n){
        df_jnc <-  tmp[1,]
        df_jnc<-do.call("rbind", replicate(3, df_jnc, simplify = FALSE))
        df_jnc$source <-'annotated_1'
        tmpSTART1 <-  tmp[jnD,'end']-n+1
        tmpEND1 <-  tmp[jnD,'end']
        tmpSTART2 <-  tmp[jnD+1,'start']
        tmpEND2 <-  tmp[jnD+1,'start']+n-1
        df_jnc[1, c('feature','start','end','transcript_name')] <- c('transcript',tmpSTART1, tmpEND2,paste0(df_jnc$seqname[1],':',tmpEND1+1,'-',tmpSTART2-1,'_',df_jnc$gene_name[1]))
        df_jnc[2, c('start','end','transcript_name')] <- c(tmpSTART1, tmpEND1,paste0(df_jnc$seqname[1],':',tmpEND1+1,'-',tmpSTART2-1,'_',df_jnc$gene_name[1]))
        df_jnc[3, c('start','end','transcript_name')] <- c(tmpSTART2, tmpEND2,paste0(df_jnc$seqname[1],':',tmpEND1+1,'-',tmpSTART2-1,'_',df_jnc$gene_name[1]))
        df_jnc$transcript_name <- gsub("[.]\\d+(?=:)", "", df_jnc$transcript_name, perl=TRUE)
        df_jnc$start_int <- tmpEND1+1
        df_jnc$end_int <- tmpSTART2-1
        string_jnc <- paste0('transcript_id',' \"',df_jnc$transcript_name,'\";',' gene_id',' "',df_jnc$gene_name,'";')
        df_jnc$string1 <- string_jnc
        list_jnc[[jnD]] <- df_jnc
      } else{
        list_jnc[[jnD]] <- NA
      }
    }
    
    list_jnc <- list_jnc[!is.na(list_jnc)]
    if(length(list_jnc)>0){
      df_jncT <- dplyr::bind_rows(list_jnc)
      df_jncT_ex <- df_jncT[df_jncT$feature=='transcript',]
      df_jncT_ex <- df_jncT_ex[,-c(12,13,14)]
      #extract seq for donor/acceptor motif
      df_jncT_bed <- apply(df_jncT_ex,1, function(j) create_bed_jnc(j,strand))
      df_jncT_bed<- do.call('rbind', df_jncT_bed)
    } else {
      df_jncT <- NA
      df_jncT_bed <- NA
    }
  } else{
    df_jncT <- NA
    df_jncT_bed <- NA
  }
  
  #remove exon with size < n and create gtf file
  y_ex2 <- y_ex2[y_ex2$size>=n, ]
  y_ex2$sizeN <-  (y_ex2$end-y_ex2$start)+1
  y_ex2 <- y_ex2[y_ex2$sizeN>=n, ]
  if(dim(y_ex2)[1]>0){
    y_ex2 <- do.call("rbind", replicate(2, y_ex2, simplify = FALSE))
    y_ex2[1:(dim(y_ex2)[1]/2),'feature'] <- 'transcript'
    y_ex2 <- y_ex2[order(y_ex2$transcript_name),]
    string1 <- paste0('transcript_id',' \"',y_ex2$transcript_name,'\";',' gene_id',' "',y_ex2$gene_name,'";')
    y_ex2 <- y_ex2[,-c(9:13)]
    y_ex2$string1 <- string1
    y_ex2$source <- 'split_exon'
  } else {
    y_ex2 <- NA
  }
  
  #create gtf files if size >= n
  y_ex <- y_ex[y_ex$size >= n, ]
  y_ex$sizeN <-  (y_ex$end-y_ex$start)+1
  y_ex <- y_ex[y_ex$sizeN>=n, ]
  if(dim(y_ex)[1]>0){
    y_ex_bed <- apply(y_ex,1, function(j) create_bed(j,strand))
    y_ex_bed<- do.call('rbind', y_ex_bed)
  } else {
    y_ex_bed <- NA
  }
  
  list_res <- list(GTF_ex = y_ex2, bed_ex = y_ex_bed, 
                   GTF_jnc = df_jncT, bed_jnc = df_jncT_bed, 
                   GTFcdna = y_cdna)
  return(list_res)
}


output_res <-  mclapply(x_split, function(y) create_new_gtf(y,n, type = typeS), mc.cores = 16)


#############
### exons ###
#############
output_res_GTF_exon <- lapply(output_res, function(x) x[[1]])
output_res_GTF_exon <- output_res_GTF_exon[!is.na(output_res_GTF_exon)]
output_res_GTF_exon <- dplyr::bind_rows(output_res_GTF_exon)
#gtf

write.table(output_res_GTF_exon, file= paste0(outputD,'/gtf_tmp/homo_sapiens_gencode_exons.gtf'), 
            col.names = F, row.names = F, quote = F, sep='\t')
#bed
output_res_bed_exon <- lapply(output_res, function(x) x[[2]])
output_res_bed_exon <- output_res_bed_exon[!is.na(output_res_bed_exon)]
output_res_bed_exon <- dplyr::bind_rows(output_res_bed_exon)

write.table(output_res_bed_exon, file= paste0(outputD,'/bed_tmp/homo_sapiens_gencode_exons.bed'), 
            col.names = F, row.names = F, quote = F, sep='\t')

############
### junc ###
############
#gtf
output_res_GTF_jnc <- lapply(output_res, function(x) x[[3]])
output_res_GTF_jnc <- output_res_GTF_jnc[!is.na(output_res_GTF_jnc)]
output_res_GTF_jnc <- dplyr::bind_rows(output_res_GTF_jnc) 
output_res_GTF_jnc <- output_res_GTF_jnc[!duplicated(output_res_GTF_jnc[c('transcript_name','feature','start','end')]),]
output_res_GTF_jnc <- output_res_GTF_jnc[,-which(colnames(output_res_GTF_jnc) %in% c('transcript_name', 'start_int', 'end_int','exon_number'))]
output_res_GTF_jnc <- output_res_GTF_jnc[,-c(9:12)]

write.table(output_res_GTF_jnc, file= paste0(outputD,'/gtf_tmp/homo_sapiens_gencode_jnc.gtf'), 
            col.names = F, row.names = F, quote = F, sep='\t')
#bed
output_res_bed_jnc <- lapply(output_res, function(x) x[[4]])
output_res_bed_jnc <- output_res_bed_jnc[!is.na(output_res_bed_jnc)]
output_res_bed_jnc <- dplyr::bind_rows(output_res_bed_jnc)
output_res_bed_jnc <- output_res_bed_jnc[!duplicated(output_res_bed_jnc[c('start','end','transcript_name')]),]

write.table(output_res_bed_jnc, file= paste0(outputD,'/bed_tmp/homo_sapiens_gencode_jnc.bed'), 
            col.names = F, row.names = F, quote = F, sep='\t')


################
####  cdna  ####
################

#gtf
output_res_GTFcdna <- lapply(output_res, function(x) x[[5]])
output_res_GTFcdna <- output_res_GTFcdna[!is.na(output_res_GTFcdna)]
output_res_GTFcdna <- dplyr::bind_rows(output_res_GTFcdna)
output_res_GTFcdna <- output_res_GTFcdna[,-9]
write.table(output_res_GTFcdna, file= paste0(outputD,'/gtf_tmp/HS_cdna_gencode.gtf'), 
            col.names = F, row.names = F, quote = F, sep='\t')