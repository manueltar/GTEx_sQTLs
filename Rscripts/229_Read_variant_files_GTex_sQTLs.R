
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("liftOver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("R.utils", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BSgenome.Hsapiens.NCBI.GRCh38", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

data_wrangling = function(option_list)
{

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")

  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform Tissue_oI ----
  
  Tissue_oI = opt$Tissue_oI
  
  cat("Tissue_oI_\n")
  cat(sprintf(as.character(Tissue_oI)))
  cat("\n")
  
  Tissue_oI_adapted<-gsub("_"," ",Tissue_oI)
  
  
  cat("Tissue_oI_adapted_0\n")
  cat(sprintf(as.character(Tissue_oI_adapted)))
  cat("\n")
  
  #### READ and transform sQTL_genes_array ----
  
  sQTL_genes_array = unlist(strsplit(opt$sQTL_genes_array, split=','))
  
  cat("sQTL_genes_array_0\n")
  cat(sprintf(as.character(sQTL_genes_array)))
  cat("\n")
  
  
  #### LOOP to read input files ----
  
  Results<-data.frame()
  
  
  DEBUG<-0
  
  for(i in 1:length(sQTL_genes_array)){
    
    sQTL_genes_array_sel<-sQTL_genes_array[i]
    
    cat("---------------->\t")
    cat(sprintf(as.character(sQTL_genes_array_sel)))
    cat("\n")
    
    #### READ input_file ----
    
    filename<-paste('/group/soranzo/manuel.tardaguila/GTEX_analysis/',sQTL_genes_array_sel,'_GTEx_Portal.csv',sep='')
    
    input_file<-as.data.frame(fread(file=filename, sep=",", header=T), stringsAsFactors=F)
    
    # Gencode Id,Gene Symbol,Variant Id,SNP Id,Phenotype Id,Intron Id,P-Value,NES,Tissue
    # ENSG00000205639.9,MFSD2B,chr2_24018699_A_T_b38,rs114340515,chr2:24017588:24021648:clu_37714:ENSG00000205639.9,24017588:24021648:clu_37714,2.0e-17,0.83,Whole Blood
    # 
    
    colnames(input_file)[which(colnames(input_file) == 'Variant Id')]<-'VAR'
    colnames(input_file)[which(colnames(input_file) == 'SNP Id')]<-'rsID'
    
    
    if(DEBUG == 1){
      
      cat("input_file_0\n")
      cat(str(input_file))
      cat("\n")
      cat(str(unique(input_file$rsID)))
      cat("\n")
      cat(str(unique(input_file$VAR)))
      cat("\n")
    }
    
    
    input_file$Chr<-gsub("_.+$","",input_file$VAR)
    input_file$Pos_GRCh38<-gsub("^[^_]+_","",input_file$VAR)
    input_file$Pos_GRCh38<-gsub("_.+$","",input_file$Pos_GRCh38)
    input_file$Ref<-gsub("^[^_]+_[^_]+_","",input_file$VAR)
    input_file$Ref<-gsub("_.+$","",input_file$Ref)
    input_file$Alt<-gsub("^[^_]+_[^_]+_[^_]+_","",input_file$VAR)
    input_file$Alt<-gsub("_.+$","",input_file$Alt)
    
    if(DEBUG == 1){
      
      cat("input_file_1\n")
      cat(str(input_file))
      cat("\n")
      cat(str(unique(input_file$rsID)))
      cat("\n")
      cat(str(unique(input_file$VAR)))
      cat("\n")
    }
    
    
    
    input_file_sel<-input_file[which(input_file$Tissue == Tissue_oI_adapted),]
    
    
    if(DEBUG == 1){
      
      cat("input_file_sel_0\n")
      cat(str(input_file_sel))
      cat("\n")
      cat(str(unique(input_file_sel$rsID)))
      cat("\n")
      cat(str(unique(input_file_sel$VAR)))
      cat("\n")
    }
  
    
    if(dim(input_file_sel)[1] >0){
      
      cat(sprintf(as.character("Hello_world")))
      cat("\n")
      
      Results<-rbind(input_file_sel,Results)
      
    }# dim(input_file_sel)[1] >0
    
  }#i in 1:length(sQTL_genes_array)
    
  
if(dim(Results)[1] >0){
    
    indx.int<-c(which(colnames(Results) == 'Chr'),which(colnames(Results) == 'rsID'),which(colnames(Results) == 'Pos_GRCh38'),which(colnames(Results) == 'Ref'),
                which(colnames(Results) == 'Alt'))
    
    
    Results_subset<-unique(Results[,indx.int])
    
    cat("Results_subset_0\n")
    cat(str(Results_subset))
    cat("\n")
    
    gr_VARS <- GRanges(
      seqnames = as.character(gsub("^chr","",Results_subset$Chr)),
      ranges=IRanges(
        start=as.numeric(Results_subset$Pos_GRCh38),
        end=as.numeric(Results_subset$Pos_GRCh38),
        name=Results_subset$rsID))
    
    
    df<-as.data.frame(getSeq(BSgenome.Hsapiens.NCBI.GRCh38, gr_VARS))
    colnames(df)<-"ref"
    
    cat("df_0\n")
    cat(str(df))
    cat("\n")
    
    
    Results_subset$ref<-df$ref
    
    cat("Results_subset_1\n")
    cat(str(Results_subset))
    cat("\n")
    
    indx.int<-which(Results_subset$ref != Results_subset$Ref)
    
    cat("indx.int_0\n")
    cat(str(indx.int))
    cat("\n")
    
    Results_subset$alt<-Results_subset$Alt
    
    if(length(indx.int) >0)
    {
      Results_subset$alt[indx.int]<-Results_subset$Ref[indx.int]
      
    }else{
      
      
    }#length(indx.int) >0
    
    cat("Results_subset_2\n")
    cat(str(Results_subset))
    cat("\n")
    
    
    
    ####vcf ----
    
    indx.int<-c(which(colnames(Results_subset) == 'Chr'),which(colnames(Results_subset) == 'rsID'),which(colnames(Results_subset) == 'Pos_GRCh38'),
                which(colnames(Results_subset) == 'Ref'),which(colnames(Results_subset) == 'Alt'))
    
    cat("indx.int_0\n")
    cat(str(indx.int))
    cat("\n")
    
    vcf_df<-unique(Results_subset[,indx.int])
    
    cat("vcf_df_0\n")
    cat(str(vcf_df))
    cat("\n")
    
    Results<-merge(Results,
                          vcf_df,
                          by=c('Chr','rsID','Pos_GRCh38','Ref','Alt'),
                          all.x=T)
    
    
    
    cat("Results_1\n")
    cat(str(Results))
    cat("\n")
    
    
    
    indx.int<-which(Results$ref != Results$Ref)
    
    cat("indx.int_0\n")
    cat(str(indx.int))
    cat("\n")
    
    
    if(length(indx.int) >0)
    {
      Results$b[indx.int]<--1*(Results$b[indx.int])
      Results$bJ[indx.int]<--1*(Results$bJ[indx.int])
      
    }else{
      
      
    }#length(indx.int) >0
    
    cat("Results_2\n")
    cat(str(Results))
    cat("\n")
    
    
    vcf_df$QUAL<-100
    vcf_df$FILTER<-'PASS'
    vcf_df$INFO<-NA
    vcf_df$FORMAT<-NA
    
    colnames(vcf_df)[which(colnames(vcf_df) == 'rsID')]<-'rs'
    colnames(vcf_df)[which(colnames(vcf_df) == 'Pos_GRCh38')]<-'pos'
    colnames(vcf_df)[which(colnames(vcf_df) == 'Chr')]<-'chr'
    colnames(vcf_df)[which(colnames(vcf_df) == 'Ref')]<-'ref'
    colnames(vcf_df)[which(colnames(vcf_df) == 'Alt')]<-'alt'
    
    
    cat("vcf_df_1\n")
    cat(str(vcf_df))
    cat("\n")
    
    
    
    vcf_df$chr<-paste('chr',vcf_df$chr,sep='')
    
    vcf_df$VAR<-paste(vcf_df$chr,vcf_df$pos,vcf_df$ref,vcf_df$alt, sep="_")
    
    vcf_df$chr<-as.character(gsub("chr","",vcf_df$chr))
    
    cat("vcf_df_2\n")
    cat(str(vcf_df))
    cat("\n")
    
    #### matrix GRCh38 ----
    
    indx.int<-c(which(colnames(vcf_df) == 'chr'),which(colnames(vcf_df) == 'pos'),which(colnames(vcf_df) == 'rs'),which(colnames(vcf_df) == 'ref'),which(colnames(vcf_df) == 'alt'),
                which(colnames(vcf_df) == 'QUAL'),which(colnames(vcf_df) == 'FILTER'),which(colnames(vcf_df) == 'INFO'),which(colnames(vcf_df) == 'FORMAT'))
    
    
    vcf_df_sel_subset<-unique(vcf_df[,indx.int])
    
    cat("vcf_df_sel_subset_0\n")
    cat(str(vcf_df_sel_subset))
    cat("\n")
    
    vcf_df_sel_subset$chr<-factor(vcf_df_sel_subset$chr,
                                  levels=c("1","2","3","4","5","6","7","8","9","10","11",
                                           "12","13","14","15","16","17","18","19","20","21",
                                           "22","23","X","Y"), ordered=T)
    
    cat("vcf_df_sel_subset_0\n")
    cat(str(vcf_df_sel_subset))
    cat("\n")
    cat(str(unique(vcf_df_sel_subset$rs)))
    cat("\n")
    
    Results_vcf_df_sel_subset<-vcf_df_sel_subset[order(vcf_df_sel_subset$chr, vcf_df_sel_subset$pos),]
    
    
    
   
    #### rbind vcf files ----
    
    DEF<-unique(Results_vcf_df_sel_subset)
    DEF<-DEF[order(DEF$chr, DEF$pos),]
    
    cat("DEF_0\n")
    cat(str(DEF))
    cat("\n")
    cat(str(unique(DEF$rs)))
    cat("\n")
    
    
    #### export vcf GRCh38 ----
    
    filename<-'GTEx_sQTLs_GRCh38.vcf'
    setwd(out)
    
    cat(c("##fileformat=VCFv4.0",
          "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT"),
        file=filename,sep="\n")
    write.table(DEF, 
                file=filename, append = TRUE, sep="\t", 
                quote=F,col.names = F, row.names = F, eol="\n")
    
    
    
    #### Annotation file ----
    
    colnames(Results)[which(colnames(Results) == 'rsID')]<-'rs'
    colnames(Results)[which(colnames(Results) == 'Pos_GRCh38')]<-'pos'
    colnames(Results)[which(colnames(Results) == 'Chr')]<-'chr'
    colnames(Results)[which(colnames(Results) == 'Ref')]<-'ref'
    colnames(Results)[which(colnames(Results) == 'Alt')]<-'alt'
    
    cat("Results_1\n")
    cat(str(Results))
    cat("\n")
    
    
    

    Results$VAR<-paste(Results$chr,Results$pos,Results$ref,Results$alt, sep="_")
    
    Results$chr<-as.character(gsub("chr","",Results$chr))
    
    
    Results$CLASS<-'GTEx_sQTLs'
    
    cat("Results_FINAL\n")
    cat(str(Results))
    cat("\n")
    cat(str(unique(Results[,7])))
    cat("\n")
    
    
    setwd(out)
    write.table(Results, file='GTEx_sQTLs.tsv',sep="\t",quote=F, row.names = F)
    
  }#dim(Results)[1] >0)
}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--sQTL_genes_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Tissue_oI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  data_wrangling(opt)
  
  
}


###########################################################################

system.time( main() )