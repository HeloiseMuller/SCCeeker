#!/usr/bin/env Rscript
if(!"optparse" %in% installed.packages()[, "Package"])install.packages("optparse")
library(optparse)

print("SCCeeker version v.0.1")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input file path", metavar="character"),
  make_option(c("-i", "--minID"), type="numeric", default=90,  
              help="minimum percentage of identity [default: 90]", metavar="numeric"),
  make_option(c("-c", "--minCov"), type="numeric", default=60, 
              help="minimum coverage [default: 60]", metavar="numeric"),
  make_option(c("-d", "--db_genes"), type="character", default="mec_database_20171117,single_genes_database_20171117",  
              help="database(s) names containing genes used for typing. If several, separate with a comma. [Default: mec_database_20171117,single_genes_database_20171117]", metavar="character"),
  make_option(c("-w", "--db_whole"), type="character", default=NULL,  
              help="database(s) names containing whole cassettes. If several, separate with a comma", metavar="character"),
  make_option(c("-k", "--kmer"), type="character", default=NULL,  
              help="file for the k-mer approach", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="path where to save the outputs. If none is providing, ouputs will be save where the input is located", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Charge those library only now so does not print anything about them when one just want to print the help
packages <- c("data.table", "dplyr", "stringr", "stringi", "ggplot2", "gridExtra", "ggpubr")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library(data.table)
library(dplyr)
library(stringr)
library(stringi)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#Function I need as soon as arguments
splitToColumns <- function(vect,
                           split, #pattern at which we want to split the vector
                           maxColumns, #optional, stop splitting once reach this number of columns --> the last column contains everything that left
                           columns, #optional: return just this column
                           empty = "",
                           mode = "character") {
  # split a character vector and returns a matrix. Characters are split upon each
  # occurence of the split argument (a character of length one). empty defines
  # what the empty cell should be, and mode is the mode of the returned matrix.
  # This basically does what stri_split_fixed() does, but may sometime be more
  # convenient or faster
  
  vect <- as.character(vect)
  nc <- stri_length(split)
  out <- NULL
  
  if (!missing(maxColumns)) {
    maxCols <- as.integer(maxColumns)
  } else {
    maxCols <- NA
  }
  
  if (!missing(columns)) {
    cols <- as.integer(columns)
    cols <- cols[cols >= 1]
    out <- matrix(empty, length(vect), length(cols))
  } else {
    cols <- 1:1000
  }
  
  col <- 1
  
  repeat {
    temp <- vect
    pos <- as.vector(regexpr(split, vect, fixed = T)) - 1
    f <- pos >= 0
    temp[f] <- stri_sub(vect[f], 0, pos[f])
    
    if (!missing(columns)) {
      if (col %in% cols) {
        out[, match(col, cols)] <- temp
      }
    } else {
      out <- cbind(out, temp)
    }
    
    col <- col + 1
    
    if (col > max(cols) | !any(f)) {
      break
    }
    
    vect[f] <- stri_sub(vect[f], pos[f] + nc + 1, stri_length(vect[f]))
    vect[!f] <- empty
    
    if (!is.na(maxCols)){
      if (col == maxCols){
        out <- cbind(out, vect[f])
        break
      }
    }
  }
  
  if (ncol(out) == 1) {
    out <- as.vector(out)
  }
  
  storage.mode(out) <- mode
  
  out
}

#Proceed to verification regarding arguments

if(is.null(opt$file)){
  print_help(opt_parser)
  stop("No input file was provided", call.=FALSE)
} else {
  file <- opt$file
}

if(is.null(opt$out)){
  out <-  dirname(opt$file)
} else {
  out <- opt$out
}
print(paste0("Outputs will be saved at: ", out, "/"))

print(paste0("Minimum percentage of identity set at ", opt$minID, "%"))
if(opt$minID>100 | opt$minID<0){
  stop("The value for --minID does not look like a percentage", call.=FALSE)
} else if(opt$minID>=0 & opt$minID<=1){
  print("WARNING: are you sure you provided a percenage for minID, not a proportion?")
}
minID <- opt$minID

print(paste0("Minimum coverage of identity set at ", opt$minCov, "%"))
if(opt$minCov>100 | opt$minCov<0){
  stop("The value for --minCov does not look like a percentage", call.=FALSE)
} else if(opt$minCov>=0 & opt$minCov<=1){
  print("WARNING: are you sure you provided a percenage for minCov, not a proportion?")
}
minCov <- opt$minCov

DB_genes_typing <- splitToColumns(opt$db_genes, split = ",")
print("SCCeeker will use the genes of the following database(s) for typing:")
for(i in 1:length(DB_genes_typing)){print(paste0("    ", DB_genes_typing[i]))}

if (is.null(opt$db_whole)){
  print("No verification on database of whole cassette")
  do_whole = F
} else {
  do_whole = T
  db_whole <- splitToColumns(opt$db_whole, split = ",")
  print("Verification on whole cassettes; using DB at:")
  for(i in 1:length(db_whole)){print(paste0("    ", db_whole[i]))}
}

if (is.null(opt$kmer)){
  print("k-mer approach not activated")
  do_kmer = F
} else {
  print(paste("k-mer apprach activated; using DB at", opt$kmer))
  do_kmer = T
}


inputPattern = sub(".out", "", basename(file))

dir.create(paste0(out, "/detailed_outputs"))

#File in which class and types are defined
source("~/Software/SCCeeker/typing.R")

####################################################
# STEP 1: select hits from raw output of blast
####################################################

print("Step 1...")

#Read output of SCCeeker.sh
dt <- fread(file,
            col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "slen", "evalue", "bitscore", "out", "database" ))

#Add columns of interest
dt[, gene:= splitToColumns(sseqid, split = ":", columns = 1)]
dt[, geneNumber:= splitToColumns(sseqid, split = ":", columns = 2)]

####################################################
#Select hits that pass thresholds

#Calculate cov
dt[, scov := length/slen]

dt_filtered <- dt[pident>=minID & scov>=(minCov/100)]

#Calculate same score as in bitbucket perl scripts of SCCmecFinder (master branch)
dt_filtered[, score_pl := slen - length + gapopen +1]


####################################################
## Keep best hits

#Function to identify what hits overlap
genes_overlap <- function(dt){

  setorder(dt, qstart)
  
  #If only one gene, well 1 region
  if(nrow(dt)==1){
    dt[ , "region"] <- 1
  } else {
    
    #Create a column to write the region
    dt[ , "region"] <- 0
    #Init region 
    region_i = 1
    dt[1,]$region = 1
    
    #Go through each line (each gene)
    for(i in 2:nrow(dt)){
      #If next gene overlap the previous one, they are in the same region
      if(dt[i,]$qstart >= dt[region==region_i][1,]$qstart & dt[i,]$qstart < dt[region==region_i][1,]$qend){
        dt[i,]$region = dt[(i-1),]$region
        #Otherwise, create a new region for that gene
      } else {
        region_i = region_i + 1 
        dt[i,]$region = region_i
      }
    }
  }
  return(dt)
}

#Those are all the pattern of genes used for typing
geneFamilies =  c("mecA", "mecC", "mecR1", "mecI", "IS1272", "ccrA", "ccrB", "ccrC", "subtyp", "mec-class")
  #dmecR1 is contained in pattern mecR1
  #ccrA1, ccrA2, etc are contained in pattern ccrA
#Note: I will deal with whole_DB separately. Because I cannot use 'SCCmec_type_" as pattern since also contain pattern "subtype" for some
  #Also, let the option not use that DB if deal with it separately

#Only keep hits that corresponds to genes used for typing
dt_filtered_whole <- dt_filtered[database %in% db_whole,]
dt_filtered_forTyping <- dt_filtered[database %in% DB_genes_typing,]

#Get gene family
dt_filtered_forTyping$i_fam <- as.integer(sapply(X = dt_filtered_forTyping$gene,  FUN = str_which, pattern = geneFamilies))
dt_filtered_forTyping[, fam:=geneFamilies[i_fam]]

#For each fam on each contig, give same region to genes that overlap
dt_filtered_forTyping <- rbindlist(lapply(split(dt_filtered_forTyping, paste(dt_filtered_forTyping$out, dt_filtered_forTyping$qseqid, dt_filtered_forTyping$fam)), genes_overlap))
#Give different regionID to different fam and contig
dt_filtered_forTyping[, regionID := paste(qseqid, region, fam, sep = ":")] 

#We follow the same ranking as SCCmecFinder, ie the best hit is the one with the lowest score_pl
  #In case of equality, the best hit is the one with the highest bitscore
setorder(dt_filtered_forTyping, out, score_pl, -bitscore, geneNumber)

#For each overlapping hits, we only keep the best hit per fam and contig
dt_filtered_forTyping <- dt_filtered_forTyping[!duplicated(dt_filtered_forTyping$regionID)]

#We don't need these columns anymore
dt_filtered_forTyping <- dt_filtered_forTyping[, -c("regionID", "region", "i_fam", "geneNumber")]

fwrite(dt_filtered_forTyping, paste0(out, "/detailed_outputs/", inputPattern, "_filtered.out"), sep = "\t") 

####################################################
# STEP 2: Identify complexes ()
####################################################
print("Step 2...")

################################################
# Identify regions <100kb in which genes of interest
################################################

### Make regions, defined as genes in the same area
 #By same area, I mean that 
  #1) the total length of the region (from 1st gene to last) has to be <maxLength
  #2) We can accept a longer region if the distance between the gene and the previous one is <lengthNext
      
genes_close <- function(dt, maxLength, lengthNext){
  
  setorder(dt, qstart)
  
  #If only one gene, well 1 region
  if(nrow(dt)==1){
    dt[ , "region"] <- 1
  } else {
    
    #Order coordinates: will look at 1st first
    setorder(dt, qstart)
    #Create a column to write the region
    dt[ , "region"] <- 0
    #Init region 
    region_i = 1
    dt[1,]$region = 1
    
    #Go through each line (each gene)
    for(i in 2:nrow(dt)){
      #If gene i is close (<100000 bp) of the first region of the current region,
      #it belongs to the same region
      #I also authorize a region >100000 if less than 10kb with the previous gene
      if(dt[i,]$qend-dt[region==region_i][1,]$qstart<maxLength | dt[i,]$qend-dt[(i-1),]$qstart<lengthNext ){
        dt[i,]$region = dt[(i-1),]$region
        #Otherwise, create a new region for that gene
      } else {
        region_i = region_i + 1 
        dt[i,]$region = region_i
      }
    }
  }
  return(dt)
}

#Apply the function on each assembly (and on each contig separately)
dt_selected_withRegions <- rbindlist(lapply(split(dt_filtered_forTyping, paste(dt_filtered_forTyping$out, dt_filtered_forTyping$qseqid)), genes_close, 100000, 10000))

#Get unique ID for each region
dt_selected_withRegions[, regionID := paste0(qseqid, ":", region)]

################################################
#Identify the ccrComplex & mecComplex
################################################

#For identification I remove mecA only, otherwise all my pattern will map
dt_class_sub = dt_class[comb!="mecA"]

get_preciseType <- function(dt){
  
  #Function to identify the ccr complex
  ccrComplex_f <- function(dt){
    tmp <- dt[str_detect(gene, "ccr")]
    
    if(nrow(tmp)==0){
      ccrComplex <- "none"
    } else {
      tmp[, gene:= ifelse(str_detect(gene, "ccrC1"), "ccrC1", ifelse(str_detect(gene, "ccrC2"), "ccrC2", gene))]
      ccrComplex <- gsub("ccr", "", paste(unique(tmp$gene), collapse = ""))
      #if BxAx: same direction strand as in figure
      #if AxBx: opposite direction. So becarfull, mec complex will be written with rev
    }
    return(ccrComplex)
  }
  
  #Function to identify the mec complex
  mecComplex_f <- function(dt){
    dt <- dt[!str_detect(gene, "ccr") & !str_detect(gene, "subtyp")]
    mecComplex <- paste(dt$gene, collapse = ":")
    if(mecComplex==""){
      mecComplex="none"
    }
    return(mecComplex)
  }
  
  print(unique(dt$regionID))
  
  #Put aside hits on C1 or C2 in case need this identification
  possible_mecC = dt[str_detect(gene, "mec-class")]
  
  #Put aside hits on subtypes
  possible_subtype = dt[str_detect(gene, "subtyp") ]
  if(length(unique(possible_subtype$gene))==1){
    possible_subtype = unique(possible_subtype$gene)
  } else if(nrow(possible_subtype)==0) {
    possible_subtype = ""
  } else {
    possible_subtype = "complex"
  }
  
  #Remove from tables genes we don't need for initial identification
  dt = dt[!str_detect(gene, "subtyp") &  !str_detect(gene, "mec-class")]
  #Don't need to know the allele of C1 nor C2
  dt$gene <- sub("-allele.*", "", dt$gene)
  dt$gene <- sub("_allele.*", "", dt$gene)
  #Collapse all remaining genes together
  genes = paste(dt$gene, collapse = ":" )
  
  #Is there any ccr genes at all?
  ccr_present = ifelse(!str_detect(genes, "ccr"), 0, 
                       ifelse(str_detect(genes, "ccrC") | (str_detect(genes, "ccrA") & str_detect(genes, "ccrB")), 2, 1))
  #0 means no ccr
  #1 or 2 means at least one gene ccr
  #2 is probably functional ccr (as ccrA and ccrB need each others)
  
  #Is there any mec
  mec_present  = ifelse(!str_detect(genes,"mec"), 0, 
                        ifelse(str_detect(genes, "mecA") | str_detect(genes, "mecC"), 2, 1))
  #0 means no mecA nor regulator
  #1 means at least one mec gene but no mecA
  #2 means mecA is there
  
  #Start by identifying mec complex
  candidate_mec <- dt_class_sub[str_detect(genes, dt_class_sub$comb)]
  
  if(nrow(candidate_mec)==1){ #Here we got a single candidate, easy
    #Replace the pattern of this mec complex by the name of the mec complex
    genes <- sub(candidate_mec$comb, candidate_mec$class, genes)
  } else if(nrow(candidate_mec)>1){
    for(i in 1:nrow(candidate_mec)){
      genes <- sub(candidate_mec[i,]$comb, candidate_mec[i,]$class, genes)
    }
    
  } else if(str_detect(genes, "mecA") & !str_detect(genes, "IS1272") & !str_detect(genes, "mecI") & !str_detect(genes, "mecR")) { #Here, check whether we have mecA only
    
    #If it is the case, could be mec of class C1 or C2. For that check the variable possible_mecC we put aside earlier
    getsubC_f <- function(possible_mecC){
      if(length(unique(possible_mecC$gene))==1){
        candidate_mecC = sub("mec-class-", "", unique(possible_mecC$gene)) #Get class C1 or C2
      } else if(length(unique(possible_mecC$gene))>1) {
        candidate_mecC = "C" #Check whether this ever happen
      } else if(nrow(possible_mecC)==0){ #Even though mecA without regulators, not hit on class C1 nor C2
        candidate_mecC = "mecA" 
      }
      return(candidate_mecC)
    }
    candidate_mecC <- getsubC_f(possible_mecC)
    genes <- sub("mecA", candidate_mecC, genes)
    
    #In this case, either no mec genes at all, or a combination which was never described as a class
  } else {
    genes <- genes
  }
  
  #Now, identify ccr complex
  
  #Start by easy ones: the 5 (ccC1) and the 9 (ccrC2)
  genes <- gsub("ccrC1", "5", genes)
  genes <- gsub("ccrC2", "9", genes)
  
  #Split so can go through with a loop
  out_split <- splitToColumns(genes, ":")
  
  #Init parameters for loop
  continue=T
  i=1
  while(continue==T){
    gene_i = out_split[i]
    
    #Care only if ccrA or ccrB (the only 3 to have pattern ccr at this step)
    if(str_detect(gene_i, "ccr")){
      gene_next = out_split[i+1] #Get following gene
      comb_ccr =  paste(c(gene_i, gene_next), collapse=":")
      
      #If the combination of gene_i and gene_next was described, we get the class 
      if(comb_ccr %in% dt_complex$comb){
        ccr_complex = dt_complex[comb==comb_ccr,]
        #Replace gene_i by the name of the class
        out_split[i] = ccr_complex$complex
        #Delete the next gene as already part of the class with gene_i
        if(i<length(out_split)-1){ #If other genes after gene_next, keep them
          out_split <- out_split[c(1:i,(i+2):length(out_split))]
        } else { #Otherwise, remove everything after gene_i
          out_split <- out_split[c(1:i)]
        }
      }
    }
    i = i+1 #Want to see the next gene
    
    #Once we went through all the list of genes, we stop the while loop
    if(i>=length(out_split)){
      continue=F
    }
  }
  
  #Check whether ISonly
  if(all(out_split=="IS1272")){
    final = "IS1272 only"
  } else{
    final <- paste(c(out_split), collapse = ":")
  }
  
  #Use to return only final, but more info like that
  if(nrow(dt)>0){ #dt will be empty if only gene in reagion is subtype
    final_tb <- data.table(out = unique(dt$out),
                           contig = unique(dt$qseqid),
                           regionID = unique(dt$regionID),
                           start1st_gene = min(dt$qstart),
                           endLast_gene = max(dt$qend),
                           mec_present = mec_present,
                           ccr_present = ccr_present,
                           ccrComplex = ccrComplex_f(dt), #Needed for the heatmap
                           mecComplex = mecComplex_f(dt), #Needed for the heatmap
                           genes = paste(dt$gene, collapse = ":" ),
                           preciseTyped = final,
                           subtypeJ = possible_subtype,
                           #I add coordinates cassettes without taking IS into account. Because IS can be found outside too
                           #But if IS1272 only can't do that
                           start1st_gene_noIS = ifelse(final!="IS1272 only", min(dt[gene!="IS1272"]$qstart), NA),
                           endLast_gene_noIS = ifelse(final!="IS1272 only", max(dt[gene!="IS1272"]$qend), NA))
    return(final_tb)
  }
}


# RUN the function
preciseTyped <- rbindlist(lapply(split(dt_selected_withRegions, dt_selected_withRegions$regionID), get_preciseType))

fwrite(preciseTyped, paste0(out, "/detailed_outputs/candidateRegions.out"), sep = "\t")

####################################################
# STEP 3: Typing (Beater)
####################################################
print("Step 3...")

#For each combination of typed cassette, check if in there. 
#If only one --> Gets type
#"simple" if just that pattern and nothing else
#"composite' of other things 
#If could be several types --> NT-composite
#If none --> untyped

#We need a delineation:
typing$comb_del <- sub("^", ":", typing$comb)
typing$comb_del <- sub("$", ":", typing$comb_del)
preciseTyped$preciseTyped_del <- sub("^", ":", preciseTyped$preciseTyped)
preciseTyped$preciseTyped_del <- sub("$", ":", preciseTyped$preciseTyped_del)

#What should I do with opposite?
#Also becarful if 5:A_rev, well 5:A is in the pattern

#all(str_detect(splitToColumns("2_rev:A_rev:4_rev",split = ":" ), 'rev'))

#Prepare vector containing all gene names which are alone, so cannot contain "rev"

genes_noRev <- c("5", "9", "C1", "C2", "IS1272", "mecA", "mecC", "mecRI", "mecI")
for(i in 1:10){ #Go to known to be large, not a probrem that ccrA10 does not exist (could in the future)
  genes_noRev <- append(genes_noRev, paste0("ccrA", i))
  genes_noRev <- append(genes_noRev, paste0("ccrB", i))
}

#Function to type each region
typing_f <- function(dt){
  
  #### Before typing, need to consider the pattern "_rev"
  
  #Split each gene (or comb of gene)
  test <- splitToColumns(dt$preciseTyped,split = ":" )
  
  #Write separately those which are not in the vector genes_noRev
  test_noGenes_noRev <- test[!test %in% genes_noRev]
  
  ### CASE 1: all genes of the region are genes that cannot be written with "_rev" (because single gene, not a combination)
  
  #In this case order does not matter
  #So do alphabetic order so those in different order depending on regions are in same order
  if(length(test_noGenes_noRev)==0){
    ntest <- length(test)
    if(test[1]<test[ntest]){
      precise <- paste(test, collapse = ":")
    } else{
      test_rev = c()
      for(i in ntest:1){
        test_rev  <- append(test_rev, test[i])
        }
      precise <- paste(test_rev, collapse = ":")
    }
    #We will need the delineation
    precise_del <- sub("^", ":", precise)
    precise_del <- sub("$", ":", precise_del)
    
    ### CASE 2: all gene combination have pattern "rev" (if single genes, don't have "rev" but that's ok)
    
    } else if (length(test_noGenes_noRev)>0 & all(str_detect(test_noGenes_noRev, 'rev'))){
    tmp <- sub("_rev", "", test)
    nTmp <- length(tmp)
    #Write tmp in reversed:
    tmp_rev = c()
    for(i in nTmp:1){
      tmp_rev  <- append(tmp_rev, tmp[i])
    }
    tmp_rev <- paste(tmp_rev, collapse = ":")
    
    precise <- tmp_rev
    
    #We will need the delineation
    precise_del <- sub("^", ":", tmp_rev)
    precise_del <- sub("$", ":", precise_del)
    
    ### CASE 3: nothing to do. Just keep same order of genes/combination and don't removing any "rev" pattern
  } else {
    precise = dt$preciseTyped
    precise_del = dt$preciseTyped_del
  }
  
  #We will process to typing thanks to this column

  #Go through each possible combination of typed cassettes
  out_i = c()
  for(i in 1:nrow(typing)){
    if(str_detect(precise_del, typing[i,]$comb_del)){
      out_i = append(out_i, i)
    }
  }
  #out_type <- ifelse(length(out_i)==1, typing[out_i]$type, ifelse(length(out_i)==0, "NT", "NT-composite"))
  out_type <- ifelse(length(unique(typing[out_i]$type))==1, unique(typing[out_i]$type), ifelse(length(out_i)==0, "NT", "NT-composite"))
  out_composite <- ifelse(!out_type %in% c("NT", "NT-composite") & any(typing$comb_del==precise_del), "simple", 
                          ifelse(out_type %in% c("NT", "NT-composite"), NA, "composite"))

  dt_out <- data.table(type = out_type, composite = out_composite, preciseTyped_noRev = precise)
  dt_out <- cbind(dt, dt_out)
  return(dt_out)
}

typed <- rbindlist(lapply(split(preciseTyped, preciseTyped$regionID), typing_f))

print("Summary of types recovered:")
table(typed$type)


######## Gives broader category to regions for figures
typed[, cat := ifelse(type!="NT" & type!="NT-composite", "typed",
                      ifelse(type=="NT-composite", "NT-composite",
                             ifelse(preciseTyped_noRev=="IS1272 only", "IS1272 only",
                                    ifelse(mec_present>0 & ccr_present==0, "mec but no ccr",
                                           ifelse(ccr_present>0 & mec_present==0, "ccr but no mec", 
                                                  "both but untyped")))))]
typed$cat <- factor(typed$cat, 
                    levels = c("IS1272 only", "mec but no ccr", "ccr but no mec", "both but untyped", "typed", "NT-composite"))
fwrite(typed, paste0(out, "/detailed_outputs/candidateRegions_typed.out"), sep = "\t")

###############
### FIGURES ###
###############
print("Generating some figures...")

dir.create(paste0(out, "/Figures/"))

typed$cat <- factor(typed$cat, levels = c("IS1272 only", "ccr but no mec", "mec but no ccr", "both but untyped", "NT-composite", "typed"))

pdf(paste0(out, "/Figures/regions_summary.pdf"))
typed %>%
  mutate(., catFig = ifelse((cat=="typed" & composite=="composite") | cat=="NT-composite", "composite", 
                            ifelse(cat=="typed" & composite=="simple", "simple", NA))) %>%
  group_by(cat, catFig) %>%
  summarise(nRegion=n()) %>%
  ggplot(., aes(x=cat, y=nRegion, fill=catFig)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Number of candidate regions") +
  xlab("Category of the region") +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.5, 0.85), 
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank()) 
dev.off()

## 

typed$composite <- factor(typed$composite, levels = c("simple", "composite"))
typed$type <- factor(typed$type, levels = c(typing_rev$type, "NT-composite", "NT"))
p1 <- typed[type!="NT"] %>%
  mutate(., catFig = ifelse((cat=="typed" & composite=="composite") | cat=="NT-composite", "composite", 
                            ifelse(cat=="typed" & composite=="simple", "simple", NA))) %>%
  ggplot(., aes(x=preciseTyped_noRev, fill = catFig)) +
  geom_bar() +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=8)) +
  facet_wrap(~type,scale = "free") +
  xlab("Precise combinaiton of genes") +
  ylab("Number of regions") +
  guides(fill=guide_legend(title="Cassette"))

typed$preciseTyped_noRev <- factor(typed$preciseTyped_noRev, levels = names(sort(table(typed$preciseTyped_noRev))))
p2<- ggplot(typed[type=="NT" & cat!="IS1272 only"], aes(x=preciseTyped_noRev)) +
  geom_bar() +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=8)) +
  xlab("Precise combinaiton of genes") +
  ylab("Number of regions") +
  guides(fill="none")

pdf(paste0(out, "/Figures/typing.pdf"), height = 20, width = 12)
gt <- arrangeGrob(p1,p2, layout_matrix = rbind(c(1,1,1,2,2), c(1,1,1, 2,2)))
as_ggplot(gt)
dev.off()

###############
## Check the blast on the whole cassette
###############

if(do_whole==T){
  
  dt_filtered_whole[, typeWhole:= splitToColumns(splitToColumns(sseqid, split = "|", columns = 1), split = "_", columns = 3)]
  dt_filtered_whole[, type:=sub("[(].*", "", typeWhole)]
  
  validate_wholeBlast_f <- function(dt, whole){
    if(dt$cat=="typed"){
      
      whole_tmp <- whole[qseqid==dt$contig & qstart<dt$endLast_gene & qend>dt$start1st_gene]
      
      #What is the ranking of the type identified with combination
      ranking <- which(whole_tmp$type==dt$type)[1]
      
      coord_best = whole_tmp[1, c("typeWhole", "qstart", "qend")]
      coord_match = whole_tmp[ranking, c("typeWhole", "qstart", "qend")]
      
      colnames(coord_best) = c("typeWhole_best", "startWhole_best", "endWhole_best")
      colnames(coord_match) = c("typeWhole_matchType", "startWhole_matchType", "endWhole_matchType")
      
      dt_out <- cbind(dt, coord_best, rank_matchType=ranking, coord_match)
      
    } else {
      dt_out = cbind(dt, typeWhole_best=NA, startWhole_best=NA, endWhole_best=NA, rank_matchType=NA, typeWhole_matchType=NA, startWhole_matchType=NA, endWhole_matchType=NA)
    }
    return(dt_out)
    
  }
  
  setorder(dt_filtered_whole, out, score_pl, -bitscore)
  typed_validate_wholeBlast <- rbindlist(lapply(split(typed, typed$regionID), validate_wholeBlast_f, whole = dt_filtered_whole))
  
  fwrite(typed_validate_wholeBlast, paste0(out, "/detailed_outputs/candidateRegions_typed_validateWhole.out"), sep = "\t")
  
  typed <- typed_validate_wholeBlast
  
  #Get longest coordinates between both sets
  typed[, cas:=ifelse(is.na(startWhole_best), 2, 1)]
  typed[, startTmp:=ifelse(cas==2, NA, 
                           ifelse(cas==1 & is.na(startWhole_matchType), startWhole_best,
                                  ifelse(cas==1 & startWhole_matchType<startWhole_best, startWhole_best, startWhole_best)))]
  typed[, startRegion:=ifelse(cas==2 & is.na(start1st_gene_noIS), start1st_gene, 
                              ifelse(cas==2 & ! is.na(start1st_gene_noIS),start1st_gene_noIS, 
                              ifelse(cas==1 & start1st_gene_noIS<startTmp, start1st_gene_noIS, startTmp)))]
  typed[, endTmp:=ifelse(cas==2, NA, 
                           ifelse(cas==1 & is.na(endWhole_matchType), endWhole_best,
                                  ifelse(cas==1 & endWhole_matchType>endWhole_best, endWhole_best, endWhole_best)))]
  typed[, endRegion:=ifelse(cas==2 & is.na(endLast_gene_noIS), endLast_gene,
                            ifelse(cas==2 & !is.na(endLast_gene_noIS), endLast_gene_noIS, 
                              ifelse(cas==1 & endLast_gene_noIS>endTmp, endLast_gene_noIS, endTmp)))]
  typed <- typed[, -c("cas", "startTmp", "endTmp")]

  
} else {
  typed[, startRegion:=start1st_gene_noIS]
  typed[, endRegion:=endLast_gene_noIS]
}


###############
## k-mer validation ##
##############

if(do_kmer==T){
  print("Step k-mer...")
  
  inputPattern_kmer = sub(".tsv", "", basename(opt$kmer))
  
  print(paste0("Reading k-mer output at: ", opt$kmer))
  
  #WHen my script will be clean, I can use just that
  kmer_output <- fread(opt$kmer,
                        col.names = c("Template", "Score", "Expected", "z", "p_value", "queryCoverage", "templateCoverage", "depth", "kmers_inTemplate", "Description", "out"))
  
  #Get the type
  kmer_output[, typeWhole_kmer:= splitToColumns(splitToColumns(Template, split = "|", columns = 1), split = "_", columns = 3)]
  kmer_output[, type_kmer:=sub("[(].*", "", typeWhole_kmer)]
  
  #Get the subtype
  #Get subtype thanks to k-mer:
  kmer_output[!is.na(Template), subType_kmer := sub("SCCmec_type_", "", splitToColumns(
    splitToColumns(Template, split = "|", columns = 2),
    split = "(",  columns = 1)
    )
  ]
  #If "gb", it means there are no subtype for that type
  kmer_output[subType_kmer=="gb",]$subType_kmer = NA
  
  #Rank them by best score (use z score for that)
  setorder(kmer_output, out, -z)
  
  #For each assembly, keep only one line per type
  kmer_output <- kmer_output[!duplicated(paste(out, type_kmer))]
  fwrite(kmer_output, paste0(out, "/detailed_outputs/", inputPattern_kmer, "_filtered.out"), sep = "\t")
  
  kmer_outp_best <- kmer_output[!duplicated(out)]
  kmer_outp_best$cat = "typed"

  #Add info about k-mer in our table
  typed_validate_wholeBlast_kmer <- left_join(typed, kmer_outp_best[, c("out", "typeWhole_kmer", "type_kmer", "subType_kmer" ,"cat")], by = c("out", "cat"))

  #Harmonized subtype name from combination genes with the one kmr
  typed_validate_wholeBlast_kmer$subtypeJ = sub("[(].*", "", sub("subtyppe-", "", sub("subtype-", "", typed_validate_wholeBlast_kmer$subtypeJ)))
  typed_validate_wholeBlast_kmer[subtypeJ=="" | subtypeJ=="NT-composite"]$subtypeJ = NA


  #All regions whose type is the same one as the type identified by the k-mer approach are validated
  #We also validate the types which are not in the database of kmer
  typed_validate_wholeBlast_kmer[, validated:=ifelse(type==type_kmer | type %in% c("XII", "XIV", "XV"), TRUE, FALSE)]

    #Save final table
  fwrite(typed_validate_wholeBlast_kmer, paste0(out, "/detailed_outputs/candidateRegions_typed_validatedKmer.out"), sep = "\t")
  typed <- typed_validate_wholeBlast_kmer
  
  #Save a table with positions of each gene
  genes_regions <- left_join(dt_selected_withRegions[, c("regionID", "gene", "fam","length", "qseqid",  "qstart", "qend")], 
                             typed_validate_wholeBlast_kmer[, c("out", "regionID", "cat", "type", "preciseTyped_noRev", "validated")],
                             by = "regionID")
  
} else {
  genes_regions <- left_join(dt_selected_withRegions[, c("regionID", "gene", "fam","length", "qseqid",  "qstart", "qend")], 
                             typed[, c("out", "regionID", "cat", "type", "preciseTyped_noRev")],
                             by = "regionID") 
}

fwrite(genes_regions, paste0(out, "/detailed_outputs/allGenes_withRegions.out"), sep = "\t")

#############################
### Save a summary for user we do not need downstream analyses ###
#############################

print("Generating summaries")

#A summary for each of the candidate region
if(do_kmer==T){
  fwrite(typed[, c("regionID", "out", "contig", "cat", "type", "composite", "preciseTyped_noRev", "startRegion", "endRegion", "validated")], paste0(out, "/SCCeeker_summary_perRegion.out"), sep = "\t")
} else {
  fwrite(typed[, c("regionID", "out", "contig", "cat", "type", "composite", "preciseTyped_noRev", "startRegion", "endRegion")], paste0(out, "/SCCeeker_summary_perRegion.out"), sep = "\t")
}

#A summary for each assembly
summary <- typed %>% group_by(out) %>% summarise(
  mec = ifelse(str_detect(paste(mecComplex, collapse = "-"), "mecA") & !str_detect(paste(mecComplex, collapse = "-"), "mecC"), "mecA", 
                                               ifelse(!str_detect(paste(mecComplex, collapse = "-"), "mecA") & str_detect(paste(mecComplex, collapse = "-"), "mecC"), "mecC",
                                                      ifelse(str_detect(paste(mecComplex, collapse = "-"), "mecA") & str_detect(paste(mecComplex, collapse = "-"), "mecC"), "mecA&C", "-"))),
  nCandidateRegions = n()) %>% data.table()

tmp <- typed[cat=="typed"] %>% group_by(out) %>% summarize(types = paste(type, collapse = ";")) %>% data.table() 
summary <- left_join(summary, tmp, by = "out")

if(do_kmer==T){
  tmp2 <- typed[cat=="typed" & validated==T] %>% group_by(out) %>% summarize(types_validated = paste(type, collapse = ";")) %>% data.table() 
  summary <- left_join(summary, tmp2, by = "out")
} 

 
fwrite(summary, paste0(out, "/SCCeeker_summary_perFile.out"), sep = "\t")

print(paste0("Final outputs were saved at: ", out, "/"))
print(paste0("Detailed and intermediate outputs were saved at: ", out, "/detailed_outputs/"))






