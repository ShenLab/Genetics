Analysis <- function(filename,flag,annofilename,folderadd){
	
	# first step to get the annotation information 
	gene.anno.Table <- genelist_anno(filename,flag, annofilename,folderadd)
	write.table(gene.anno.Table,file=paste(filename,"_anno",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
		
}


genelist_anno <- function(filename,flag, annofilename,folderadd){
	
	# filename for the muatation and gene lists
	gene.Tab <- read.delim(filename,sep="\t")
	genelist <- as.matrix(gene.Tab[,"Gene"])
	geneTable19 <- gene_anno_table(genelist,folderadd)
	colnames(geneTable19) <-c("GENE","DESCRIPTION","ENTREZ_SUMMARY","UNIPROTKB_SUMMARY","HGMD_ENTRY","PMIDs","DISEASE","REFERENCE","MUTATION_DESCRIPTION","MGI_HUMAN_DISEASE","MGI_MOUSE_MUTATION","RVIS_TOLERANCE_SCORE","TOLERANCE_ALL_Daly","TOLERANCE_syn_Daly","TOLERANCE_mis_Daly","TOLERANCE_non_Daly","TOLERANCE_splice_Daly","TOLERANCE_frameshift_Daly","BETWEENNESS")
	
	geneTable <- cbind(gene.Tab,geneTable19[,2:11,drop=FALSE])
	
	# to get the annovar annotation information
	anno.Tab <- gene_anno_annovar(genelist,annofilename)
	colnames(anno.Tab) <- c('GENE','esp6500_all','cosmic_70','clinvar_2014')
	
	geneTable <- cbind(geneTable,anno.Tab[,c(4,3),drop=FALSE])
	geneTable <- cbind(geneTable,geneTable19[,12:19,drop=FALSE])
	
		
	if(flag!=0){	
		addinfo <- add_anno_info(genelist,flag,folderadd)
		geneTable <- cbind(geneTable,addinfo[,-1,drop=FALSE])
	}

	
	geneTable

}

gene_anno_annovar <- function(genelist,annofilename){
	
	anno.Tab  <- as.matrix(read.delim(annofilename,sep=","))
	tmptab <- matrix(".",length(genelist),4)
	tmptab[,1] <- genelist
	if(dim(anno.Tab)[1] > 0 ){
	otherinfo <- anno.Tab[,"Otherinfo"]
	anno.Tab <- anno.Tab[,1:8,drop=FALSE]
	anno.Tab <- as.matrix(cbind(anno.Tab,"."))
	for(i in 1:dim(anno.Tab)[1]){
		anno.Tab[i,9] <- unlist(strsplit(otherinfo[i],"\t"))[1] 
	}
	subs <- match(genelist,anno.Tab[,9],nomatch=-1)
	tmptab[subs > -1 ,2:4] <- as.matrix(anno.Tab[subs[subs > -1],6:8])
	}
	
	tmptab
	
}

gene_anno_table <- function(genelist,folderadd){
	
  tmpfile <- paste(folderadd,"/gene_annotation.txt",sep="")
	all.table <- read.delim(tmpfile,sep="\t",quote=NULL)
	subs <- match(genelist,all.table[,1],nomatch=-1)
	tmptab <- matrix(".",length(genelist),dim(all.table)[2])
	tmptab[subs > -1 ,] <- as.matrix(all.table[subs[subs > -1],])
	
	tmptab
}

add_anno_info <- function(genelist,flag,folderadd){

	addinfo <- matrix(".",length(genelist),5)
	colnames(addinfo) <- c("Gene","HHE","BCells","OBESITY","brain.FRKM")
	addinfo[,1] <- genelist	
	
	# HHE gene information: now is Y AND N :cardiomyopathy
	#tmp <- as.matrix(read.table("./rawdata/HHE_mouse_E14.5_top-quartile_genes.txt"))
	tmp <- as.matrix(read.table(paste(folderadd,"/HHEtop25.txt",sep="")))
	subs <- match(genelist,tmp[,1],nomatch=-1)
	addinfo[subs> -1,2] <- tmp[subs[subs > -1],2]

	# add the primary beta cells: CM cardiomyopathy
	tmp <- as.matrix(read.csv(paste(folderadd,"/BCells_Expression_human.csv",sep="")))
	subs <- match(genelist,tmp[,1],nomatch=-1)
	addinfo[subs> -1,3] <- tmp[subs[subs > -1],7]

	# add the Priority information and brain Hypothalamus gene Ranks: OBESITY	 		
	tmp <- as.matrix(read.csv(paste(folderadd,"/ObesityGeneList.csv",sep="")))
	subs <- match(genelist,tmp[,5],nomatch=-1)
	addinfo[subs> -1,4] <- tmp[subs[subs > -1],1]	
	
	tmp <- as.matrix(read.csv(paste(folderadd,"/brain.FPKM.csv",sep="")))
	subs <- match(genelist,tmp[,2],nomatch=-1)
	addinfo[subs> -1,5] <- tmp[subs[subs > -1],3]		

	
	flag <- flag + 1
	
	addinfo
}

Analysisbatch <- function(path,flag,folder,folderadd){
	
	pathfiles <- list.files(path = path, pattern=".tsv$", full.names = TRUE, include.dirs=TRUE, recursive = TRUE)

	for(i in 1:length(pathfiles)){
  		annofilename <- paste(folder,"/",basename(pathfiles[i]),".hg19_multianno.csv",sep="")
  		Analysis(pathfiles[i],flag,annofilename,folderadd)
  		print(i)
	}
	
}

Analysistest <- function(){
  
  source("Analysis.R")
  path <- "/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/LijiangPAH/filtered"
  flag <- 1
  folder <- "./annovar/PAH"
  folderadd <- "./rawdata"
  Analysisbatch(path,flag,folder,folderadd)
  
  path <- "/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/EmilyDiabetes/filtered"
  flag <- 1
  folder <- "./annovar/DIA"
  folderadd <- "./rawdata"
  Analysisbatch(path,flag,folder,folderadd)
  
  path <- "/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/RichardObesity/filtered"
  flag <- 1
  folder <- "./annovar/OBE"
  folderadd <- "./rawdata"
  Analysisbatch(path,flag,folder,folderadd)
  
  # Note: # there are two files with name "Fam231.denovo.tsv" under this folder. 
  
  
  path <- "/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/TeresaCM/filtered"
  flag <- 1
  folder <- "./annovar/CM"
  folderadd <- "./annovar/CMtsv"
  Analysisbatch(path,flag,folder,folderadd)
  
  
  path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/EmilyDiabetes/EmilyExtraFam"
  flag <- 1
  folder <- "./annovar/DIAadd"
  folderadd <- "./rawdata"
  Analysisbatch(path,flag,folder,folderadd)  
  
  
  path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/JuliaMendelian/filtered"
  flag <- 1
  folder <- "./annovar/MEN"
  folderadd <- "./rawdata"
  Analysisbatch(path,flag,folder,folderadd)	
  
}

Readme <- function(){
  
  ##=================================================================
  strreadme <- "This file is mainly for annotation information for the tsv files under a given directory.
  The mainly function is Analysisbatch with four parameters.
  path   -- The tsv file path.
  flag   -- A reserved parameter for additional information.
  folder   -- The annovar annotation files' directory.
  folderadd   -- The added information directory for different projects. The added information contains HHE gene, BCells_Expression_human, ObesityGenelist, brain.hypothalamus.FPKM.ranking. And also, all the annotation informations.

  Usage:
  Analysisbatch(path,flag,folder,folderadd)
  

  Examples:
  source(\"Analysis.R\")
  path <- \"/ifs/scratch/c2b2/ys_lab/yshen/WENDY/Regeneron/Filtering/LijiangPAH/filtered\"
  flag <- 1
  folder <- \"./annovar/PAH\"
  folderadd <- \"./rawdata\"
  Analysisbatch(path,flag,folder,folderadd)

  The function Analysistest is a test function for different project annotations.
  They also can be seen as several examples.


  Notes:
  The necessary files for this code is listed as follows:
  HHEtop25.txt
  BCells_Expression_human.csv
  ObesityGeneList.csv
  brain.FPKM.csv
  gene_annotation.txt
  They are all in the below directory and anyone can read and write them.
  /ifs/scratch/c2b2/ys_lab/qh2159/VCFannotation/rawdata

  For now, the annovar annotated files contain three information columns as: 
  esp6500_all 
  cosmic_70 
  clinvar_2014 
  If there is any change for the columns, the code will be changed which should be updated later.
  "
  
  cat(strreadme)
  ##=================================================================
}








