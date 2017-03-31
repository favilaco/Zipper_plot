##############################################################################################
##   Zipper plot: visualizing transcriptional activity of genomic regions                   ##
##                                                                                          ##
##   Author: Francisco Avila Cobos (Francisco.AvilaCobos@UGent.be)                          ##
##                                                                                          ##
##############################################################################################
# Script usage:
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=7){
  print("Please check that all required parameters are indicated or the paths are correct")
  print("Required parameters after Rscript (in order): Script_file, input_file; path_to_bedtools, Zipper_width, Number_permutations_for_AUZpval, mark, peak_type, output_filename")
  print("Example usage: Rscript ZP_ROADMAP_ALL.R input.txt /usr/bin/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools 5000 100 H3K4me3 narrowPeak output_1")
  stop()
} 

input_filename = args[1]
bedtools_path = args[2]
Zipper_width = as.numeric(args[3])
Nper = as.numeric(args[4])
MARK <- args[5]
TYPE <- args[6]
output_filename = args[7]

##############################################################################################

#Loading required R packages
packages <- c("data.table","ggplot2","R.utils","data.table","grid","gridExtra")
for(i in packages){require(i, character.only = TRUE)}

##############################################################################################
# Other directories, functions and parameters that are needed:

COLOUR = "darkblue"
ChromInfo = "ChromInfo.txt"
source("ZP_functions.R")
User_dir <- getwd()
Roadmap_dir <- paste(paste(paste("ACTIVATING",TYPE,sep="_"),"Roadmap",sep="_"),MARK,sep="/")
INFO <- fread("ROADMAP_info.txt")
no_gaps = paste(User_dir,"No_gaps/filtered_",sep="/")

user_params =tableGrob(rbind(c('File',output_filename),c('DataType',paste(MARK,TYPE,sep="_")),c('Tissues', 'ALL'),c('Zipper width',Zipper_width)),theme=ttheme_minimal(),rows=NULL)

##############################################################################################

setwd(User_dir)
lnc_names <- read.table(input_filename,fill=TRUE)

Input_size = countLines(input_filename)
COUNTER = 0 #To avoid dynamically growing structures

USER_INPUT <- data.frame(matrix(nrow=Input_size,ncol=3))
colnames(USER_INPUT) <- c("chr","TSS","strand")
RESULTS <- data.frame(matrix(nrow=Input_size,ncol=7))
PERMUTATIONS <- data.frame(matrix(nrow=Input_size,ncol=2*Nper))

system(paste(paste("sort -k1,2 ",input_filename,sep=""), " > Input_sorted.txt",sep=""))
system("awk \'{print >> $1; close($1)}\' Input_sorted.txt")
system("rm *Input*")


##############################################################################################

for (File in dir(User_dir,pattern="chr")){  #chr1.txt ... chrY.txt of User. CAGE split in same way!

  counter_plus = 0
  counter_minus = 0

  CHR <- File

  ### User
  User <- data.table(read.table(paste(User_dir,File,sep="/"),fill=TRUE,colClasses=c('character','numeric','character'))[,1:3])
  colnames(User) <- c("chr","TSS","strand")
  setkey(User,strand)

  if("+" %in% names(table(User[,3,with=FALSE]))){

    User_plus <- User["+"]
    setkey(User_plus,TSS)
    counter_plus = counter_plus + nrow(User_plus) 

  }

  if("-" %in% names(table(User[,3,with=FALSE]))){

    User_minus <- User["-"]
    setkey(User_minus,TSS)
    counter_minus = counter_minus + nrow(User_minus) 

  }


  ###### Keeping track of User input 
  if(exists("User_plus")){USER_INPUT[(COUNTER + 1):(COUNTER + counter_plus),] <- User_plus} 
  if(exists("User_minus")){USER_INPUT[(1 + COUNTER + counter_plus):(COUNTER + counter_plus + counter_minus),] <- User_minus} 




  ### ROADMAP
  if(length(grep("^ROADMAP$",ls()))!=0){rm(list="ROADMAP")}

  ROADMAP <- readRDS(paste(Roadmap_dir,dir(Roadmap_dir)[grep(paste(CHR,"all_outerjoin",sep=""),dir(Roadmap_dir))],sep="/"))
  ROADMAP[,width:=NULL]
  #setkey(ROADMAP,start) #Not needed anymore

  if(exists("User_plus")){RESULTS[(COUNTER + 1):(COUNTER + counter_plus),] <- closest_ROADMAP_signal(User_plus,ROADMAP,detailed=TRUE)}
  if(exists("User_minus")){RESULTS[(1 + COUNTER + counter_plus):(COUNTER + counter_plus + counter_minus),] <- closest_ROADMAP_signal(User_minus,ROADMAP,detailed=TRUE)}


  ##############################################################################################
  ##############################################################################################
  # 2) PERMUTATIONS FOR ZIPPER SHAPE PVAL:
  # Here we cannot narrow down to only 1 tissue!
  # Permutations are done across ALL tissues! This is translated into looking at ALL peaks

  RL = nrow(User)
    
  ##############################################################################################
    
  Random = paste(CHR, "random.txt",sep="")
  to_filter = paste(paste("filtered",CHR,sep="_"),".txt",sep="") #genomic regions WITHOUT gaps, heterochromatin, centromeres, telomeres and repetitive regions
  
  command1 = paste(bedtools_path, " random -g ", User_dir, "/ChromInfo.txt -l 1 -seed 4 -n ", as.integer(RL*Nper), sep="")
  command1 = paste(paste(command1, ">", sep = " "), "temp.txt" , sep=" ")
  system(command1)
  
  #No_gaps folder contained the filtered version of the genome (No gaps, centromeres, telomeres, heterochromatin, repetitive elements...)
  command2 = paste(bedtools_path," shuffle -g ", User_dir ,"/ChromInfo.txt -i temp.txt -seed 4 -chrom -incl ", paste(no_gaps,CHR,".txt",sep=""), sep="")
  command2 = paste(paste(command2, " | cut -f2 >", sep = " "), Random , sep=" ") 
  system(command2)
  
  per_tmp = fread(Random,showProgress = T)
  per_tmp = (as.vector(as.matrix(per_tmp)))
  per_tmp = matrix(per_tmp,nrow=RL,ncol=Nper) 
    
  per_tmp <- data.table(data.frame(chr = User$chr, strand = User$strand,per_tmp))
  setkey(per_tmp,strand)

  ##############################################################################################


  if("+" %in% names(table(User[,3,with=FALSE]))){

    pt_plus <- per_tmp["+"]
    colnames(pt_plus) <- c("chr","strand",rep("TSS",(ncol(pt_plus)-2)))
    Pp <- matrix(nrow=nrow(pt_plus),ncol=(ncol(per_tmp)-2)*2)

  }

  if("-" %in% names(table(User[,3,with=FALSE]))){

    pt_minus <- per_tmp["-"]
    colnames(pt_minus) <- c("chr","strand",rep("TSS",(ncol(pt_minus)-2)))
    Pm <- matrix(nrow=nrow(pt_minus),ncol=(ncol(per_tmp)-2)*2)

  }


  # This option is faster than do.call(cbind,apply...) [already checked]
  for (i in 1:(ncol(per_tmp)-2)){

    if(exists("pt_plus")){Pp[,c(2*i-1,2*i)] <- as.matrix(closest_ROADMAP_signal(pt_plus[,c(1,i+2,2),with=FALSE],ROADMAP[,1:2,with=FALSE],detailed=FALSE))}  #Only peak coordinates matter
    if(exists("pt_minus")){Pm[,c(2*i-1,2*i)] <- as.matrix(closest_ROADMAP_signal(pt_minus[,c(1,i+2,2),with=FALSE],ROADMAP[,1:2,with=FALSE],detailed=FALSE))}
    
  }

  # Needed for cases where only "+" or "-" in a given chromosome!
  if(exists("pt_plus") & exists("pt_minus")){

  	P <- rbind.data.frame(Pp,Pm)

  } else {

    if(exists("pt_plus")){

      P <- Pp
  	
    } else {

      P <- Pm
  	 
    }
  }


  colnames(P) <- rep(c("dist_start","dist_end"),(ncol(P)/2))

  #Permutations already containing distance info
  PERMUTATIONS[(COUNTER + 1): (COUNTER + counter_plus + counter_minus),] <- P

  if(exists("pt_plus")){rm(list=c("pt_plus","Pp"))}
  if(exists("pt_minus")){rm(list=c("pt_minus","Pm"))}
  	   
  	   
  ###

  if(exists("User_plus")){rm(list="User_plus")} #Needed for cases when only 1 strand features are present!!
  if(exists("User_minus")){rm(list="User_minus")}

  COUNTER =  COUNTER + counter_plus + counter_minus

  system(paste("rm","temp.txt",Random,sep=" "))
  system(paste("rm",File))
  rm(list="ROADMAP")

}


#Since some chromosomes have only features on "+" or "-" strand:
PERMUTATIONS <- data.table(na.omit(PERMUTATIONS))


RESULTS <- data.table(na.omit(RESULTS))
Z <- ORDERING_ROADMAP(RESULTS[,3:5,with=FALSE],Plot=TRUE)
ZP <- ROADMAP_ZP(Z,Zipper_width,COLOUR)


lnc_key <- do.call(paste,lnc_names[,1:3])
RESULTS_key <- do.call(paste,RESULTS[,1:3,with=FALSE])
RESULTS$name <- lnc_names[as.vector(sapply(RESULTS_key,function(x) grep(x,lnc_key))),4]
setcolorder(RESULTS,c(ncol(RESULTS),1:(ncol(RESULTS)-1)))


colnames(RESULTS) <- c("name","chr","TSS","strand","dist_ROADMAP_start","dist_ROADMAP_end","signalValue","closest_hit")

hits <- unlist(strsplit(RESULTS$closest_hit,"-"))[seq(1,(nrow(RESULTS)*2),2)]
RESULTS$Type <- INFO[match(hits,INFO[,Epigenome,]),Name,]

setkey(RESULTS,dist_ROADMAP_start)


##############################################################################################

# p-values and Summary table
# PERMUTATIONS: calculated while looping through chromosomes!
setwd(User_dir)

Rgrids <- vector()
AUZrs <- vector()
Lgrids <- vector()
AUZls <- vector()

PERMUTATIONS <- data.table(data.frame(strand=USER_INPUT$strand, PERMUTATIONS))

for (j in seq(2,ncol(PERMUTATIONS),2)){

  tmp <- ORDERING_ROADMAP(PERMUTATIONS[,c(1,j,j+1),with=FALSE])
  Rgrids <- c(Rgrids,tmp$R_grid)
  AUZrs <- c(AUZrs,tmp$AUZ_right)
  Lgrids <- c(Lgrids,tmp$L_grid)
  AUZls <- c(AUZls,tmp$AUZ_left)

}


if(Z$R_grid != 0){

  	transformed_AUZr <- AUZrs * Rgrids/Z$R_grid
  	AUZr_per_pval <- sum(transformed_AUZr <= Z$AUZ_right, na.rm=TRUE)/length(Rgrids)  

} else {

	AUZr_per_pval <- "NA"

}

if(Z$L_grid != 0){

	transformed_AUZl <- AUZls * Lgrids/Z$L_grid 
	AUZl_per_pval <- sum(transformed_AUZl <= Z$AUZ_left, na.rm=TRUE)/length(Lgrids)

} else {

	AUZl_per_pval <- "NA"

}

if(AUZr_per_pval == 0){AUZr_per_pval <- paste("<",1/Nper)}
if(AUZl_per_pval == 0){AUZl_per_pval <- paste("<",1/Nper)}


Summary <- cbind(c("ZH","AUZ_right_global","AUZ_left_global","AUZ_right_pval","AUZ_left_pval","AUZ_right_window","AUZ_left_window"),
                 c(Z$Rel_zipper_height,round(Z$AUZ_right,4), round(Z$AUZ_left,4), AUZr_per_pval, AUZl_per_pval,round(Z$AUZr_window,4), round(Z$AUZl_window,4)))
colnames(Summary) <- c("Parameter","Value")

Summary <- tableGrob(data.frame(Summary),rows=NULL)

##############################################################################################

# OUTPUT
pdf(paste(output_filename,'.pdf', sep=''),width=10, height=7)
grid.arrange(arrangeGrob(ZP),arrangeGrob(user_params,Summary),nrow=1,ncol=2,widths=2:0.5)
invisible(dev.off())

write.table(RESULTS, col.names=T, row.names=F, quote=F, sep='\t', file=paste(output_filename, '.txt', sep=''))