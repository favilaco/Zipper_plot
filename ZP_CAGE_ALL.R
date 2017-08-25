##############################################################################################
##   Zipper plot: visualizing transcriptional activity of genomic regions                   ##
##                                                                                          ##
##   Author: Francisco Avila Cobos (Francisco.AvilaCobos@UGent.be)                          ##
##                                                                                          ##
##############################################################################################
# Script usage:
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=6){
  print("Please check that all required parameters are indicated or the paths are correct")
  print("Required parameters after Rscript (in order): Script_file, input_file; path_to_bedtools, Zipper_width, Number_permutations_for_AUZpval, tpm_threshold, output_filename")
  print("Example usage: Rscript ZP_CAGE_ALL.R input.txt /usr/bin/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools 5000 100 0 output_1")
  stop()
} 

input_filename = args[1]
bedtools_path = args[2]
Zipper_width = as.numeric(args[3])
Nper = as.numeric(args[4])
tpm_threshold = as.numeric(args[5])
output_filename = args[6]

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
CAGE_input <- ("CAGE_input/")
mc_table <- fread("ManualCuration_2column.txt")
no_gaps = paste(User_dir,"No_gaps/filtered_",sep="/")

user_params =tableGrob(rbind(c('File',output_filename),c('DataType','CAGE-seq'),c('Sample', 'ALL'),c('Zipper width',Zipper_width),c('TPM threshold',tpm_threshold)),theme=ttheme_minimal(),rows=NULL)

##############################################################################################

setwd(User_dir)

system(paste("awk '!array[$1,$2,$3]++' ",input_filename, "> tmp && mv tmp ",input_filename, sep="")) #Avoid problems when users forget to remove duplicate lines

lnc_names <- read.table(input_filename,fill=TRUE,header=FALSE,check.names=FALSE)
Input_size = countLines(input_filename)
COUNTER = 0 #To avoid dynamically growing structures

RESULTS <- data.frame(matrix(nrow=Input_size,ncol=7))
PERMUTATIONS <- data.frame(matrix(nrow=Input_size,ncol=2*Nper))
USER_INPUT <- data.frame(matrix(nrow=Input_size,ncol=3))
colnames(USER_INPUT) <- c("chr","TSS","strand")


system(paste(paste("sort -k1,2 ",input_filename,sep=""), " > Input_sorted.txt",sep=""))
system("awk \'{print >> $1; close($1)}\' Input_sorted.txt")

system("rm *Input*")

##############################################################################################
set.seed(4)
setwd(User_dir)

for (File in dir(User_dir,pattern="chr")){  #chr1
  counter_plus = 0
  counter_minus = 0
  
  CHR <- File
  
  ### CAGE
  CAGE <- readRDS(paste(CAGE_input,paste(paste("BIGCAGE",CHR,sep="_"),"rds",sep="."),sep=""))
  
  
  #NARROWING CAGE TABLE BASED ON USER CHOICE
  CAGE <- CAGE[rowSums(CAGE[,4:ncol(CAGE),with=FALSE] >= tpm_threshold,na.rm=TRUE)>0,]
  
  setkey(CAGE,strand)
  
  ### User
  User <- data.table(read.table(paste(User_dir,File,sep="/"),fill=TRUE,colClasses=c('character','numeric','character'))[,1:3])
  colnames(User) <- c("chr","TSS","strand")
  setkey(User,strand)
  
  if("+" %in% names(table(User[,3,with=FALSE]))){

    User_plus <- User["+"]
    setkey(User_plus,TSS)
    counter_plus = counter_plus + nrow(User_plus) #*
    CAGE_plus <- CAGE["+"]

  }
  
  if("-" %in% names(table(User[,3,with=FALSE]))){

    User_minus <- User["-"]
    setkey(User_minus,TSS)
    counter_minus = counter_minus + nrow(User_minus) #*
    CAGE_minus <- CAGE["-"]

  }
  
  ###### Keeping track of User input 
  if(exists("User_plus")){

    USER_INPUT[(COUNTER + 1):(COUNTER + counter_plus),] <- User_plus 
    RESULTS[(COUNTER + 1):(COUNTER + counter_plus),] <- closest_CAGE_alltissues(User_plus,CAGE_plus,detailed=TRUE)

  } 
  
  if(exists("User_minus")){

    USER_INPUT[(1 + COUNTER + counter_plus):(COUNTER + counter_plus + counter_minus),] <- User_minus 
    RESULTS[(1 + COUNTER + counter_plus):(COUNTER + counter_plus + counter_minus),] <- closest_CAGE_alltissues(User_minus,CAGE_minus,detailed=TRUE)

  } 
  
  
  ##############################################################################################
  ##############################################################################################
  # 2) PERMUTATIONS FOR ZIPPER SHAPE PVAL:
  
  RL = nrow(User)
  
  ##############################################################################################
  
  Random = paste(CHR, "random.txt",sep="")
  to_filter = paste(paste("filtered",CHR,sep="_"),".txt",sep="") #genomic regions WITHOUT gaps, heterochromatin, centromeres, telomeres and repetitive regions
  
  command1 = paste(bedtools_path, " random -g ", User_dir, "/ChromInfo.txt -l 1 -seed 4 -n ", as.integer(RL*Nper), sep="")
  command1 = paste(paste(command1, ">", sep = " "), "temp.txt" , sep=" ")
  system(command1, wait=TRUE)
  
  #No_gaps folder contained the filtered version of the genome (No gaps, centromeres, telomeres, heterochromatin, repetitive elements...)
  command2 = paste(bedtools_path," shuffle -g ", User_dir ,"/ChromInfo.txt -i temp.txt -seed 4 -chrom -incl ", paste(no_gaps,CHR,".txt",sep=""), sep="")
  command2 = paste(paste(command2, " | cut -f2 >", sep = " "), Random , sep=" ") 
  system(command2, wait=TRUE)
  
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

    if(exists("pt_plus")){Pp[,c(2*i-1,2*i)] <- as.matrix(closest_CAGE_alltissues(pt_plus[,c(1,i+2,2),with=FALSE],CAGE_plus,detailed=FALSE))}
    if(exists("pt_minus")){Pm[,c(2*i-1,2*i)] <- as.matrix(closest_CAGE_alltissues(pt_minus[,c(1,i+2,2),with=FALSE],CAGE_minus,detailed=FALSE))}
  
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
  
  COUNTER =  COUNTER + counter_plus + counter_minus
  
  if(exists("User_plus")){rm(list="User_plus")} #Needed for cases when only 1 strand features are present!!
  if(exists("User_minus")){rm(list="User_minus")}


system(paste("rm","temp.txt",Random,sep=" "))
system(paste("rm",File))
rm(list="CAGE")

}

##############################################################################################

#Since some chromosomes have only features on "+" or "-" strand:
PERMUTATIONS <- data.table(na.omit(PERMUTATIONS))

colnames(RESULTS) <- c("chr","TSS","strand","dist_CAGE_start","dist_CAGE_end","tpm","closest_hit")

RESULTS <- data.table(na.omit(RESULTS))
RESULTS$Type <- mc_table[match(RESULTS[,closest_hit,],mc_table[,Name,]),Case,]
setkey(RESULTS,dist_CAGE_start)

##############################################################################################
# ZIPPER PLOT

Z <- ORDERING_CAGE(RESULTS[,3:5,with=FALSE],Plot=TRUE)
ZP <- CAGE_ZP(Z,Zipper_width,COLOUR)

##############################################################################################

# p-values and Summary table
# PERMUTATIONS: calculated while looping through chromosomes!


lnc_key <- do.call(paste,lnc_names[,1:3])
RESULTS_key <- do.call(paste,RESULTS[,1:3,with=FALSE])
RESULTS$name <- lnc_names[as.vector(sapply(RESULTS_key,function(x) grep(x,lnc_key))),4]
setcolorder(RESULTS,c(ncol(RESULTS),1:(ncol(RESULTS)-1)))

# Change working directory to create ouput files in user-specific dir
setwd(User_dir)

Rgrids <- vector()
AUZrs <- vector()
Lgrids <- vector()
AUZls <- vector()

PERMUTATIONS <- data.table(data.frame(strand=USER_INPUT$strand, PERMUTATIONS))

for (j in seq(2,ncol(PERMUTATIONS),2)){

  tmp <- ORDERING_CAGE(PERMUTATIONS[,c(1,j,j+1),with=FALSE])
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
