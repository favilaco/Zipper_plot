General Information
===================
This code was created by Francisco Avila Cobos.

There is a webtool implementation available at http://zipperplot.cmgg.be

Feel free to provide feedback/inform about issues by sending an email to: Francisco.AvilaCobos *(at)* UGent *(dot)* be

Citing the Zipper plot
======================
If used, please cite our work:

Avila Cobos, F. et al. Zipper plot: visualizing transcriptional activity of genomic regions. bioRxiv 077073; doi: https://doi.org/10.1101/077073 

Step 1: Installing software requirements
========================================
BEDTools v2.26.0:
----------------- 

Option 1) 

* First: go to your browser and download:

```
https://github.com/arq5x/bedtools2/archive/6f9c61fa34c077082ca9f8785992f3804c210e5d.zip
```

* Second: don't forget to **move the .zip file you have just downloaded to the folder where you want to run the analysis**.

* Third: Go to the folder in which you want to run your analysis and type:
```
$ unzip 6f9c61fa34c077082ca9f8785992f3804c210e5d.zip
$ cd bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d
$ make
```

*NOTE*: Once this is complete, the **path to bedtools should be**:
```
/your_folder/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools
```
Where "/your_folder/" could be, for example, "/Users/favilaco/Desktop/"



Option 2) 

* Go to the folder in which you want to run your analysis and type on the command line:

```
$ wget https://github.com/arq5x/bedtools2/archive/6f9c61fa34c077082ca9f8785992f3804c210e5d.zip
$ unzip https://github.com/arq5x/bedtools2/archive/6f9c61fa34c077082ca9f8785992f3804c210e5d.zip
$ cd bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d
$ make

```
*NOTE*: Once this is complete, the **path to bedtools should be**:
```
/your_folder/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools
```
Where "/your_folder/" could be, for example, "/Users/favilaco/Desktop/"


R packages: 
---------- 
*NOTE*: R version 3.3.1 should be already installed

Required packages: ggplot2_2.2.1, R.utils_2.5.0, data.table_1.10.4, gridExtra_2.2.1   

* Paste and run the following commands **in R**: 

```
packages <- c("ggplot2","R.utils","data.table","gridExtra")
for (i in packages){ install.packages(i, character.only = TRUE)}
```


Step 2: Cloning the Zipper plot Github repository to your computer/server
=========================================================================
Go to the folder in which you want to run your analysis and type **on the command line**:

```
$ git clone https://github.com/favilaco/Zipper_plot
```


Step 3: Download the filtered and pre-processed input files:
============================================================
**Move to the folder** Zipper_plot and, once you are there, download the following file by running the following **on the command line**:

```
$ cd Zipper_plot
$ wget https://www.dropbox.com/s/t6g3fb8atcwjp8c/InputFiles_ZP.zip
```

Unzip "InputFiles_ZP.zip" from the **command line**, by running the commands:

```
$ unzip InputFiles_ZP.zip
```

Step 4: Preparing the input file with the correct format:
=========================================================
Please create a text file containing **4 tab separated columns** as in this box (1 line per genomic feature):

```
chr1	17383802	-	lnc-name1
chr2	1936171	+	lnc-name2
...
```

**NOTE**: Please do not forget to include "chr" as part of your chromosome name.

As an example, you can download an example file (input.txt) by running the following **on the command line**: 
```
wget https://www.dropbox.com/s/cjvh9b414qpoadg/input.txt
```

Example usage:
==============
Once you completed all the above steps and have **in the same folder**: *the unzipped version of "InputFiles_ZP", your input file and the R scripts*, you can run one of the following **four options [see our manuscript for detailed information]**:

FANTOM5 **(CAGE-seq)** with **ALL** sample types:
-------------------------------------------------

```
Rscript ZP_CAGE_ALL.R input.txt /usr/bin/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools 5000 100 0 output_1
```
Where the parameters after "Rscript" (in order): 
**Script_file, input_file; path_to_bedtools, Zipper_width, Number_permutations_for_AUZpval, tpm_threshold** and **output_filename**



FANTOM5 **(CAGE-seq)** with **1 sample types** (e.g.):
------------------------------------------------------
```
Rscript ZP_CAGE_1tissue.R input.txt /usr/bin/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools 5000 100 0 CNhs10723.10400-106A4 OFF output_1
```
Where the parameters after "Rscript" (in order): **Script_file, input_file; path_to_bedtools, Zipper_width, Number_permutations_for_AUZpval, tpm_threshold, sample_type, Tissue_pval** and **output_filename**

*NOTE*: ALL possible sample types can be found in the file **"ManualCuration_2column.txt"**
	
For example:  **"CNhs10625.10019-101D1"** *(lung, adult, pool1)*


Roadmap **(ChIP-seq/DNase-seq)** with **ALL** sample types:
-----------------------------------------------------------
```
Rscript ZP_ROADMAP_ALL.R input.txt /usr/bin/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools 5000 100 H3K4me3 narrowPeak output_1"
```

Where the parameters after "Rscript" (in order): **Script_file, input_file; path_to_bedtools, Zipper_width, Number_permutations_for_AUZpval, mark, peak_type** and **output_filename**

You can select:

1) **narrowPeak, broadPeak** or **gappedPeak**

2) H3K4me1, H3K4me2, H3K4me3, H3K36me3, H3K79me2, H4K20me1, H3K9ac, H3K14ac, H3K27ac or DNaseI


Roadmap **(ChIP-seq/DNase-seq)** with **1 sample type**:
--------------------------------------------------------
```
Rscript ZP_ROADMAP_1tissue.R input.txt /usr/bin/bedtools2-6f9c61fa34c077082ca9f8785992f3804c210e5d/bin/bedtools 5000 100 E096-H3K4me3 narrowPeak OFF output_1
```

Where the parameters after "Rscript" (in order): **Script_file, input_file; path_to_bedtools, Zipper_width, Number_permutations_for_AUZpval, tissue-mark, peak_type, Tissue_pval** and **output_filename**


You can select:

1) **narrowPeak, broadPeak** or **gappedPeak**

2) H3K4me1, H3K4me2, H3K4me3, H3K36me3, H3K79me2, H4K20me1, H3K9ac, H3K14ac, H3K27ac or DNaseI

3) Any tissue that can be found in the file **"ROADMAP_info.txt"**.

* For example: **"E096"** *(lung)*
