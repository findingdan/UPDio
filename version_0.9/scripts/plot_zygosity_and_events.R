#!/usr/bin/Rscript
# This is a script to intake vcf.gz genotypes and output a png zygosity map


## Subroutines

load_library <- function() { 
	library(quantsmooth)
	library(ggplot2)
}

get_sample_name <- function () {
	after_last_slash<-tail(unlist(strsplit(as.character(events_list_path),split='/')),1)
	name<-gsub(".events_list","",after_last_slash)
	return(name)	
}

use_zcat_or_cat <- function (vcf) {
	split_by_period<-as.vector(cbind(unlist(strsplit(vcf,split="\\."))))
	last_two_fields<-tail(split_by_period,n=2)
	if (last_two_fields[1]=="vcf" & last_two_fields[2] =="gz") {
		return ("zcat") 
	} 
	else if (length(last_two_fields) == 2 & last_two_fields[2] == "vcf") {
		return ("cat") 
	} 
	else { stop() }
}	

get_genotype <- function(genotypes) {
	genotype<-grep("/",unlist(strsplit(genotypes,":")),value=T)
	genotype<-gsub("/","",genotype)
#	genotype<-gsub("1/1",0,genotype); 
#	genotype<-gsub("0/1",1,genotype); 
#	genotype<-gsub("1/0",1,genotype); 
#	genotype<-gsub("0/0",0,genotype); 
	genotype<-gsub("11",0,genotype); 
	genotype<-gsub("01",1,genotype); 
	genotype<-gsub("10",1,genotype); 
	genotype<-gsub("00",0,genotype); 
	return(genotype)
}

convert_event_to_numeric <- function(event) {

	event<-sub("B",0,event); 
	event<-sub("MI_S",0.5,event); 
	event<-sub("MI_D",0.5,event); 
	event<-sub("AI_P",1,event); 
	event<-sub("UI_P",1.5,event); 
	event<-sub("AI_M",2,event); 
	event<-sub("UI_M",2.5,event); 

#	event<-sub("BPI",0,event); 
#	event<-sub("MI_S",0.5,event); 
#	event<-sub("MI_D",0.5,event); 
#	event<-sub("hUPI_P,iUPI_P",1,event); 
#	event<-sub("iUPI_P",1.5,event); 
#	event<-sub("hUPI_M,iUPI_M",2,event); 
#	event<-sub("iUPI_M",2.5,event); 
	return(event)
}

load_vcf <- function(vcf) {
	awk_command<-paste(sep=" ", "awk" ,"'$1 ~/[XY]|(chr)?[1-9][0-9]?/","&&","$4 ~/^[ATCG]$/", "&&", "$5 ~/^[ATCG]$/", "&&","$7 ~/PASS|^[.]$/","&&", "$10 !~ /indel/" , "&&", "$10 ~ /[01][/]?[01]/", "{print $1,$2,$10}'")	
	sed_command<-"sed -e s/[Cc]hr//"	
	chrom_pos_zyg<-data.matrix(do.call(rbind,strsplit((system(intern=T, command=paste(use_zcat_or_cat(vcf), vcf, "|", awk_command, "|", sed_command)))," ")))
	chrom_pos_zyg[,1]<-numericCHR(chrom_pos_zyg[,1])
	chrom_pos_zyg[,3]<-get_genotype(chrom_pos_zyg[,3])
	chrom_pos_zyg<-data.matrix(data.frame(chrom_pos_zyg,stringsAsFactors=F))
	colnames(chrom_pos_zyg)<-c("CHR","MapInfo","zygosity")
	return(chrom_pos_zyg)
}

load_events <- function(in_file) {
	# awk_command<-paste(sep=" ", "awk" ,"'{print $1,$2,$8}'")	
	awk_command<-paste(sep=" ", "awk" ,"'$8 !~ /uninformative/", "{print $1,$2,$8}'")	
	sed_command<-"sed -e s/[Cc]hr//"	
	chrom_pos_event<-data.matrix(do.call(rbind,strsplit((system(intern=T, command=paste("cat", in_file, "|", awk_command, "|", sed_command)))," ")))
	chrom_pos_event[,1]<-numericCHR(chrom_pos_event[,1])
	chrom_pos_event[,3]<-convert_event_to_numeric(chrom_pos_event[,3])
	chrom_pos_event<-data.matrix(data.frame(chrom_pos_event,stringsAsFactors=F))
	colnames(chrom_pos_event)<-c("CHR","MapInfo","event")
	return(chrom_pos_event)
}

get_matching_chr_pos <- function (items_to_match, bigger_set) {
	index1<-paste(items_to_match[,1],items_to_match[,2])
	index2<-paste(bigger_set[,1], bigger_set[,2])
	return(bigger_set[which(index2%in%index1), ])
 }

numericCHR <- function(CHR) {
	CHR <- as.character(CHR)
	CHR[CHR== "X"] <- "98"
	CHR[CHR==23] <- "98"
	CHR[CHR== "Y"] <- "99"
	as.numeric(CHR)
}

characterCHR <- function(CHR) {
    CHR <- as.character(CHR)
	CHR[CHR==23] <- "98"
    CHR[CHR == "98"] <- "X"
    CHR[CHR == "99"] <- "Y"
    CHR[CHR == "100"] <- "XY"
    CHR[CHR == "101"] <- "MT"
    CHR
}

convert_chr_to_match_GenomePlot <-function (df_to_convert,sexChromosomes=plot_sex) {
	chrom.n <-22
	chrs2 <- factor(numericCHR(df_to_convert[, "CHR"]), levels = c(1:chrom.n, if (sexChromosomes) c("X", "Y") else NULL))
	#chrs2 <- factor(numericCHR(df_to_convert[, "CHR"]), levels = c(1:chrom.n, if (sexChromosomes) c("X") else NULL))
	lens <- lengthChromosome(levels(chrs2),units="bases")
	names(lens) <- characterCHR(names(lens))

	dwidth <- NULL
	for (i in 1:(chrom.n%/%2)) {
		dwidth[i] <- lens[i] + lens[chrom.n + 1 - i]
	}
	if (chrom.n%%2 == 1) {
		dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
	}
	if (sexChromosomes) {
		dwidth <- c(dwidth, lens["X"] + lens["Y"])
	}

	maxdwidth <- max(dwidth) * 1.05
	leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 1)%/%2):1)
	rightrow <- c(if (sexChromosomes) "Y" else NULL, if (chrom.n%%2 == 1) "" else NULL, ((chrom.n + 1)%/%2 + 1):chrom.n)

	dchrompos <- matrix(0, nrow = nrow(df_to_convert), ncol = 3, dimnames = list(rownames(df_to_convert), c("CHR", "MapInfo","zygosity")))
	for (i in 1:length(rightrow)) if (rightrow[i] != "") {
		probes <- characterCHR(df_to_convert[, "CHR"]) == rightrow[i]
		dchrompos[probes, 2] <- df_to_convert[probes, "MapInfo"] + maxdwidth - lens[rightrow[i]]
		dchrompos[probes, 1] <- i
		dchrompos[probes, 3] <- df_to_convert[probes, 3]
	}
	for (i in 1:length(leftrow)) {
		probes <- characterCHR(df_to_convert[, "CHR"]) == leftrow[i]
		dchrompos[probes, 2] <- df_to_convert[probes, "MapInfo"]
		dchrompos[probes, 1] <- i
		dchrompos[probes, 3] <- df_to_convert[probes, 3]
	}
	return(dchrompos)
}

# why can't i put the subroutines at the end of the script, like in perl?

#MAIN

load_library() 

events_list_path<-commandArgs(trailingOnly = TRUE)[[1]]
vcf<-commandArgs(trailingOnly = TRUE)[[2]]
output_dir<-commandArgs(trailingOnly = TRUE)[[3]]
plot_sex<-commandArgs(trailingOnly=TRUE)[[4]]

sample_name <- get_sample_name()

# Load events
chrom_pos_events<-load_events(events_list_path)
# Only take positions from vcf if they are present in the events file; disabling for now
chrom_pos_zyg<-get_matching_chr_pos(chrom_pos_events,load_vcf(vcf))
# chrom_pos_zyg<-chrom_pos_events
# Convert plotting coordinates to match the quantsmooth paradigm
chrom_pos_zyg_converted<-convert_chr_to_match_GenomePlot(chrom_pos_zyg)
chrom_pos_events_converted<-convert_chr_to_match_GenomePlot(chrom_pos_events)

#prepare data for graphing
# prepare zygosity
hom<-subset(chrom_pos_zyg_converted,chrom_pos_zyg_converted[,3]==0)
het<-subset(chrom_pos_zyg_converted,chrom_pos_zyg_converted[,3]==1)

# prepare Genotype events
BPI<-subset(chrom_pos_events_converted,chrom_pos_events_converted[,3]==0)
MI<-subset(chrom_pos_events_converted,chrom_pos_events_converted[,3]==0.5)
hUPI_P<-subset(chrom_pos_events_converted,chrom_pos_events_converted[,3]==1)
iUPI_P<-subset(chrom_pos_events_converted,chrom_pos_events_converted[,3]==1.5)
hUPI_M<-subset(chrom_pos_events_converted,chrom_pos_events_converted[,3]==2)
iUPI_M<-subset(chrom_pos_events_converted,chrom_pos_events_converted[,3]==2.5)



# Plotting

sample_name_events_plot_png<-paste(sep="", as.character(sample_name), ".events_plot.png")
file_out_png<-paste(output_dir,sample_name_events_plot_png,sep="/")
png(file=file_out_png,width=1200,height=700,res=75)

# Plot using the vcf coordinates, since the event coordinates are a subset of the vcf coordiantes
x<-prepareGenomePlot(chrom_pos_zyg,paintCytobands = TRUE, organism="hsa",main=sample_name,topspace=0,sexChromosomes=plot_sex)
#legend(x="bottom",legend=rev(c("Hom","Het","BPI","MI","iUPI_M","hUPI_M","iUPI_P","hUPI_P")),fill=rev(c("black","red", "black","grey","darkgreen","lightgreen","darkblue","lightblue")))
legend(x="bottom",legend=rev(c("Hom","Het","B","MI","UI_M","AI_M","UI_P","AI_P")),fill=rev(c("black","red", "black","grey","darkgreen","lightgreen","darkblue","lightblue")))

#Plot Zygosity of Informative Events
points(x=hom[,2],y=hom[,1]+0.1,pch=".",cex=1,col="black")
points(x=het[,2],y=het[,1]+0.18,pch=".",cex=1,col="red")

#Plot Genotype Events
points(x=BPI[,2],y=BPI[,1]			+ 0.26,pch=".",cex=2,col="black")
points(x=MI[,2],y=MI[,1]			+ 0.32,pch=".",cex=2,col="grey")
points(x=hUPI_M[,2],y=hUPI_M[,1] 	+ 0.4,pch=".",cex=2,col="darkgreen")
points(x=iUPI_M[,2],y=iUPI_M[,1] 	+ 0.48,pch=".",cex=2,col="lightgreen")
points(x=hUPI_P[,2],y=hUPI_P[,1] 	+ 0.56,pch=".",cex=2,col="darkblue")
points(x=iUPI_P[,2],y=iUPI_P[,1] 	+ 0.64,pch=".",cex=2,col="lightblue")

dev.off()

## rm(chrom_pos_zyg) <-this does not seem to release the memory :-( , 
## and gc() didn't seem to work. why R, rhy.


#for testing
#events_list_path<-"/lustre/scratch107/projects/ddd/users/dk6/projects/upd/output/ddd/exome/DDD_MAIN5250962.events_list"
#vcf<-"/lustre/scratch107/projects/ddd/users/dk6/data/exomes/indiv_samples/DDD_MAIN5250962.vcf.gz"
#output_dir<-"~/programs/scripts/upd"
