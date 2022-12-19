#!/usr/bin/Rscript
# This is a script to intake vcf.gz genotypes and output a png zygosity map

hom_geno <- function (sample,output_dir) { 
	library(quantsmooth)
	library(ggplot2)
	sample_name<-tail(unlist(strsplit(sample,split='/')),1)
	input<-read.table(sample,sep="\t",head=T)
#	input<-(system(intern=T, command=paste("zcat", sample, "|","awk","'$1 ~/(chr?)[1-9][0-9]?/ && $7 ~/PASS|^[.]$/ && $10 ~ /[01]\\/[01]/'", "|", "cut -f 1,2,10")))
	chrs<-input[,1];
	chrs<-gsub("chr","",chrs)
	positions<-input[,2]
	genotype<-grep("/",unlist(strsplit(as.character(input[,10]),":")),value=T)
#	genotype<-grep("/",unlist(strsplit(as.character(input[,3]),":")),value=T)
	genotype<-gsub("1/1",1,genotype); 
	genotype<-gsub("0/1",0.5,genotype); 
	genotype<-gsub("1/0",0.5,genotype); 
	genotype<-gsub("0/0",0,genotype); 

	zygosity<-as.numeric(genotype)
	hom_df<-data.frame(chrs,positions,zygosity)

	CHR<-hom_df[,1]
	MapInfo<-hom_df[,2]
	plot_title<-as.character(sample_name)
	plot_title<-paste(output_dir, plot_title,sep="/");
	file_out_png<-gsub("vcf.gz","zygosity.plot.png",plot_title)
	
	png(file=file_out_png,width=11,height=8,units="in",res=72)
	chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, organism="hsa",main=plot_title,topspace=0)
	chrompos_zyg<-cbind(chrompos,hom_df[,3])
	hom_ref<-subset(chrompos_zyg,chrompos_zyg[,3]==0)
	het<-subset(chrompos_zyg,chrompos_zyg[,3]==0.5)
	hom_alt<-subset(chrompos_zyg,chrompos_zyg[,3]==1)

	#plot hom_alt
		points(x=hom_alt[,2],y=hom_alt[,1]+hom_alt[,3]/4+0.125,pch=".",cex=0.001,col="black")
	#plot het
		points(x=het[,2],y=het[,1]+het[,3]/4+0.125,pch=".",cex=0.001,col="red")
	#plot hom_ref
		points(x=hom_ref[,2],y=hom_ref[,1]+hom_ref[,3]/4+0.125,pch=".",cex=0.001,col="black")

	dev.off()
	
}
args = commandArgs(trailingOnly = TRUE);
my_sample<-args[[1]]
output_dir<-args[[2]]

hom_geno(my_sample,output_dir)


#hom_geno("WTCCC113674.homo.geno.vcf.gz")
