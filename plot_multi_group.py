#!/usr/bin/env python
#coding:utf-8
__description__ = '''This Script is disigned to compare different otus or taxas from different groups by using
		     group_significance.py in QIIME 1.9.1, with -s Test: nonparametric_t_test, bootstrap_mann_whitney_u, 
		     ANOVA, kruskal_wallis, g_test, parametric_t_test, mann_whitney_u [default: kruskal_wallis]
                     initial otu_taxa_table in the first column and group of this 
                     sample in the second column; a taxonomy level which limit the
                     output to otu_taxa_table_L1-"level";'''

import os,argparse
parser = argparse.ArgumentParser(description='group bar ')
parser.add_argument("-i", "--input", help='a txt file from output of group_significance.py',required=True)
parser.add_argument("-w", "--width", default = "8",help='figure width')
parser.add_argument("-f_h", "--height", default = "8",help='figure height')

args = parser.parse_args( )
infile = args.input
width = args.width
height = args.height

rscript='''
df<-read.table(file="'''+infile+'''",header=T,check.names=FALSE,sep="\\t",comment.char=\"\",quote=\"\\t\")
head.keep<-colnames(df)[1:5]
df<-df[df$P<0.05,]
library(dplyr)
library(reshape2)

df.melt<-melt(df,id.vars=head.keep,
              value.name = "Relative_abundance",variable.name = 'Treatment')
df.melt$mean_std<-ifelse(grepl( "_mean",df.melt$Treatment) ,'mean','std')
#df.melt$Treatment<-gsub("_.*","",df.melt$Treatment)
df.melt$Treatment<-gsub("_mean","",df.melt$Treatment)
df.melt$Treatment<-gsub("_std","",df.melt$Treatment)
df.melt_mean_std<-split( df.melt , f = df.melt$mean_std )
#str(df.melt_mean_std)
df.melt.1<-merge(df.melt_mean_std[[1]][,c(head.keep,'Treatment','Relative_abundance')],
                 df.melt_mean_std[[2]][,c(head.keep,'Treatment','Relative_abundance')],
                 by=c(head.keep,'Treatment'),
                 suffixes = c(".mean",".std"),
                 all = TRUE)
levels <- c(-Inf, 0.001, 0.01, 0.05)
labels <- c("***","**","*")
df.melt.1$P_star <- cut(df.melt.1$P, levels,labels)
df.melt.1$P= round(df.melt.1$P, digits = 4)
df.melt.2 <- as.data.frame(df.melt.1 %>% arrange(OTU))
 

library(ggplot2)
pdf("multi_comparison.pdf",width = '''+width+''', height = '''+height+''')

ggplot(data=df.melt.2, aes(x=OTU, y=Relative_abundance.mean, fill=Treatment)) +
  geom_bar(stat="identity",position = position_dodge())+
  geom_errorbar(aes(ymin=Relative_abundance.mean-Relative_abundance.std, 
                    ymax=Relative_abundance.mean+Relative_abundance.std), width=1,
                position=position_dodge(.9))+
  geom_text(aes(x=OTU,y=max(Relative_abundance.mean)*1.05,label=P_star), vjust=0.7, color="red", 
            size=3.5)+
  geom_text(aes(x=OTU,y=max(Relative_abundance.mean)*1.2,label=P), vjust=0.7, color="red", 
          size=3)+
  scale_fill_brewer(palette="Paired")+coord_flip()+
  theme_minimal()
dev.off()

'''
with open('group_bar.r','w') as wR:
        wR.writelines(rscript)
os.system("Rscript group_bar.r")
