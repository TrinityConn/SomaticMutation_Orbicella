###Script for Analyzing Orbicella Mutation Data ##



##load in VAF Files ##


#load in required packages 
library(ggplot2)
library(tidyverse)
library(data.table)
library(cowplot)
library(vegan)
library(TSA)
library(wesanderson)

#loading in finalized VAF file. These files were exported with the example column names "# [1]CHROM", "[2]POS" "[3]145_1:AF",
# "[4]145_1:DP", "[5]145_3:AF", "[6]145_3:DP"

#set working directory 

setwd("/Users/trini/Documents/PSU/Orbicella Project/Text_Files")


#read in files in directory 
filelist=list.files(pattern="*_refaltsamp*")

datalist=lapply(filelist, function(x)read.table(x, header=FALSE))

#add in column for filename 
for (i in 1:length(datalist)){datalist[[i]]<-cbind(datalist[[i]],filelist[i])}

#create function to label columns based on if column 5 sum is greater than column 3 
tumor2<-function(x) {if (sum(x$V7) > sum(x$V5)) {
  colnames(x)[5]="Normal_AF"
  colnames(x)[6]="Normal_DP"
  colnames(x)[7]="Tumor_AF"
  colnames(x)[8]="Tumor_DP"
  x
}
}

#create function to label columns based on if column 3 is greater than column 5 
tumor3<-function(x) {if (sum(x$V5) > sum(x$V7)) {
  colnames(x)[5]="Tumor_AF"
  colnames(x)[6]="Tumor_DP"
  colnames(x)[7]="Normal_AF"
  colnames(x)[8]="Normal_DP"
  x
}
}


#apply function tumor3 to data list
data2=lapply(datalist, tumor3) 


#apply function tumor2 to data list 
data3=lapply(datalist, tumor2)

#combine the two lists 
data4<-c(data2,data3)

#remove duplicate 0 length files 
data5<-data4[lapply(data4, length)>0]


#combine rows of all imported files
data_old=do.call("rbind", data5)

#add in column names 

colnames(data_old)[1]="CHRO" 
colnames(data_old)[2]="POS"
colnames(data_old)[3]="REF"
colnames(data_old)[4]="ALT"


#add in normal v tumor classifications 
#load in metadata

#load out names to add to metadata sheet 
data_names<-as.data.frame(unique(data_old$`filelist[i]`))
colnames(data_names)[1]="FileNames"
write.csv(data_names, file="Data_Names.csv")

#load in metadata 
metadata<-read.csv("./Mutect_Metadata .csv")


colnames(metadata)[1]="filelist[i]"

#combine metadata and mutect results 

data<-merge(metadata, data_old, by=c("filelist[i]"))

data_2<-data%>%
  filter(Normal!="20007")
data_3<-data_2%>%
  filter(Tumor!="20007")
data_4<-data_3%>%
  filter(Normal!="20021")
data_5<-data_4%>%
  filter(Tumor!="20021")
data_6<-data_5%>%
  filter(Normal!="20026")
data_7<-data_6%>%
  filter(Tumor!="20026")
#get total mutation load 

#filter by depth 
whole_filter<-data_7%>%
  filter(Tumor_AF<=0.7)

#substitution burden per tumor sample 

burdenwholeu<-unique(whole_filter)

#filter to fixed mutations 

burdenfixed<-burdenwholeu%>%
  filter(Tumor_AF>=0.5)

burdenfixed2<-unique()

burdenfixed_count<-burdenfixed%>%
  group_by(Tumor)%>%
  summarise(countfixed=n_distinct(CHRO,POS))

burdenmosaic<-burdenwholeu%>%
  filter(Tumor_AF<=0.49)%>%
  filter(Tumor_AF>0.1)

burdenmosaic_count<-burdenmosaic%>%
  group_by(Tumor)%>%
  summarise(countmosaic=n_distinct(CHRO,POS))

burdenlowfreq<-burdenwholeu%>%
  filter(Tumor_AF<=0.1)

burdenlowfreq_count<-burdenlowfreq%>%
  group_by(Tumor)%>%
  summarise(count_lowfreq=n_distinct(CHRO,POS))

burden_all<-merge(burdenwhole1, burdenlowfreq_count, by=c("Tumor"))
burden_all<-merge(burden_all, burdenfixed_count, by=c("Tumor"))
burden_all<-merge(burden_all, burdenmosaic_count, by=c("Tumor"))

burden_all$propfixed<-burden_all$countfixed/burden_all$count_whole
burden_all$propmosaic<-burden_all$countmosaic/burden_all$count_whole
burden_all$proplowfreq<-burden_all$count_lowfreq/burden_all$count_whole

write.csv(burden_all2, file="allelefreq_stats.csv")

burden_all2<-merge(metadata, burden_all, by=c("Tumor"))

propfixed<-burden_all2%>%
  group_by(Colony)%>%
  summarise(mean=mean(propfixed))

proplowfreq<-burden_all2%>%
  group_by(Transect)%>%
  summarise(mean=mean(proplowfreq))

anova1<-aov(propmosaic~Colony*Transect, data=burden_all2)

burdenall2<-unique(burden_all2$Tumor)

burdenwhole1<-burdenwholeu%>%
  group_by(Tumor)%>%
  summarise(count_whole=n_distinct(CHRO,POS))
burdenfixed_2$Colony<-as.character(burdenfixed_2$Colony)

burdenfixed_2<-merge(burdenfixed_count, metadata, by=c("Tumor"))

ggplot(burdenfixed_2, aes(x=Colony, y=count, fill=Transect))+geom_boxplot()+theme_bw()

#vaf by tumor site and colony 

#plot 
ggplot(MutationClass_Proportions, aes(x=Mutation_Class, y=100*Proportion, fill=Transect))+
  geom_boxplot()+stat_boxplot(geom = 'errorbar')+theme_bw()+ylab("Proportion of Mutations")+xlab("Mutation Class")

MutationClass_Proportions$Mutation_Class<-factor(MutationClass_Proportions$Mutation_Class, levels=c("Low Frequency", "Mosaic", "Fixed"))


#k-s tests 

#start with within colony 3 
upper<-burdenwholeu%>%
  filter(Transect=="Upper")
lower<-burdenwholeu%>%filter(Transect=="Lower")
colony3<-upper%>%
  filter(Colony=="3")
colony4<-upper%>%
  filter(Colony=="4")
colony5<-upper%>%
  filter(Colony=="5")

ggplot(D_values_colony, aes(x=Comparison, y=D, color=Comparison))+geom_point(size=6)+theme_bw()
D_values_colony$Comparison<-factor(D_values_colony$Comparison, levels=c("Three-Four", "Three-Five", "Four-Five"))


ks.boot(colony3$Tumor_AF, colony4$Tumor_AF, nboots=1000)
burdenwholeu$mutid<-paste(burdenwholeu$CHRO, burdenwholeu$POS)
burdenfixed$mutid<-paste(burdenfixed$CHRO, burdenfixed$POS)
#upset plot of fixed mutations 

colony5_fixed<-burdenwholeu

#starting with M145
C5_1<-colony5_fixed%>%
  filter(Tumor=="20023")

C5_1_n<-C5_1%>%
  count(mutid, Tumor)

C5_1_counts<-pivot_wider(C5_1_n, names_from = Tumor, values_from = n)


C5_2<-colony5_fixed%>%
  filter(Tumor=="20034")

C5_2_n<-C5_2%>%
  count(mutid, Tumor)
C5_2_counts<-pivot_wider(C5_2_n, names_from = Tumor, values_from = n)

C5_3<-colony5_fixed%>%
  filter(Tumor=="10019")

C5_3_n<-C5_3%>%
  count(mutid, Tumor)
C5_3_counts<-pivot_wider(C5_3_n, names_from = Tumor, values_from = n)

C5_4<-colony5_fixed%>%
  filter(Tumor=="20024")

C5_4_n<-C5_4%>%
  count(mutid, Tumor)
C5_4_counts<-pivot_wider(C5_4_n, names_from = Tumor, values_from = n)

C5_5<-colony5_fixed%>%
  filter(Tumor=="20022")
C5_5_n<-C5_5%>%
  count(mutid, Tumor)
C5_5_counts<-pivot_wider(C5_5_n, names_from = Tumor, values_from = n)

C5_6<-colony5_fixed%>%
  filter(Tumor=="20018")

C5_6_n<-C5_6%>%
  count(mutid, Tumor)
C5_6_counts<-pivot_wider(C5_6_n, names_from = Tumor, values_from = n)

C5_7<-colony5_fixed%>%
  filter(Tumor=="20029")
C5_7_n<-C5_7%>%
  count(mutid, Tumor)
C5_7_counts<-pivot_wider(C5_7_n, names_from = Tumor, values_from = n)


C5_upset<-full_join(C5_1_counts,C5_2_counts)
C5_upset2<-full_join(C5_upset, C5_3_counts)
C5_upset3<-full_join(C5_upset2, C5_4_counts)
C5_upset4<-full_join(C5_upset3, C5_5_counts)
C5_upset5<-full_join(C5_upset4, C5_6_counts)
C5_upset6<-full_join(C5_upset5, C5_7_counts)

#replace NAs with 0 
C5_upset6[is.na(C5_upset6)] <- 0
C5_upset6<-as.data.frame(C5_upset6)
#rename columsn to remove numeric constant 

colnames(C5_upset6)[2]="S20023"
colnames(C5_upset6)[3]="S20034"
colnames(C5_upset6)[4]="S10019"
colnames(C5_upset6)[5]="S20024"
colnames(C5_upset6)[6]="S20022"
colnames(C5_upset6)[7]="S20018"
colnames(C5_upset6)[8]="S20029"

#convert all numbers greater than 1 to 1 
C5_upset6$S20023<-ifelse(C5_upset6$S20023,1,C5_upset6$S20023)
C5_upset6$S20034<-ifelse(C5_upset6$S20034,1,C5_upset6$S20034)
C5_upset6$S10019<-ifelse(C5_upset6$S10019,1,C5_upset6$S10019)
C5_upset6$S20024<-ifelse(C5_upset6$S20024,1,C5_upset6$S20024)
C5_upset6$S20022<-ifelse(C5_upset6$S20022,1,C5_upset6$S20022)
C5_upset6$S20018<-ifelse(C5_upset6$S20018,1,C5_upset6$S20018)
C5_upset6$S20029<-ifelse(C5_upset6$S20029,1,C5_upset6$S20029)

#load upset package
set_vars<-c("S_1", "M145_2", "M145_3M", "M145_4", "M145_5")
C5_6<-upset(C5_upset6, keep.order=T, sets=c("S20023", "S20034", "S10019", "S20024", "S20022", "S20018", "S20029", "S20027"))


##create co-occurence matrix to build tree 


#extract sets 
Data145<-as.data.frame(M145_list$New_data)
Data145<-Data145%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_M145<-rowSums(Data145[,2:6])==1
private_M145<-as.data.frame(private_M145)
private_M145<-private_M145%>%
  mutate(ID=row_number())

M145_Identi<-merge(Data145, private_M145, by=c("ID"))
M145_private<-M145_Identi%>%
  filter(private_M145=="TRUE")
View(M145_private)

#add in colony and mutation type 

M145_private<-M145_private%>%
  add_column(Intersection="Private")

#merge with original dataframe
private_mutations<-merge(whole_filter, M145_private, by=c("mutid"))
private_mutations<-private_mutations%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)

ggplot(private_mutations, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)

#extract shared with two 
#mutations where the row sum is 2
private_M1452<-rowSums(Data145[,2:6])==2
private_M1452<-as.data.frame(private_M1452)
private_M1452<-private_M1452%>%
  mutate(ID=row_number())
M145_Identi2<-merge(Data145, private_M1452, by=c("ID"))
M145_private2<-M145_Identi2%>%
  filter(private_M1452=="TRUE")

#add in colony and mutation type 


M145_private2<-M145_private2%>%
  add_column(Intersection="Two")

twoshare<-merge(whole_filter, M145_private2, by=c("mutid"))
twoshare<-twoshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)+theme_bw()
ggplot(twoshare, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)
M145_Upset_AF<-rbind(private_mutations, twoshare)

#extract shared with three
#mutations where the row sum is 3
private_M1453<-rowSums(Data145[,2:6])==3
private_M1453<-as.data.frame(private_M1453)
private_M1453<-private_M1453%>%
  mutate(ID=row_number())
M145_Identi3<-merge(Data145, private_M1453, by=c("ID"))
M145_private3<-M145_Identi3%>%
  filter(private_M1453=="TRUE")
M145_private3<-M145_private3%>%
  add_column(Intersection="Three")
threeshare<-threeshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)
threeshare<-merge(whole_filter, M145_private3, by=c("mutid"))
ggplot(threeshare, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)
M145_Upset_AF<-rbind(M145_Upset_AF, fourshare)



##starting with colony 3 -- colony 3 smaller, more scalloped in the layering ###
#colony 3 samples 20009, 20014, 20028, , 20033
#colony 4
#colony 4 samples 20032, 20010, 20027, 20015, 20013, 20016, 20011

#colony 5 samples 20023, 20034, 10019, 20024, 20022, 20018, 20029
S20033<-burdenwholeu%>%
  filter(Tumor=="20033")

ggplot(S20033, aes(x=Tumor_AF))+geom_density()+geom_histogram()+theme_bw()

S20028<-burdenwholeu%>%
  filter(Tumor=="20028")

#vaf by tumor site 
S10019<-burdenwholeu%>%
  filter(Tumor=="10019")

ggplot(S10019, aes(x=Tumor_AF))+geom_density()+geom_histogram(binwidth =0.03)+theme_bw()+xlab("Mutation Allele Frequency, Sample 10019 (Top)")

S20018<-burdenwholeu%>%
  filter(Tumor=="20018")

ggplot(S20018, aes(x=Tumor_AF))+geom_density()+geom_histogram(binwidth =0.03)+theme_bw()+xlab("Mutation Allele Frequency, Sample 20018 (Bottom)")

S20022<-burdenwholeu%>%
  filter(Tumor=="20022")

ggplot(S20022, aes(x=Tumor_AF))+geom_density()+geom_histogram(binwidth =0.03)+theme_bw()+xlab("Mutation Allele Frequency, Sample 20022 (Bottom)")

S20024<-burdenwholeu%>%
  filter(Tumor=="20024")

ggplot(S20024, aes(x=Tumor_AF))+geom_density()+geom_histogram(binwidth =0.03)+theme_bw()+xlab("Mutation Allele Frequency, Sample 20024 (Top))")

S20029<-burdenwholeu%>%
  filter(Tumor=="20029")

ggplot(S20029, aes(x=Tumor_AF))+geom_density()+geom_histogram(binwidth =0.03)+theme_bw()+xlab("Mutation Allele Frequency, Sample 20029 (Bottom)")

S20034<-burdenwholeu%>%
  filter(Tumor=="20034")

ggplot(S20034, aes(x=Tumor_AF))+geom_density()+geom_histogram(binwidth =0.03)+theme_bw()+xlab("Mutation Allele Frequency, Sample 20034 (Top)")


##plotting mutation load by location 

ggplot(Colony5_MutationLoad, aes(x=Location, y=Mutation_Burden, fill=Location))+theme_bw()+geom_violin()+geom_point()+scale_fill_manual(values=wes_palette("GrandBudapest2", n=2
                              
                                                                                                                                                           
#upset plots 

whole_filter$mutid<-paste(whole_filter$CHRO, whole_filter$POS)

C5<-whole_filter%>%filter(Colony=="Five")

#starting with M145
C5_1<-C5%>%
  filter(Tumor=="20034")

C5_1_n<-C5_1%>%
  count(mutid, Tumor)
C5_1_counts<-pivot_wider(C5_1_n, names_from = Tumor, values_from = n)


C5_2<-C5%>%
  filter(Tumor=="20024")
C5_2_n<-C5_2%>%
  count(mutid, Tumor)
C5_2_counts<-pivot_wider(C5_2_n, names_from = Tumor, values_from = n)

C5_3<-C5%>%
  filter(Tumor=="20022")
C5_3_n<-C5_3%>%
  count(mutid, Tumor)
C5_3_counts<-pivot_wider(C5_3_n, names_from = Tumor, values_from = n)

C5_4<-C5%>%
  filter(Tumor=="20018")

C5_4_n<-C5_4%>%
  count(mutid, Tumor)
C5_4_counts<-pivot_wider(C5_4_n, names_from = Tumor, values_from = n)

C5_5<-C5%>%
  filter(Tumor=="20029")
C5_5_n<-C5_5%>%
  count(mutid, Tumor)
C5_5_counts<-pivot_wider(C5_5_n, names_from = Tumor, values_from = n)

C5_6<-C5%>%
  filter(Tumor=="10019")
C5_6_n<-C5_6%>%
  count(mutid, Tumor)
C5_6_counts<-pivot_wider(C5_6_n, names_from = Tumor, values_from = n)

C5_upset<-full_join(C5_1_counts, C5_2_counts)
C5_upset2<-full_join(C5_upset, C5_3_counts)
C5_upset3<-full_join(C5_upset2, C5_4_counts)
C5_upset4<-full_join(C5_upset3, C5_5_counts)
C5_upset5<-full_join(C5_upset4, C5_6_counts)

#replace NAs with 0 
C5_upset5[is.na(C5_upset5)] <- 0
C5_upset5<-as.data.frame(C5_upset5)
#rename columsn to remove numeric constant 

colnames(C5_upset5)[2]="C5_1"
colnames(C5_upset5)[3]="C5_2"
colnames(C5_upset5)[4]="C5_3"
colnames(C5_upset5)[5]="C5_4"
colnames(C5_upset5)[6]="C5_5"
colnames(C5_upset5)[7]="C5_6"
#convert all numbers greater than 1 to 1 
C5_upset5$C5_1<-ifelse(C5_upset5$C5_1,1,C5_upset5$C5_1)
C5_upset5$C5_2<-ifelse(C5_upset5$C5_2,1,C5_upset5$C5_2)
C5_upset5$C5_3<-ifelse(C5_upset5$C5_3,1,C5_upset5$C5_3)
C5_upset5$C5_4<-ifelse(C5_upset5$C5_4,1,C5_upset5$C5_4)
C5_upset5$C5_5<-ifelse(C5_upset5$C5_5,1,C5_upset5$C5_5)
C5_upset5$C5_6<-ifelse(C5_upset5$C5_6,1,C5_upset5$C5_6)
#load upset package
set_vars<-c("C5_1", "C5_2", "C5_3", "C5_4", "C5_5", "C5_6")
C5_list<-upset(C5_upset5, keep.order=T, sets=c("C5_1", "C5_2", "C5_3", "C5_4", "C5_5", "C5_6"))

#extract sets 
DataC5<-as.data.frame(C5_list$New_data)
DataC5<-DataC5%>%
  mutate(ID=row_number())

#mutations where the row sum is 1 
private_C5<-rowSums(DataC5[,2:7])==1
private_C5<-as.data.frame(private_C5)
private_C5<-private_C5%>%
  mutate(ID=row_number())

C5_Identi<-merge(DataC5, private_C5, by=c("ID"))
C5_private<-C5_Identi%>%
  filter(private_C5=="TRUE")
View(M145_private)

#add in colony and mutation type 

C5_private<-C5_private%>%
  add_column(Intersection="Private")

#merge with original dataframe
private_mutations<-merge(whole_filter, C5_private, by=c("mutid"))
private_mutations<-private_mutations%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)

ggplot(private_mutations, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)

#extract shared with two 
#mutations where the row sum is 2
private_C52<-rowSums(DataC5[,2:7])==2
private_C52<-as.data.frame(private_C52)
private_C52<-private_C52%>%
  mutate(ID=row_number())
C5_Identi2<-merge(DataC5, private_C52, by=c("ID"))
C5_private2<-C5_Identi2%>%
  filter(private_C52=="TRUE")

#add in colony and mutation type 


C5_private2<-C5_private2%>%
  add_column(Intersection="Two")

twoshare<-merge(whole_filter, C5_private2, by=c("mutid"))
twoshare<-twoshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)
ggplot(twoshare, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)
M145_Upset_AF<-rbind(private_mutations, twoshare)

#extract shared with three
#mutations where the row sum is 3
private_C53<-rowSums(DataC5[,2:7])==3
private_C53<-as.data.frame(private_C53)
private_C53<-private_C53%>%
  mutate(ID=row_number())
C5_Identi3<-merge(DataC5, private_C53, by=c("ID"))
C5_private3<-C5_Identi3%>%
  filter(private_C53=="TRUE")
C5_private3<-C5_private3%>%
  add_column(Intersection="Three")
threeshare<-threeshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)
threeshare<-merge(whole_filter, C5_private3, by=c("mutid"))
ggplot(threeshare, aes(x=Tumor_AF))+geom_histogram()+theme_bw()+xlim(0,1.0)
M145_Upset_AF<-rbind(M145_Upset_AF, fourshare)

#extract shared with four 
#mutations where the row sum is 4
private_M1454<-rowSums(Data145[,2:6])==4
private_M1454<-as.data.frame(private_M1454)
private_M1454<-private_M1454%>%
  mutate(ID=row_number())
M145_Identi4<-merge(Data145, private_M1454, by=c("ID"))
M145_private4<-M145_Identi4%>%
  filter(private_M1454=="TRUE")
M145_private4<-M145_private4%>%
  add_column(Intersection="Four")


fourshare<-merge(whole_filter, M145_private4, by=c("mutid"))
fourshare<-fourshare%>%
  select(mutid, Tumor_AF, Tumor_DP, Colony, Tumor, Intersection)

M145_Upset_AF1<-M145_Upset_AF%>%
  filter(Intersection=="Private")
M145_Upset_AF4<-M145_Upset_AF_filter%>%
  filter(Intersection=="Four")

ggplot(M145_Upset_AF1, aes(x=Tumor_AF))+geom_histogram(bins=100
)+xlim(0,1.0)+theme_bw()
M145_Upset_AF_filter<-M145_Upset_AF%>%
  filter(Tumor_DP>=60)


