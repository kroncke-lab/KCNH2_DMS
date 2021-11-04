library(ggplot2)
library(data.table)
library(gtools)
library(dplyr)

setwd("C:/Users/KRONCKE/Box Sync/Kroncke_Lab/kcnh2/sequencing/VANTAGE/Tile1/Library1-DM278/")

# all barcode-related r script removed. I am assuming the barcode was correctly 
# processed.

a=read.csv('barcodeVariant-combined-wCodon.csv',header=TRUE,stringsAsFactors=FALSE)
a["resnum"]<-sub(".*_[A-Z]([0-9]+)[A-Z]_.*","\\1",a$mutation)
a["mutAA"]<-sub(".*_[A-Z][0-9]+([A-Z]).*_.*","\\1",a$mutation)
a["wtAA"]<-sub(".*_([A-Z])[0-9]+[A-Z]_.*","\\1",a$mutation)
a$mutAA[a$wtAA==a$resnum]<-NA
a$resnum[a$wtAA==a$resnum]<-NA
a$wtAA[is.na(a$resnum)]<-NA
a$resnum<-as.numeric(a$resnum)
a["cwt"]<-sub(".*_([A-Z]{3})[0-9]+[A-Z]{3}.*","\\1",a$mutation)
a["cmut"]<-sub(".*_[A-Z]{3}[0-9]+([A-Z]{3}).*","\\1",a$mutation)

# read barcode counts from experiment 
# negative AF647: NO KCNH2 surface expression
neg<-read.table("AF-sort-exp/5239-2/AF-neg",header=FALSE,stringsAsFactors=FALSE)
colnames(neg)<-c("CountNeg","barcode") # added another column name because of silly proc.sh script (hopefully won't need to do in future)

# Low AF647: Low KCNH2 surface expression
low<-read.table("AF-sort-exp/5239-2/AF-low",header=FALSE,stringsAsFactors=FALSE)
colnames(low)<-c("CountLow","barcode") # added another column name because of silly proc.sh script (hopefully won't need to do in future)

# medium AF647: medium KCNH2 surface expression
med<-read.table("AF-sort-exp/5239-2/AF-med",header=FALSE,stringsAsFactors=FALSE)
colnames(med)<-c("CountMed","barcode")

# High AF647: high KCNH2 surface expression
high<-read.table("AF-sort-exp/5239-2/AF-high",header=FALSE,stringsAsFactors=FALSE)
colnames(high)<-c("CountHigh","barcode")

# merge all barcode counts from subassembly and from 
b<-merge(a,neg,all.x = T,all.y = T)
b<-b[!is.na(b$mutation),]
c<-merge(b,low,all.x = T,all.y = T)
c<-c[!is.na(c$mutation),]
d<-merge(c,med,all.x = T,all.y = T)
d<-d[!is.na(d$mutation),]
e<-merge(d,high,all.x = T,all.y = T)
e=e[!is.na(e$mutation),]
e=e[!(is.na(e$CountNeg) & is.na(e$CountLow) & is.na(e$CountMed) & is.na(e$CountHigh)),]
e<-unique(e)

# trying analysis strategies:
e[is.na(e$CountNeg),"CountNeg"]<-0
e[is.na(e$CountLow),"CountLow"]<-0
e[is.na(e$CountMed),"CountMed"]<-0
e[is.na(e$CountHigh),"CountHigh"]<-0
e$CountAll<-e$CountHigh+e$CountMed+e$CountLow+e$CountNeg
e<-e[e$CountAll>4,]

tot_countNeg<-sum(e$CountNeg, na.rm = TRUE)
tot_countLow<-sum(e$CountLow, na.rm = TRUE)
tot_countMed<-sum(e$CountMed, na.rm = TRUE)
tot_countHigh<-sum(e$CountHigh, na.rm = TRUE)
e$CountTotal<-e$CountNeg/tot_countNeg+e$CountLow/tot_countLow+e$CountMed/tot_countMed+e$CountHigh/tot_countHigh
#e$r.max<-pmax(e$CountHigh/e$CountTotal, e$CountLow/e$CountTotal, e$CountMed/e$CountTotal)#, e$CountNeg/e$CountTotal)
#e<-e[e$r.max<1,]
#e<-e[log10(e$CountTotal)<4,]
e$score<-((1*e$CountNeg/tot_countNeg+2*e$CountLow/tot_countLow+3*e$CountMed/tot_countMed+4*e$CountHigh/tot_countHigh)/e$CountTotal)
e$res<-paste(e$resnum+1,e$wtAA)

# somewhat logical ordering of amino acids
aas<-c("A","V","I","L","M","F","Y","W","P","G","C","S","T","N","Q","H","R","K","D","E","X")

x<-data.frame(0,0)
i<-1
colnames(x)<-c("res","mutAA")
for (r in unique(e$res)){
  for (aa in aas){
    x[i,c("mutAA","res")]<-c(aa,r)
    i=i+1
  }
}

e<-bind_rows(e,x)
e<-e[!is.na(e$mutation),]
e$type<-"missense"
e[e$mutAA=="X" & !is.na(e$wtAA),"type"]<-"nonsense"
e[e$mutAA==e$wtAA & !is.na(e$wtAA),"type"]<-"synonymous"
e[e$mutation=="wt","type"]<-"WT"
ave<-mean(e[e$type=="WT","score"])
e$score<-100*(e$score-1)/(ave-1)
write.csv(e,file = "/Users/KRONCKE/Box Sync/Kroncke_Lab/kcnh2/sequencing/VANTAGE/Tile1/Library1-DM278/DMS_DM278_5239-2_e.csv")


# PLOTTING
mat = e[!is.na(e$mutAA),c("res","mutAA","score","CountTotal","count")]
f<-aggregate(cbind(mat$score,mat$count), list(res=mat$res,Mut=mat$mutAA), median, na.rm =T)
colnames(f)[3:4]<-c("score","subcount")

f$score<-100*(f$score-1+norm)/(norm-1)
write.csv(f,file = "/Users/KRONCKE/Box Sync/Kroncke_Lab/kcnh2/sequencing/VANTAGE/Tile1/Library1-DM278/DMS_DM278_5239-2.csv")

f.m<-melt(f,c("res","Mut"),"score")
p<-ggplot(f.m, aes(ordered(Mut, levels = aas),ordered(res, levels = rev(mixedsort(unique(res))))))+
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradientn(colours = c("#C65911","#FFC000","#FFFFA8", "#FFFFFF","#BDD7EE"), values = c(0,0.25,0.5,0.70,1), na.value="grey50", guide = "colourbar") 
base_size<-9
p+ theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none", axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) + 
  scale_x_discrete(position = "top")




#barcode counts
f.m<-melt(f,c("res","Mut"),"subcount")
f.m[f.m$value==0,"value"]<-NA
f.m$value<-log10((f.m$value))
p<-ggplot(f.m, aes(ordered(Mut, levels = aas),ordered(res, levels = rev(mixedsort(unique(res))))))+ geom_tile(aes(fill = value), color = "white") + scale_fill_gradient(low = "white", high = "steelblue") 
base_size<-9
p+ theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "none", axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) + scale_x_discrete(position = "top")



