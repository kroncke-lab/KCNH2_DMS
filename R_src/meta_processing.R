


library(reshape2)
library(ggplot2)
library(gtools)
library(stringr)
library(dplyr)
library(data.table)

f.out <- function(dt){
  #dt<-dt[dt$CountAll>100,]
  #dt<-dt[dt$CountHigh/dt$CountAll<1 & dt$CountMed/dt$CountAll<1 & dt$CountLow/dt$CountAll<1 & dt$CountNeg/dt$CountAll<1,]
  mat = dt[!is.na(dt$mutAA),c("res","mutAA","score","CountTotal","count","resnum")]
  f<-aggregate(cbind(mat$score,mat$count), 
               list(res=mat$res,Mut=mat$mutAA,resnum=mat$resnum), median, na.rm =T)
  colnames(f)[3:5]<-c("resnum", "score","subcount")
  return(f)
}

setwd("C:/Users/KRONCKE/Box Sync/Kroncke_Lab/kcnh2/sequencing/VANTAGE/Tile1/")
a<-read.csv("Library1-DM278/DMS_DM278_5239-2_e.csv", stringsAsFactors = F)
a<-a[!is.na(a$mutation),]
a.f<-f.out(a)

b<-read.csv("Library1-DM278/DMS_DM278_5239-1_e.csv", stringsAsFactors = F)
b<-b[!is.na(b$mutation),]
b.f<-f.out(b)

c<-read.csv("Library1-DM278/DMS_DM278_5125_e.csv", stringsAsFactors = F)
c<-c[!is.na(c$mutation),]
c.f<-f.out(c)

e<-read.csv("Library2-LV310/DMS_LV310_5282_e.csv", stringsAsFactors = F)
e<-e[!is.na(e$mutation),]
e.f<-f.out(e)

g<-read.csv("Library2-LV310/DMS_LV310_5635_e.csv", stringsAsFactors = F)
g<-g[!is.na(g$mutation),]
g.f<-f.out(g)

h<-read.csv("Library2-LV310/DMS_LV310_5572_e.csv", stringsAsFactors = F)
h<-h[!is.na(h$mutation),]
h.f<-f.out(h)

m<-read.csv("Library3-DM344/DMS_wCodons_e-RU-6326.csv", stringsAsFactors = F)
m<-m[!is.na(m$mutation),]
m<-m[m$resnum>78 & m$resnum<86,]
m.f<-f.out(m)

n<-read.csv("Library2-LV310/DMS_LV310_5635_2_e.csv", stringsAsFactors = F)
n<-n[!is.na(n$mutation),]
n.f<-f.out(n)

o<-read.csv("Library2-LV310/DMS_LV310_5572_2_e.csv", stringsAsFactors = F)
o<-o[!is.na(o$mutation),]
o.f<-f.out(o)

p<-read.csv("Library2-LV310/DMS_LV310_6114_e.csv", stringsAsFactors = F)
p<-p[!is.na(p$mutation),]
p.f<-f.out(p)

q<-read.csv("Library2-LV310/DMS_LV310_6210_e.csv", stringsAsFactors = F)
q<-q[!is.na(q$mutation),]
q.f<-f.out(q)

r<-read.csv("Library2-LV310/DMS_LV310_6276-1_e.csv", stringsAsFactors = F)
r<-r[!is.na(r$mutation),]
r.f<-f.out(r)

s<-read.csv("Library2-LV310/DMS_LV310_6276-2_e.csv", stringsAsFactors = F)
s<-s[!is.na(s$mutation),]
s.f<-f.out(s)

# removed n.f and g.f
i<-rbind(a.f,b.f,c.f,e.f,h.f,m.f,o.f,p.f,r.f,s.f)#,q.f)

#i$mutation<-as.character(i$mutation)
#i$res<-as.character(i$res)
i<-i[i$resnum>25 & i$resnum<103,]
#i<-i[!is.na(i$mutation),]

#i<-i[i$ i$CountAll>199,]

mat = i[!is.na(i$Mut),c("res","Mut","score","subcount")]
f<-aggregate(cbind(mat$score,mat$subcount), list(res=mat$res,Mut=mat$Mut), mean, na.rm =T)
colnames(f)[3:4]<-c("score","subcount")

# testing ways to remove outliers
medians<-f[,c(1,2,3)]
names(medians)<-c("res","Mut","median")
q<-merge(i,medians,all = T)

q$residual<-(q$median-q$score)
q$stdev<-q$residual/sd(q$residual)
i<-q[abs(q$stdev)<2.5,]

aas<-c("A","V","I","L","M","F","Y","W","P","G","C","S","T","N","Q","H","R","K","D","E","X")
x<-data.frame(0,0)
k<-0
colnames(x)<-c("res","Mut")
for (r in c(unique(i$res), "79 A", "80 A", "81 Q" )){
  for (aa in aas){
    x[k,c("Mut","res")]<-c(aa,r)
    k=k+1
  }
}

i<-bind_rows(i,x)
mat = i[!is.na(i$Mut),c("res","Mut","score","subcount")]
f<-aggregate(cbind(mat$score,mat$subcount), list(res=mat$res,Mut=mat$Mut), median, na.rm =T)
colnames(f)[3:4]<-c("score","subcount")

d<-f

t<-str_split_fixed(d$res, "[0-9]+ ",2)
d$wt<-t[,2]

j<-str_split_fixed(d$res,"[A-Z]",3)
d$resnum<-as.integer(j[,1])
#d<-d[d$resnum>25 & d$resnum<104,]

d$type<-"missense"
d[d$Mut=="X","type"]<-"nonsense"
d[d$Mut==d$wt,"type"]<-"synonymous"

d$type<-factor(d$type, levels = c( "synonymous","nonsense", "missense"))

ave<-median(d[!is.na(d$score) & d$type=="synonymous","score"], na.rm = T)
d$score<-100*d$score/ave

p <- ggplot(d, aes(x=type,y=score)) + 
  geom_violin(scale = "width")
p + geom_boxplot(width=0.1) + geom_jitter(shape=16, position=position_jitter(0.1))

### HEAT MAP ###
#f.dt <- data.table(d)
#d<-melt(f,c("res","Mut"),"score")
d$category<-NA
d[!is.na(d$score) & d$score<25,"category"]<-"Dark"
d[!is.na(d$score) & d$score>=25 & d$score<50,"category"]<-"DarkLight"
d[!is.na(d$score) & d$score>=50 & d$score<75,"category"]<-"Medium"
d[!is.na(d$score) & d$score>=75 & d$score<150,"category"]<-"White"
d[!is.na(d$score) & d$score>=150,"category"]<-"Blue"
d$category<-factor(d$category, levels=c("Dark","DarkLight","Medium","White","Blue"))

p<-ggplot(d, aes(ordered(Mut, levels = aas),ordered(res, levels = rev(mixedsort(unique(res))))))+
  geom_tile(aes(fill = factor(category)), color = "white") +
  scale_colour_manual(values = c("#C65911","#FFC000","#FFFFA8","#FFFFFF","#BDD7EE"), 
                      aesthetics = c('color','fill'),
                      na.value = "grey50") 
p




# Save output
write.csv(d[,c("resnum", "wt", "res","Mut","score", "type")],"C:/Users/KRONCKE/Dropbox/sat paper/DMS_all.csv", row.names = F)
write.csv(d[,c("resnum", "wt", "res","Mut","score", "type")],"DMS_all.csv", row.names = F)


