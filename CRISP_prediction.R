############## Prediction of eGFR Decline in Autosomal Dominant Polycystic Kidney Disease (ADPDK) ##############
############## Jiawei Meng
############## Mar 2020
############## R version 3.6.3

# load library
library(data.table)
library(zoo)
library(broom)
library(ggplot2)
# for animate the plot to present one line at a time
library(gganimate)
library(dplyr)
library(magrittr)
library(gifski)
library(png)
# to add time sequence
library(splitstackshape)
# for lasagnar plot
library(devtools)                                                                                                                                
library(fields)
library(lasagnar)   
library(RColorBrewer)
library(colorspace)
library(gplots)
# for PCA and machine learning
#load package
library(e1071)
library(dplyr)
library(mice)
library(pROC)
library(tidyverse)
library(glmnet)
library(factoextra)
library(caret)
library(ranger)
library(nnet)
library(e1071)
library(SuperLearner)
library(MASS)   
library(MLmetrics) # for MAPE
library(car)
library(caret)
library(randomForest)
library(ranger)
library("tidyr")
# for sigmoidal fit
library(minpack.lm)
library(drc)
library(lattice)
library(table1)

######################## Table 1 #######################
#####read data
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","basedate","dvdate","birthdate",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
summary(unique(raw1$visc))
#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
#don't know why delete "trunc_grp"
raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid")
fortableone<-d[,head(.SD,1),by=pkdid]  # records Y0 and Y1
fortableone$gender <- factor(fortableone$gender ,levels=c(1,2), labels=c("male","female"))
fortableone$race4 <- factor(fortableone$race4 ,levels=c(1,2,3,4), labels=c("Caucasian", "African American", "Hispanic", "Asian"))
table1( ~ age + educ + gender + race4+ vis, data = fortableone)
##########################################################


############################# visualization ##############################
########### Plot the original data eGFR to years for each pkdid###########
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc","birthdate",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]



#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
str(raw1)
#factor the pkdid
raw1$pkdid <- as.factor(raw1$pkdid)
#convert data to data table
d <- data.table(raw1, key="pkdid")
time<-ts(d,frequency = 1, start = s)

#the base date
s<-as.Date(d$basedate,format = "%m/%d/%Y")
#the visit data
e<-as.Date(d$dvdate,format = "%m/%d/%Y")
#convert the basedate and visitdate to the m/d/y format
d$basedate<-s
d$dvdate<-e

#to plot, make the datafame contain the visit data, eGFR and group
#x axis in the plot, the visit date
x <- d$dvdate
#y axis in the plot, eGFR
y <- d$ckd_epi
#check the eGFR missing value
which(is.na(y))
#remove the patients's date and eGFR when eGFR is missing
x1 <- x[-which(is.na(y))]   
y1<-y[-which(is.na(y))]  
#add the pkdid indicate the group
z<-d$pkdid
z1<-z[-which(is.na(y))]  
#combine visit date, eGFR and pkdid and convert to data frame
dd<-cbind(x1,y1,z1)
ddd<-as.data.frame(dd)
ddd$x1<-as.Date(x1)

### Plot the original data eGFR over the actual years for each patient 
# plot the original data eGFR over the actual years for each patient and add small point
plot_egfr<-ggplot(data=ddd, aes(x1, y1, group=z1))+scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month")
plot_egfr_jitter <- plot_egfr + geom_point(size=0.09,color = z1 ,position = "jitter")+labs(x = "years",y = "eGFR")
plot_egfr_addline <- plot_egfr_jitter + geom_line(aes(group = z1), color=z1, alpha = 0.5) #add lines



###########  make plot animate and present one line at a time ########### 
#to generate animate pdf. get library
devtools::install_github('thomasp85/gganimate')
setwd("/Users/mengjiawei/Desktop/animation")
as.Date(x1)

#one line add to previous one and repeat
#anim2<- plot_egfr_addline + transition_reveal(z1) 
#anim2 + enter_fade() + exit_shrink()

#present one line each time
anim2<- plot_egfr_addline + labs(title = "PKDID: {current_frame}")+
  transition_manual(z1,cumulative = FALSE) 
animate(anim2 + enter_fade() + exit_shrink() , renderer = av_renderer('animation.mp4'))
anim_save(filename = 'plot.mpg',
          path = '/Users/mengjiawei/Desktop/animation')


#################### create a lasagna plot to visualization #####################
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
#don't know why delete "trunc_grp"
raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid")

#transform to date type
d$dvdate <- as.Date(d$dvdate,format = "%m/%d/%Y")
d$basedate <- as.Date(d$basedate,format = "%m/%d/%Y")
d$birthdate <- as.Date(d$birthdate,format = "%m/%d/%Y")

#calculate the time sequence for each patient
dd2<-getanID(data = d, id.vars = "pkdid")
d$yrs <- dd2$.id

#x is the age
x <- rep(NA,length(d$age))

#y is the ckd_epi for each visit for each patient
y <- d$ckd_epi

#z is the patient group number 1 is for first patient..
z<-d$pkdid

#q<-d$basedate
w<- d$dvdate
e<- d$birthdate
h<- d$yrs

#conbine these three to a data frame
dd<-cbind(x,y,z,w,e,h)
ddd<-as.data.frame(dd)
ddd$e<-as.Date(e)
ddd$w<-as.Date(w)

#add for missing birthdate using same group people exist birthdate
for (s in 1:max(ddd$z)){
  ddd$e[ddd$z == s][which(is.na(ddd$e[ddd$z == s]))] <- (ddd$e[ddd$z == s])[1]
}

#add age using birthdate and dvdate
for (f in 1:max(ddd$z)){
  diff <- round(as.numeric(difftime(ddd$w[ddd$z == f], ddd$e[ddd$z==f]))/365.25,1)
  ddd$x[ddd$z == f] <- c(diff)
}

#remove the age na row in data (vecause that has no dvdata)
ddd1<-ddd[-which(is.na(ddd$x)),]

#remove the same age in one group( wrong data)
for (f in 1:max(ddd1$z)){
  ddd1$x[ddd1$z == f][duplicated(ddd1$x[ddd1$z == f],fromLast=TRUE)] <- NA
}
ddd2<-ddd1[-which(is.na(ddd1$x)),]

#create the H matrix for these data to plot
H.mat <- matrix(NA, nrow=max(ddd$z), ncol=max(ddd$h))
j<-c(sort(unique(ddd$h)))

for (n in 1:length(sort(unique(ddd$h)))){
  for (i in (1:max(ddd$z))){
    if (j[n] %in% (ddd$h[ddd$z == i])){
      H.mat[i,n] <- ddd$y[ddd$h == j[n] & ddd$z == i]
    }
  }
}
# add rownames and colnames for H matrix
rownames(H.mat)<-c(1:max(ddd$z))
colnames(H.mat)<-c(sort(unique(ddd$h)))

#plot
lasagna(H.mat,col=palette,legend=TRUE, cex.axis=0.75)
title(xlab= 'time sequence for visit', ylab = 'subject', line = 2)


########################### cluster heatmap  ############################

rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
#don't know why delete "trunc_grp"
raw1
raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid")

d$yr<-d$vis/10
unique(d$yr)

#transform to date type
d$dvdate <- as.Date(d$dvdate,format = "%m/%d/%Y")
d$basedate <- as.Date(d$basedate,format = "%m/%d/%Y")
d$birthdate <- as.Date(d$birthdate,format = "%m/%d/%Y")


#calculate the time sequence for each patient
library(splitstackshape)
dd2<-getanID(data = d, id.vars = "pkdid")
d$yrs <- dd2$.id
d

#x is the age
x <- rep(NA,length(d$age))
x
#y is the ckd_epi for each visit for each patient
y <- d$ckd_epi
y
#z is the patient group number 1 is for first patient..
z<-d$pkdid
z
#q<-d$basedate
w<- d$dvdate
e<- d$birthdate
h<- d$yr
j<- d$yrs
#conbine these three to a data frame
dd<-cbind(x,y,z,w,e,h,j)
dd
ddd<-as.data.frame(dd)
ddd
ddd$e<-as.Date(e)
ddd$w<-as.Date(w)
ddd
library(dplyr)
ddd<- filter(ddd, !(h %in% c("4.5", "6.5","7.5","8.5","9.5","9.6","9.9")))
ddd
#add for missing birthdate using same group people exist birthdate
for (s in 1:max(ddd$z)){
  ddd$e[ddd$z == s][which(is.na(ddd$e[ddd$z == s]))] <- (ddd$e[ddd$z == s])[1]
}
ddd$e
ddd

#add age using birthdate and dvdate
for (f in 1:max(ddd$z)){
  diff <- round(as.numeric(difftime(ddd$w[ddd$z == f], ddd$e[ddd$z==f]))/365.25,1)
  ddd$x[ddd$z == f] <- c(diff)
}
ddd$x

#remove the age na row in data (vecause that has no dvdata)
ddd1<-ddd[-which(is.na(ddd$x)),]
ddd1


ddd1
g<- ddd1[ddd1$z == 36,]
library(data.table)

#for (f in 1:max(ddd1$z)){
ddd2<-unique(setDT(ddd1[ddd1$z ==36,])[order(x, -y)], by = "x")
#}


for (f in 1:max(ddd1$z)){
  ddd1$x[ddd1$z == f][duplicated(ddd1$x[ddd1$z == f],fromLast=TRUE)] <- NA
}
ddd2<-ddd1[-which(is.na(ddd1$x)),]

ddd2<-ddd1[-which(is.na(ddd1$x)),]
library(devtools)                                                                                                                                
install_github("swihart/lasagnar")                                                               
library(fields)
library(lasagnar)   
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(colorspace)

#get the H matrix for these data to plot
H.mat <- matrix(NA, nrow=max(ddd$z), ncol=max(ddd$x))
for(i in 1:max(ddd$z)){
  H.mat[i,min(ddd$x):max(ddd$x)]
}
H.mat

library(gplots)
heatmap.2(H.mat)
rownames(H.mat)<-c(1:max(ddd$z))
colnames(H.mat)<-seq(ncol(H.mat))
library(RColorBrewer)
display.brewer.all()
palette <- brewer.pal(length(unique(c(x))), "Greys")[-2]
palette
#lasagna plot
lasagna(H.mat,col=palette,legend=TRUE, cex.axis=0.75)
title(xlab= 'time sequence for visit', ylab = 'subject', line = 2)


ddd2
length(unique(ddd2$h))

HH<-matrix(NA,nrow=max(ddd2$z),ncol=length(sort(unique(ddd2$h))))
dim(HH)
rownames(HH)<-c(1:max(ddd2$z))
colnames(HH)<-c(sort(unique(ddd2$h)))

ddd2h

k<-sort(unique(ddd$h))
length(k)
table(ddd$h)
j<-c(sort(unique(ddd$h)))
j


ddd$h

for (n in 1:length(sort(unique(ddd2$h)))){
  for (i in (1:242)){
    if (j[n] %in% (ddd2$h[ddd2$z == i])){
      HH[i,n] <- ddd2$y[ddd2$h == j[n] & ddd2$z == i]
    }
  }
}

HH
HH1 <- HH[,colSums(is.na(HH))<nrow(HH)]
HH1
dim(HH1)

library(devtools)                                                                                                                                
install_github("swihart/lasagnar")                                                               
library(fields)
library(lasagnar)   
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(colorspace)
library(gplots)

dist_no_na <- function(HH) {
  edist <- dist(HH)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

library(RColorBrewer)
my_palette <- colorRampPalette(c("grey0","grey50", "grey100"))(n = 299)
heatmap.2(HH, col=my_palette, scale = "none",trace="none",dendrogram="row", Colv="NA", na.color="white", distfun=dist_no_na, xlab = "visit",
          ylab = "participant number",keysize=0.75, key.par = list(cex=0.5))



############### plot the original data serum creatinine to actual years  ##############
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc","birthdate",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
#don't know why delete "trunc_grp"

#factor the pkdid
raw1$pkdid <- as.factor(raw1$pkdid)
#convert data to data table
d <- data.table(raw1, key="pkdid")

time<-ts(d,frequency = 1, start = s)

#the base date
s<-as.Date(d$basedate,format = "%m/%d/%Y")
#the visit data
e<-as.Date(d$dvdate,format = "%m/%d/%Y")
#convert the basedate and visitdate to the m/d/y format
d$basedate<-s
d$dvdate<-e

#to plot, make the datafame contain the visit data, eGFR and group
#x axis in the plot, the visit date
x <- d$dvdate
#y axis in the plot, eGFR
y <- d$serumcreat
#check the eGFR missing value
which(is.na(y))
#remove the patients's date and eGFR when eGFR is missing
x1 <- x[-which(is.na(y))]   
y1 <- y[-which(is.na(y))]  
#add the pkdid indicate the group
z<-d$pkdid
z1<-z[-which(is.na(y))]  
#combine visit date, eGFR and pkdid and convert to data frame
dd<-cbind(x1,y1,z1)
ddd<-as.data.frame(dd)
ddd$x1<-as.Date(x1)


#plot
plot_egfr <- ggplot(data=ddd, aes(x1, y1, group=z1))+scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month")
plot_egfr_jitter <- plot_egfr + geom_point(size=0.09,color = z1 ,position = "jitter")+labs(x = "years",y = "eGFR")
plot_egfr_addline <- plot_egfr_jitter + geom_line(aes(group = z1), color=z1, alpha = 0.5) #add lines


#check minimum serum creatinine
min(y1) #0
#check maximum serum creatinine
max(y1) #2.2
#need to deal with serum creatnine=0
ddd$y1==0
#find the serum creatinine = 0 is z1=38
ddd[ddd$z1==38,]
ddd$y1
#replace the 0 with the near two value's mean:0.92
cc<-ddd[ddd$z1==38,]$y1
cc2<-replace(cc, cc==0, 0.92)
#recheck now
ddd[ddd$z1==38,]$y1<-cc2
ddd[ddd$z1==38,]$y1
min(ddd$y1)
#replot again
plot_egfr<-ggplot(data=ddd, aes(x1, y1, group=z1))+scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month")
plot_egfr_jitter <- plot_egfr + geom_point(size=0.09,color = z1 ,position = "jitter")+labs(x = "years",y = "serum creatinine")
plot_egfr_addline <- plot_egfr_jitter + geom_line(aes(group = z1), color=z1, alpha = 0.5) #add lines


#hard to see the pattern so try log transform on serum creatinine
plot_egfr <- ggplot(data=ddd, aes(x1, log2(y1), group=z1))+scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month")
plot_egfr_jitter <- plot_egfr + geom_point(size=0.09,color = z1 ,position = "jitter")+labs(x = "years",y = "log2(serum creatinine)")
plot_egfr_addline <- plot_egfr_jitter + geom_line(aes(group = z1), color=z1, alpha = 0.5) #add lines


################## plot eGFR versus serum creatinine according to the formula ################
creat<-seq(0.49, 8.6, length.out = 100)
egfrfemaleblack<- 175 * (creat)^(-1.154) * (55)^(-0.203) * (0.742) * (1.212)
egfrfemaleother<- 175 * (creat)^(-1.154) * (55)^(-0.203) * (0.742) 
egfrmaleblack<- 175 * (creat)^(-1.154) * (55)^(-0.203) * (1.212)
egfrmaleother<- 175 * (creat)^(-1.154) * (55)^(-0.203) 
egfrfemaleblack


plot(creat,egfrfemaleblack,xlab="serum creatinine", ylab="eGFR",xlim=range(0:9), ylim =range(0:180),col="black",type="l")
lines(creat,egfrfemaleother,col="red",type="l")
lines(creat,egfrmaleblack,col="green",type="l")
lines(creat,egfrmaleother,col="blue",type="l")
legend("topright", c("FemaleBlack","FemaleOther","MaleBlack","MaleOther"),  
       col=c("black", "red", "green", "blue"), lty=1, cex=0.6, pt.cex = 1.2, lwd=1.5)
dev.off()



################################# PCA analysis for CRISP data #################3###########
#####read data
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc","basedate","dvdate","birthdate",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
#don't know why delete "trunc_grp"
raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid")
fortableone<-d[,head(.SD,1),by=pkdid]  # records Y0 and Y1
fortableone$gender <- factor(fortableone$gender ,levels=c(1,2), labels=c("male","female"))
fortableone$race4 <- factor(fortableone$race4 ,levels=c(1,2,3,4), labels=c("Caucasian", "African American", "Hispanic", "Asian"))
table1( ~ age + educ + gender + race4, data = fortableone)



raw3 <- d[,head(.SD,2),by=pkdid]  # records Y0 and Y1
dim(raw3) #482*62
#add y1 eGFR as ckd_epi_c to the y0 data
raw3$ckd_epi_c <- raw3$ckd_epi
raw3$ckd_epi_c[1] <- 0
for (i in 1:nrow(raw3))
{
  raw3$ckd_epi_c[i] <- raw3$ckd_epi[i+1]
}                                       
dim(raw3) #482*63 add y1 eGFR to Y0 data
rawt3 <- setDT(raw3)[,.SD[1L],by=pkdid] 
dim(rawt3) #242*63 derive the y0 data

## Data imputation of missing values
# delete ids and variables with a lot of NA
sort(colSums(is.na(rawt3))) 
sort(rowSums(is.na(rawt3)))
rawt3 <- rawt3[,-c("sldl","uprotc_ca","uprote_ca","esode_cc")] 
rawt3 <- rawt3[-which(rowSums(is.na(rawt3))>10),]  

# imputation for variables with NA using MICE, variables without NA just keep
rawt3imp <- rawt3 %>% select_if(colSums(is.na(rawt3))>0)
rawt3keep <- rawt3 %>% select_if(colSums(is.na(rawt3))==0)
imp <- mice(rawt3imp,maxit=5)
rawt3f <- complete(imp,5)
rawt3ff <- cbind(rawt3f,rawt3keep) 
sort(colSums(is.na(rawt3ff)))  # checked no missing values, this is the final data of rawt3 for use
dim(rawt3ff) #235*59


#####pca analysis for remove genetype, gender, race4, age, educ and tkv for just doing PCA on biomarkers
rawt3fff <- rawt3ff[,!(names(rawt3ff) %in% c("pkdid","ckd_epi","ckd_epi_c","genetype","gender","race4","age","educ","tkv"))]
pca_fit <- prcomp(as.matrix(rawt3fff),scale =T)
biomark_PC <- pca_fit$x
# apply pca on biomarker data
loadings <- pca_fit$rotation %>% as.data.frame()
raw.var <- round(pca_fit$sdev^2/sum(pca_fit$sdev^2)*100) #percent explained variance 
cumsum(raw.var)
par(mfrow=c(1,1))
#sd of PC plot
plot(pca_fit$sdev[1:30],type="l",ylab="SD of PC",xlab="PC number",main="Scree Plot")
#cumulative SD of PC plot
plot(cumsum(raw.var)[1:30],type="l",ylab="cumulative SD of PC",xlab="PC number",main="Scree Plot")
# find out 27 PCs explained for 90%


#####10-fold cross validation to find the properly numbers of PCs to use
#use lasso
rsqqq1<- vector() # create vector for r square 
preddd1<-rep(NA,10) # create vector for predicted value
obs1<-rep(NA,10) # create vector for observed value
# loop for all the PC, using 10-fold cross validation and leave-one-out method to do cross validation to pick a PC which is small and have good R-squared
for (m in 1:50){
  bioPC <- biomark_PC[,c(1:m)]
  raw_df <- cbind(rawt3ff$ckd_epi,rawt3ff$ckd_epi_c,rawt3ff$genetype,rawt3ff$gender,rawt3ff$race4,rawt3ff$age,rawt3ff$educ,rawt3ff$tkv,
                  data.frame(bioPC));names(raw_df)[1]<- 'ckd_epi';
  names(raw_df)[2]<- 'ckd_epi_c';names(raw_df)[3]<-'genetype';names(raw_df)[4]<-'gender';names(raw_df)[5]<-'race4';names(raw_df)[6]<-'age';names(raw_df)[7]<-'educ';names(raw_df)[8]<-'tkv'
  n_folds<-10
  folds_i <- sample(rep(1:n_folds,length.out=nrow(raw_df)))
  for (k in 1:n_folds)
  {
    set.seed(123456) # set a seed for SuperLearner model
    test_i <- which(folds_i == k) # test_i: find i = k fold 
    train <- raw_df[-test_i, ] # train: the rest obs
    test <- raw_df[test_i, ] # test: the left out obs (folds_i)
    trainx <- train[,!(names(train) %in% c("pkdid","ckd_epi_c"))] # trainx: predictors for model building, "pkdid", "ckd_epi_c" should not be used
    trainy <- train[,c("ckd_epi_c")] # trainy: outcome for model building
    testx <- test[,!(names(test) %in% c("pkdid","ckd_epi_c"))] # testx:  predictors for the left out obs prediction
    testy <- test[,c("ckd_epi_c")] # trainy: outcome for model building
    sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
    preddd1[k]<-data.table(sl_lasso$SL.predict) ## the predicted outcome value of ckd_epi_c for each k fold, after the loop get a data table
    obs1[k]<-data.table(testy) # the data table for observed value of ckd_epi_c for each k fold
  }
  pred1<-unlist(preddd1[1:k]) # combine all the fold's predicted outcome value 
  obsy1<-unlist(obs1[1:k]) # combine all the fold's observed outcome value
  rsqqq1[m]<-cor(obsy1,pred1)^2  # R squared calculation for each number's of used PC

}
#plot to see the pattern for r squared and number of used PC 
y<-rsqqq1[1:m] # r squared
x<-1:m # number of PC are used 
#plot
plot(x,y,type="o",xlab="n_PCs",ylab="rsq",main="10-fold cross validation")
axis(side=1, at=seq(0,40,5))
#table the result
rlt<-data.frame(1:m, rsqqq1[1:m])

############## use first 12 PCs to prediction #############
#lasso
rsqu11 <- vector()
rsq11<- vector()
mse11<-vector()
preddd11<-rep(NA,10)
obs11<-rep(NA,10)
rsqu11<-rep(NA,10)
bioPC <- biomark_PC[,c(1:12)]
raw_df <- cbind(rawt3ff$ckd_epi,rawt3ff$ckd_epi_c,rawt3ff$genetype,rawt3ff$gender,rawt3ff$race4,rawt3ff$age,rawt3ff$educ,rawt3ff$tkv,
                data.frame(bioPC));names(raw_df)[1]<- 'ckd_epi';
names(raw_df)[2]<- 'ckd_epi_c';names(raw_df)[3]<-'genetype';names(raw_df)[4]<-'gender';names(raw_df)[5]<-'race4';names(raw_df)[6]<-'age';names(raw_df)[7]<-'educ';names(raw_df)[8]<-'tkv'
n_folds<-10
folds_i <- sample(rep(1:n_folds,length.out=nrow(raw_df)))
for (k in 1:n_folds)
{
  set.seed(123456)
  test_i <- which(folds_i == k)
  train <- raw_df[-test_i, ]
  test <- raw_df[test_i, ]
  trainx <- train[,!(names(train) %in% c("pkdid","ckd_epi_c"))]
  trainy <- train[,c("ckd_epi_c")]
  testx <- test[,!(names(test) %in% c("pkdid","ckd_epi_c"))]
  testy <- test[,c("ckd_epi_c")]
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet")
  rsqu11[k] <- cor(testy,sl_lasso$SL.predict)^2
  mse11[k] <- sum((testy-sl_lasso$SL.predict)^2) / length(testy)
  preddd11[k]<-data.table(sl_lasso$SL.predict)
  obs11[k]<-data.table(testy)
}
pred11<-unlist(preddd11[1:k])
obsy11<-unlist(obs11[1:k])
cor(obsy11,pred11)^2 # 0.6993328
sum((obsy11-pred11)^2) / length(obsy11)
#164.9528
MAE(pred11, obsy11) ## MAE calculation
#9.65814
MAPE(pred11, obsy11) #  mean absolute percentage error calculation
#0.1160113
testx

#random forest
rsqu3 <- vector()
rsq3<- vector()
mse3<-vector()
preddd3<-rep(NA,10)
obs3<-rep(NA,10)
rsqu3<-rep(NA,10)
bioPC <- biomark_PC[,c(1:12)]
raw_df <- cbind(rawt3ff$ckd_epi,rawt3ff$ckd_epi_c,rawt3ff$genetype,rawt3ff$gender,rawt3ff$race4,rawt3ff$age,rawt3ff$educ,rawt3ff$tkv,
                data.frame(bioPC));names(raw_df)[1]<- 'ckd_epi';
names(raw_df)[2]<- 'ckd_epi_c';names(raw_df)[3]<-'genetype';names(raw_df)[4]<-'gender';names(raw_df)[5]<-'race4';names(raw_df)[6]<-'age';names(raw_df)[7]<-'educ';names(raw_df)[8]<-'tkv'
n_folds<-10
folds_i <- sample(rep(1:n_folds,length.out=nrow(raw_df)))
for (k in 1:n_folds)
{
  set.seed(123456)
  test_i <- which(folds_i == k)
  train <- raw_df[-test_i, ]
  test <- raw_df[test_i, ]
  trainx <- train[,!(names(train) %in% c("pkdid","ckd_epi_c"))]
  trainy <- train[,c("ckd_epi_c")]
  testx <- test[,!(names(test) %in% c("pkdid","ckd_epi_c"))]
  testy <- test[,c("ckd_epi_c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest")
  rsqu3[k] <- cor(testy,sl_rf$SL.predict)^2
  mse3[k] <- sum((testy-sl_rf$SL.predict)^2) / length(testy)
  preddd3[k]<-data.table(sl_rf$SL.predict)
  obs3[k]<-data.table(testy)
}
pred3<-unlist(preddd3[1:k])
obsy3<-unlist(obs3[1:k])
cor(obsy3,pred3)^2 #0.7173774
sum((obsy3-pred3)^2) / length(obsy3)
#162.1528
MAE(pred3, obsy3) ## MAE calculation
#9.851923
MAPE(pred3, obsy3) #  mean absolute percentage error calculation
#0.1177857

#svm

rsqu4 <- vector()
rsq4<- vector()
mse4<-vector()
preddd4<-rep(NA,10)
obs4<-rep(NA,10)
rsqu4<-rep(NA,10)
bioPC <- biomark_PC[,c(1:12)]
raw_df <- cbind(rawt3ff$ckd_epi,rawt3ff$ckd_epi_c,rawt3ff$genetype,rawt3ff$gender,rawt3ff$race4,rawt3ff$age,rawt3ff$educ,rawt3ff$tkv,
                data.frame(bioPC));names(raw_df)[1]<- 'ckd_epi';
names(raw_df)[2]<- 'ckd_epi_c';names(raw_df)[3]<-'genetype';names(raw_df)[4]<-'gender';names(raw_df)[5]<-'race4';names(raw_df)[6]<-'age';names(raw_df)[7]<-'educ';names(raw_df)[8]<-'tkv'
n_folds<-10
folds_i <- sample(rep(1:n_folds,length.out=nrow(raw_df)))
raw_df$genetype <- nnet::class.ind(raw_df$genetype)  # Genetype need to be used as ind
raw_df$genetype <- raw_df$genetype[,2:4]
raw_df$genetype
for (k in 1:n_folds)
{
  set.seed(123456)
  test_i <- which(folds_i == k)
  train <- raw_df[-test_i, ]
  test <- raw_df[test_i, ]
  trainx <- train[,!(names(train) %in% c("pkdid","ckd_epi_c"))]
  trainy <- train[,c("ckd_epi_c")]
  testx <- test[,!(names(test) %in% c("pkdid","ckd_epi_c"))]
  testy <- test[,c("ckd_epi_c")]
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  rsqu4[k] <- cor(testy,sl_svm$SL.predict)^2
  mse4[k] <- sum((testy-sl_svm$SL.predict)^2) / length(testy)
  preddd4[k]<-data.table(sl_svm$SL.predict)
  obs4[k]<-data.table(testy)
}
pred4<-unlist(preddd4[1:k])
obsy4<-unlist(obs4[1:k])
cor(obsy4,pred4)^2 #0.6086237
sum((obsy4-pred4)^2) / length(obsy4)
#214.0377
MAE(pred4, obsy4) ## MAE calculation
#11.19481
MAPE(pred4, obsy4) #  mean absolute percentage error calculation
# 0.1380251


#linear regression
rsqu2 <- vector()
rsq2<- vector()
mse2<-vector()
preddd2<-rep(NA,10)
obs2<-rep(NA,10)
rsqu2<-rep(NA,10)
bioPC <- biomark_PC[,c(1:12)]
raw_df <- cbind(rawt3ff$ckd_epi,rawt3ff$ckd_epi_c,rawt3ff$genetype,rawt3ff$gender,rawt3ff$race4,rawt3ff$age,rawt3ff$educ,rawt3ff$tkv,
                data.frame(bioPC));names(raw_df)[1]<- 'ckd_epi';
names(raw_df)[2]<- 'ckd_epi_c';names(raw_df)[3]<-'genetype';names(raw_df)[4]<-'gender';names(raw_df)[5]<-'race4';names(raw_df)[6]<-'age';names(raw_df)[7]<-'educ';names(raw_df)[8]<-'tkv'
raw_df$genetype <- class.ind(raw_df$genetype)
raw_df$gender <- class.ind(raw_df$gender)
raw_df$race4 <- class.ind(raw_df$race4)
n_folds<-10
folds_i <- sample(rep(1:n_folds,length.out=nrow(raw_df)))
for (k in 1:n_folds)
{
  set.seed(123456)
  test_i <- which(folds_i == k)
  train <- raw_df[-test_i, ]
  test <- raw_df[test_i, ]
  trainx <- train[,!(names(train) %in% c("pkdid","ckd_epi_c"))]
  trainy <- train[,c("ckd_epi_c")]
  testx <- test[,!(names(test) %in% c("pkdid","ckd_epi_c"))]
  testy <- test[,c("ckd_epi_c")]
  sl_lm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.lm")
  #rsqu2[k] <- cor(testy,sl_lm$SL.predict)^2
  #mse2[k] <- sum((testy-sl_lm$SL.predict)^2) / length(testy)
  preddd2[k]<-data.table(sl_lm$SL.predict)
  obs2[k]<-data.table(testy)
}
pred2<-unlist(preddd2[1:k])
obsy2<-unlist(obs2[1:k])
(cor(obsy2,pred2))^2 #0.6833897
sum((obsy2-pred2)^2) / length(obsy2)
#173.6954
MAE(pred2, obsy2) ## MAE calculation
#10.12158
MAPE(pred2, obsy2) #  mean absolute percentage error calculation
# 0.1208312


# for lasso
par(mfrow=c(2,2))
op<- par(cex = 0.5)
plot(obsy11,pred11,lty="solid",col="black",cex=0.5, ylab="Predicted eGFR",xlab="Fitted eGFR",xlim=range(10:160),
     ylim=range(10:160),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.3,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy3,pred3,lty="solid",col="black",cex=0.5,
     ylab="Predicted eGFR",xlab="Fitted eGFR",xlim=range(10:160),ylim=range(10:160),cex.axis=1.5, cex.lab=1.5)
op <- par(cex = 0.5)
legend("topleft",legend=c("RandomForest"),
       col=c("black"),lty=1,cex=2.3,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy4,pred4,lty="solid",col="black",cex=0.5,
     ylab="Predicted eGFR",xlab="Fitted eGFR",xlim=range(10:160),ylim=range(10:160),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("black"),lty=1,cex=2.3,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for linear regression
plot(obsy2,pred2,lty="solid",col="black",cex=0.5,
     ylab="Predicted eGFR",xlab="Fitted eGFR",xlim=range(10:160),ylim=range(10:160),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("LinearRegression"),
       col=c("black"),lty=1,cex=2.3,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots



############ generate a matrix of R2 n year by n year for the pairwise correlation ############
#####read data
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")
#delete the data which will not be used for PCA 
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc","basedate","dvdate","birthdate",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
#delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 
#don't know why delete "trunc_grp"
raw1
raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid")
sort(colSums(is.na(d))) 
d <- d[,-c("sldl","uprotc_ca","uprote_ca","esode_cc")]
# esode_cc is deleted due to collinear
d <- d[-which(rowSums(is.na(d))>10),]   # 466*58
# imputation 
dimp <- d %>% select_if(colSums(is.na(d))>0)
dkeep <- d %>% select_if(colSums(is.na(d))==0)
imp <- mice(dimp,maxit=5)
rawf <- mice::complete(imp,5)
rawff <- cbind(rawf,dkeep) 
sort(colSums(is.na(rawff)))#466*58
#derive the pkdid and eGFR column
newraw <-rawff[, c("pkdid","ckd_epi")]
#get the year column
newraw2<-getanID(data = newraw, id.vars = "pkdid")
newraw$yrs <- newraw2$.id
newraw
sort(colSums(is.na(newraw)))
newrawspr<-spread(newraw,yrs,ckd_epi)
#calculate the pairewise correlation
corr<-cor(newrawspr[sapply(newrawspr, is.numeric)], use='pairwise')
rsquare<-corr^2
rsquare
# Note:
#1         2         3         4
#1 1.0000000 0.6898540 0.7273942 0.6671557
#2 0.6898540 1.0000000 0.8085315 0.7524516
#3 0.7273942 0.8085315 1.0000000 0.7932669
#4 0.6671557 0.7524516 0.7932669 1.0000000
#The R squared value (R2) seems good enough, but for the prediction power the eGFR have on themselves, it’s pretty closed to the machine learning R squared value (R2)


############################# sigmoidal fitting #######################
################ Using clinical information predict inflection point parameter (c) 
######## original data cleaning ###########
######## import data 
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata_Xinhang.csv")

######## delete the data which will not be used ########
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp"))]
# delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 

raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid") # transform to data table

######## transform the date data to date type ########
d$dvdate <- as.Date(d$dvdate,format = "%m/%d/%Y")
d$basedate <- as.Date(d$basedate,format = "%m/%d/%Y")
d$birthdate <- as.Date(d$birthdate,format = "%m/%d/%Y")


######## calculate the visit time sequence for each patient ########
dd2<-getanID(data = d, id.vars = "pkdid")
d$yrs <- dd2$.id


######## calculating age ########
# x is the age
x <- rep(NA,length(d$age))

# y is the ckd_epi for each visit for each patient
y <- d$ckd_epi

# z is the patient group (number 1 is for first patient..)
z<- d$pkdid

# w is the visit date
w<- d$dvdate

# e is the birthdate for each patient
e<- d$birthdate

# h is the visit time sequence
h<- d$yrs

# conbine these to a data frame dd
dd<-cbind(x,y,z,w,e,h)
dd<-as.data.frame(dd)


# add this 'dd' data frame to original data frame 'd'
newd<- cbind(dd,d)
ddd<-as.data.frame(newd)
ddd$e<-as.Date(e)
ddd$w<-as.Date(w)


# fill missing birthdate using same birthdate from the same person
for (s in 1:max(ddd$z)){
  ddd$e[ddd$z == s][which(is.na(ddd$e[ddd$z == s]))] <- (ddd$e[ddd$z == s])[1]
}

# get each age data for patients using birthdate and dvdate
for (f in 1:max(ddd$z)){
  diff <- round(as.numeric(difftime(ddd$w[ddd$z == f], ddd$e[ddd$z==f]))/365.25,1)
  ddd$x[ddd$z == f] <- c(diff)
}

# remove the missing age row in data (because that one doesn't have dvdata)
ddd1<-ddd[-which(is.na(ddd$x)),]

# remove junk
for (f in 1:max(ddd1$z)){
  ddd1$x[ddd1$z == f][duplicated(ddd1$x[ddd1$z == f],fromLast=TRUE)] <- NA
}
ddd2<-ddd1[-which(is.na(ddd1$x)),]
ddd2

# remove uncompleted age column from original clinical information data frame, keep the completed age
ddd3 <- subset(ddd2, select = -age)
ddd3 #ddd3 is the clinical information for each patient which has no missing age data


##################################### nlsLM function in minpack #####################################
######## write a function for plotting ########
m3 = function(data){
  x = data$x
  y = data$y
  
  # find the best starting point for c
  c.start = x[which(abs(y-60)==min(abs(y-60),na.rm=T))]
  
  # fit model
  fitmodel <- nlsLM(y ~(5+((110-5)/(1 + (x/c)^b))),start=c(b=0,c=c.start))
  b = round(summary(fitmodel)$parameters[1,1],2)
  c = round(summary(fitmodel)$parameters[2,1],1)
  plot(x, y, ylim=c(0,130), xlab="Age", ylab="eGFR", main=paste("SubjectID = ", eval(i)))
  legend("bottomleft", c( paste("slope = ", eval(b)),paste("Age_50% = ", eval(c)) ), bty="n")
  lines(x[!is.na(y)], fitted(fitmodel),col = 2, lwd = 2)  
}


######## make some plots ########
par(mfrow=c(3,3))
for (i in 41:49){
  m3(ddd3[which(ddd3$z == i),])
}
## Notes:  a  small portion still do not fit well  due to flat eGFR over time.


######## fit on everyone and get c,b and p values ######## 
n = max(ddd3$z)
bcEst = matrix(NA,n, 5)
names(bcEst) =c("b","pval_b","c","pval_c", "cor_coef")

for(i in (1:n)){
  print(i) # monitor the progress
  data = ddd3[which(ddd3$z == i),]
  x = data$x
  y = data$y
  c.start = x[which(abs(y-60)==min(abs(y-60),na.rm=T))]
  
  # fit model (Errors are generated for some subjects.  Use tryCatch to skip them)
  fitmodel <- tryCatch(nlsLM(y ~(5+((110-5)/(1 + (x/c)^b))),start=c(b=0,c=c.start)),error = function(e) NULL)
  # get the results for the subjects that works for fitting the model
  if (!is.null(fitmodel)){
    # get the output
    bcEst[i,1] = summary(fitmodel)$parameters[1,1]
    bcEst[i,2] = summary(fitmodel)$parameters[1,4]
    bcEst[i,3] = summary(fitmodel)$parameters[2,1]
    bcEst[i,4] = summary(fitmodel)$parameters[2,4]
    bcEst[i,5] = cor(y[!is.na(y)], predict(fitmodel)) 
  }
} 

bcEst = matrix(NA,n, 6)
names(bcEst) =c("b","pval_b","c","pval_c", "cor_coef", "ESRD_age")

for(i in (1:n)){
  print(i) # monitor the progress
  data = ddd3[which(ddd3$z == i),]
  x = data$x
  y = data$y
  c.start = x[which(abs(y-60)==min(abs(y-60),na.rm=T))]
  # fit model (Errors are generated for some subjects.  Use tryCatch to skip them)
  fitmodel <- tryCatch(nlsLM(y ~(5+((110-5)/(1 + (x/c)^b))),start=c(b=0,c=c.start)),error = function(e) NULL)
  # get the results for the subjects that works for fitting the model
  if (!is.null(fitmodel)){
    # get the output
    bcEst[i,1] = summary(fitmodel)$parameters[1,1]
    bcEst[i,2] = summary(fitmodel)$parameters[1,4]
    bcEst[i,3] = summary(fitmodel)$parameters[2,1]
    bcEst[i,4] = summary(fitmodel)$parameters[2,4]
    bcEst[i,5] = cor(y[!is.na(y)], predict(fitmodel)) 
    b= summary(fitmodel)$parameters[1,1]
    c= summary(fitmodel)$parameters[2,1]
    bcEst[i,6] = c*((110-15)/(15-5))^(1/b)
  }
} 

# Note: For some reason, the output matrix has some extra stuff at the end.  It is likely related to the matrix defination.  Try to change it to data frame in the future.

# round the coef and pval outputs
bcEstimate = cbind(round(bcEst[1:242,1],2 ), round(bcEst[1:242,2],3 ),round(bcEst[1:242,3],1 ), round(bcEst[1:242,4],3 ),round(bcEst[1:242,5],3 ),round(bcEst[1:242,6], 1 ) ) 
colnames(bcEstimate) =c("b","pval_b","c","pval_c", "cor_coef", "ESRD_age")

# remove the extra stuff
BCestimate = bcEstimate[1:242,]
write.csv(BCestimate, file="/Users/mengjiawei/Documents/CRISP-PKD/DATA/OutputData/sigmoidalParameterEst.csv") 


### Look at the R2 of the subjects with successful fit
idx = (which(BCestimate[,4]<0.05)) #176 (73%)
par(mfrow=c(1,1))
hist(BCestimate[idx,5])
hist(BCestimate[idx,5]^2, main="Histogram of R square",
     xlab="R square of subjects with successful fit") # 116 with R2 >0.7

# note: will use these successful fit data for below analysis

############################# machine learning to predict c #####################
################### combine estimated c from fit and clinical info 
######## remove the birthdate, dvdate etc which are used to calculate missing age ########
colnames(ddd3)[which(names(ddd3) == "x")] <- "age"
ddd3 <- subset(ddd2, select = -age)
ddd3 <- subset(ddd3, select = -y)  
ddd3 <- subset(ddd3, select = -w)
ddd3 <- subset(ddd3, select = -e)
ddd3 <- subset(ddd3, select = -h)

######## import estimated 'c' data ########
rawdata2 <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/OutputData/sigmoidalParameterEst.csv")
rawdata2
######## add c and pvalue for c to the clinical information data ########
ddd3$c<-NA
ddd3$cpval<-NA
for(i in (1:nrow(rawdata2))){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
}

ddd3[ddd3$z == 215,]
rawdata2[215,]
# Note: find out z=215 is empty

# remove the no.215
ddd3<-ddd3[-which((ddd$z == 215)),]
ddd3


######add c and pvalue for c for each patient
for(i in (1:214)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
for(i in (216:242)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
ddd4<-ddd3[-which(is.na(ddd3$c)),] # remove the NA in exist c
ddd5 <- ddd4 %>% filter(cpval <= 0.05) # keep the c which p_value < 0.05


# delete ids and variables with a lot of NA
ddd5<-data.table(ddd5, key="z")
sort(colSums(is.na(ddd5))) 
sort(rowSums(is.na(ddd5)))
ddd5 <- ddd5[,-c("sldl","uprotc_ca","uprote_ca","esode_cc")] 
ddd5 <- ddd5[-which(rowSums(is.na(ddd5))>10),]  


# imputation for variables with NA using MICE, variables without NA just keep
ddd5imp <- ddd5 %>% select_if(colSums(is.na(ddd5))>0)
ddd5keep <- ddd5 %>% select_if(colSums(is.na(ddd5))==0)
imp <- mice(ddd5imp,maxit=5)
ddd55 <- complete(imp,5)
ddd555 <- cbind(ddd55,ddd5keep) 
sort(colSums(is.na(ddd555)))  # checked no missing values, this is the final data for use

#change the name of column for age
colnames(ddd555)[which(names(ddd555) == "x")] <- "age"


###################################### machine learning prediction ###############################
######################## prediction 1):using everything else + visit 1 biomarker #################
# select visit 1 biomarkers
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd66 <- ddd6[,head(.SD,1),by=pkdid]  
# ddd66 is the visit 1 biomarkers dataset
ddd66
# convert to data frame and use this finalized data for prediction
final1 <- as.data.frame(ddd66)
# checked no missing values, this is the final data for use
sort(rowSums(is.na(final1)))
sort(colSums(is.na(final1))) 
dim(final1)

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy <- final1[,"c"] # observed true outcome value
mse <- sum((obsy-testy_pre)^2) / nrow(final1)  ## MSE calculation
mse # 46.50625
rsq <- (cor(obsy,testy_pre))^2 ## R squared calculation
rsq # 0.5648081
MAE(testy_pre, obsy) ## MAE calculation
# 4.708493
MAPE(testy_pre, obsy) #  mean absolute percentage error calculation
#0.1077987


# Method = Random Forest to predict c
testy_pre1.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre1.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy1.2 <- final1[,"c"] # observed true outcome value
mse1.2 <- sum((obsy1.2-testy_pre1.2)^2) / nrow(final1) ## MSE calculation
mse1.2 # 59.8971
rsq1.2 <- (cor(obsy1.2,testy_pre1.2))^2 ## R squared calculation
rsq1.2 # 0.4693433
MAE(testy_pre1.2, obsy1.2) ## MAE calculation
#5.45751
MAPE(testy_pre1.2, obsy1.2) #  mean absolute percentage error calculation
#0.1335129

# Method = SVM to predict c
final1$genetype <- nnet::class.ind(final1$genetype)  # Genetype need to be used as ind
final1$genetype <- final1$genetype[,2:4]
final1$genetype
testy_pre1.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre1.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy1.3 <- final1[,"c"] # observed true outcome value
mse1.3 <- sum((obsy1.3-testy_pre1.3)^2) / nrow(final1)  ## MSE calculation
mse1.3 #62.93364
rsq1.3 <- (cor(obsy1.3,testy_pre1.3))^2 ## R squared calculation
rsq1.3 #0.4439044
MAE(testy_pre1.3, obsy1.3) ## MAE calculation
#5.61435
MAPE(testy_pre1.3, obsy1.3) #  mean absolute percentage error calculation
#0.1330551

#####plot the c values (predicted from machine learning model vs predicted from sigmoidal curve fitting). 
# for lasso
par(mfrow=c(2,2))
op<- par(cex=0.5)
plot(obsy,testy_pre,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy1.2,testy_pre1.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy1.3,testy_pre1.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

######################## prediction 2):using everything else + visit 2 biomarker #################
#derive visit 2 biomarkers
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd77 <- ddd6[,.SD[2L],by=pkdid]  
# ddd77 is the visit 2 biomarkers dataset

final2<-as.data.frame(ddd77)
# checked no missing values, this is the final data for use
sort(rowSums(is.na(final2)))
sort(colSums(is.na(final2))) 
# note: have two missing rows for biomarkers

final2 <- na.omit(final2) ## delete the missing values

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre2<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final2$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final2[-i, ] # train: the rest obs
  test <- final2[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre2[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy2 <- final2[,"c"] # observed true outcome value
mse2 <- sum((obsy2-testy_pre2)^2) / nrow(final2)  ## MSE calculation
mse2 # 51.76693
rsq2 <- (cor(obsy2,testy_pre2))^2 ## R squared calculation
rsq2 # 0.5200685
MAE(testy_pre2, obsy2) ## MAE calculation
#5.137601
MAPE(testy_pre2, obsy2) #  mean absolute percentage error calculation
#0.1184768


# Method = Random Forest to predict c
testy_pre2.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final2$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- final2[-i, ] # train: the rest obs
  test <- final2[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre2.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy2.2 <- final2[,"c"] # observed true outcome value
mse2.2 <- sum((obsy2.2-testy_pre2.2)^2) / nrow(final2) ## MSE calculation
mse2.2 # 61.98482
rsq2.2 <- (cor(obsy2.2,testy_pre2.2))^2 ## R squared calculation
rsq2.2 # 0.4495015
MAE(testy_pre2.2, obsy2.2) ## MAE calculation
# 5.5667
MAPE(testy_pre2.2, obsy2.2) #  mean absolute percentage error calculation
#0.1354317


# Method = SVM to predict c
final2$genetype <- nnet::class.ind(final2$genetype)  # Genetype need to be used as ind
final2$genetype <- final2$genetype[,2:4]
final2$genetype
testy_pre2.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final2$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final2[-i, ] # train: the rest obs
  test <- final2[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre2.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy2.3 <- final2[,"c"] # observed true outcome value
mse2.3 <- sum((obsy2.3-testy_pre2.3)^2) / nrow(final2)  ## MSE calculation
mse2.3 #  62.84017
rsq2.3 <- (cor(obsy2.3,testy_pre2.3))^2 ## R squared calculation
rsq2.3 # 0.4594261
MAE(testy_pre2.3, obsy2.3) ## MAE calculation
#5.636193
MAPE(testy_pre2.3, obsy2.3) #  mean absolute percentage error calculation
#0.1326531


##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy,testy_pre,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy2,testy_pre2,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy1.2,testy_pre1.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy2.2,testy_pre2.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy1.3,testy_pre1.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy2.3,testy_pre2.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots





# Note: we found that compared with other c, when c > 60, the predicted value is far from the fitted value.
# Note: We want to plot the sigmoidal when c >60 and check these plots to see if they have some characteristics. 
#       We may also want to check on the starting point for these c > 60. Also, we will plot some c = 50, c=40 and c=30 to see whether the c > 60 plots have some problems with fitting.



########## Plot c > 60, check on fit, check on starting point, and plot some c = 50, c=40, c=30 ############
###################plot c > 60 #######################
id2 = (which(BCestimate[,3]>60 & BCestimate[,4]<0.05)) 
id2
par(mfrow=c(3,4))
for (i in id2){
  m3(ddd3[which(ddd3$z == i),])
}

###################plot c ~= 50 #######################
id50 = (which(BCestimate[,3] > 50 & BCestimate[,3] < 51 & BCestimate[,4]<0.05)) 
id50
par(mfrow=c(3,4))
for (i in id50){
  m3(ddd3[which(ddd3$z == i),])
}

###################plot c ~= 40 #######################
id40 = (which(BCestimate[,3] > 40 & BCestimate[,3] < 41 & BCestimate[,4]<0.05)) 
id40
par(mfrow=c(3,4))
for (i in id40){
  m3(ddd3[which(ddd3$z == i),])
}

###################plot c ~= 30 #######################
id30 = (which(BCestimate[,3] > 30 & BCestimate[,3] < 31 & BCestimate[,4]<0.05)) 
id30
par(mfrow=c(3,4))
for (i in id30){
  m3(ddd3[which(ddd3$z == i),])
}


##note: for now, run machine learning without the c>60 on on visit to see how it compares with the previous results in terms of R2

################### try just run machine learning without the c>60 on visit ###########################
######## remove the birthdate, dvdate etc which are used to calculate missing age ########
colnames(ddd3)[which(names(ddd3) == "x")] <- "age"
ddd3 <- subset(ddd2, select = -age)
ddd3 <- subset(ddd3, select = -y)  
ddd3 <- subset(ddd3, select = -w)
ddd3 <- subset(ddd3, select = -e)
ddd3 <- subset(ddd3, select = -h)

######## import estimated 'c' data ########
rawdata2 <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/OutputData/sigmoidalParameterEst.csv")

######## add c and pvalue for c to the clinical information data ########
ddd3$c<-NA
ddd3$cpval<-NA
for(i in (1:nrow(rawdata2))){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
}

ddd3[ddd3$z == 215,]
rawdata2[215,]
# Note: find out z=215 is empty

# remove the no.215
ddd3<-ddd3[-which((ddd$z == 215)),]

######add c and pvalue for c for each patient
for(i in (1:214)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
for(i in (216:242)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
ddd4<-ddd3[-which(is.na(ddd3$c)),] # remove the NA in exist c
ddd5 <- ddd4 %>% filter(cpval <= 0.05) # keep the c which p_value < 0.05
# keep the c which c < 60
ddd5 <- ddd5  %>% filter(c < 60) #

# delete ids and variables with a lot of NA
ddd5<-data.table(ddd5, key="z")
sort(colSums(is.na(ddd5))) 
sort(rowSums(is.na(ddd5)))
ddd5 <- ddd5[,-c("sldl","uprotc_ca","uprote_ca","esode_cc")] 
ddd5 <- ddd5[-which(rowSums(is.na(ddd5))>10),]  

# imputation for variables with NA using MICE, variables without NA just keep
ddd5imp <- ddd5 %>% select_if(colSums(is.na(ddd5))>0)
ddd5keep <- ddd5 %>% select_if(colSums(is.na(ddd5))==0)
imp <- mice(ddd5imp,maxit=5)
ddd55 <- complete(imp,5)
ddd555 <- cbind(ddd55,ddd5keep) 
sort(colSums(is.na(ddd555)))  # checked no missing values, this is the final data for use

#change the name of column for age
colnames(ddd555)[which(names(ddd555) == "x")] <- "age"

######################################## machine learning prediction ###################################
######################## prediction 1):using everything else + visit 1 biomarker #################
# select visit 1 biomarkers
ddd555
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd66 <- ddd6[,head(.SD,1),by=pkdid]  
# ddd66 is the visit 1 biomarkers and everything else

# checked no missing values, this is the final data for use
sort(rowSums(is.na(ddd66)))
sort(colSums(is.na(ddd66))) 
# convert to dataframe to run machine learning
final1 <- as.data.frame(ddd66) 

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre3[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy3 <- final1[,"c"] # observed true outcome value
mse3 <- sum((obsy3-testy_pre3)^2) / nrow(final1)  ## MSE calculation
mse3 #26.0397
rsq3 <- (cor(obsy3,testy_pre3))^2 ## R squared calculation
rsq3 # 0.6196209
MAE(testy_pre3, obsy3) ## MAE calculation
#3.737644
MAPE(testy_pre3, obsy3) #  mean absolute percentage error calculation
# 0.09410657

# Method = Random Forest to predict c
testy_pre3.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre3.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy3.2 <- final1[,"c"] # observed true outcome value
mse3.2 <- sum((obsy3.2-testy_pre3.2)^2) / nrow(final1) ## MSE calculation
mse3.2 # 36.8668
rsq3.2 <- (cor(obsy3.2,testy_pre3.2))^2 ## R squared calculation
rsq3.2 # 0.501132
MAE(testy_pre3.2, obsy3.2) ## MAE calculation
#4.7293
MAPE(testy_pre3.2, obsy3.2) #  mean absolute percentage error calculation
#0.1216615

# Method = SVM to predict c
final1$genetype <- nnet::class.ind(final1$genetype)  # Genetype need to be used as ind
final1$genetype <- final1$genetype[,2:4]
final1$genetype
testy_pre3.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre3.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy3.3 <- final1[,"c"] # observed true outcome value
mse3.3 <- sum((obsy3.3-testy_pre3.3)^2) / nrow(final1)  ## MSE calculation
mse3.3 #40.0301
rsq3.3 <- (cor(obsy3.3,testy_pre3.3))^2 ## R squared calculation
rsq3.3 #0.4585742
MAE(testy_pre3.3, obsy3.3) ## MAE calculation
#4.939079
MAPE(testy_pre3.3, obsy3.3) #  mean absolute percentage error calculation
#0.1277152



##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy,testy_pre,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy3,testy_pre3,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy1.2,testy_pre1.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy3.2,testy_pre3.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy1.3,testy_pre1.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy3.3,testy_pre3.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


#####plot the c values (predicted from machine learning model vs predicted from sigmoidal curve fitting). 
# for lasso
par(mfrow=c(2,2))
op<- par(cex = 0.5)
plot(obsy3,testy_pre3,lty="solid",col="black",cex=0.5,main="Fitted c  vs Predicted c",
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=1.2,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy3.2,testy_pre3.2,lty="solid",col="red",cex=0.5,main="Fitted c vs Predicted c",
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
op <- par(cex = 0.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=1.2,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy3.3,testy_pre3.3,lty="solid",col="green",cex=0.5,main="Fitted c vs Predicted c",
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=1.2,bty="n",
       seg.len=0.1,x.intersp=0.3)
abline(1,1)  # adding a diagonal line to show dispersion of dots


######################## prediction 2):using everything else + visit 2 biomarker #################
#derive visit 2 biomarkers
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd77 <- ddd6[,.SD[2L],by=pkdid]  
# ddd77 is the visit 2 biomarkers dataset

# hecked no missing values, this is the final data for use
ddd88<-as.data.frame(ddd77)
sort(rowSums(is.na(ddd88)))
sort(colSums(is.na(ddd88))) 
# note: have two missing rows for biomarkers

# delete the missing values
ddd88 <- na.omit(ddd88) 

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre4<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(ddd88$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- ddd88[-i, ] # train: the rest obs
  test <- ddd88[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre4[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy4 <- ddd88[,"c"] # observed true outcome value
mse4 <- sum((obsy4-testy_pre4)^2) / nrow(ddd88)  ## MSE calculation
mse4 # 25.62294
rsq4 <- (cor(obsy4,testy_pre4))^2 ## R squared calculation
rsq4 #0.6279293
MAE(testy_pre4, obsy4) ## MAE calculation
#3.578099
MAPE(testy_pre4, obsy4) #  mean absolute percentage error calculation
#0.09040555

# Method = Random Forest to predict c
testy_pre4.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(ddd88$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- ddd88[-i, ] # train: the rest obs
  test <- ddd88[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre4.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy4.2 <- ddd88[,"c"] # observed true outcome value
mse4.2 <- sum((obsy4.2-testy_pre4.2)^2) / nrow(ddd88) ## MSE calculation
mse4.2 # 37.20502
rsq4.2 <- (cor(obsy4.2,testy_pre4.2))^2 ## R squared calculation
rsq4.2 # 0.5108225
MAE(testy_pre4.2, obsy4.2) ## MAE calculation
#4.76788
MAPE(testy_pre4.2, obsy4.2) #  mean absolute percentage error calculation
#0.1230746


# Method = SVM to predict c
ddd88$genetype <- nnet::class.ind(ddd88$genetype)  # Genetype need to be used as ind
ddd88$genetype <- ddd88$genetype[,2:4] #delete all 0 column for SVM to running
testy_pre4.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(ddd88$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- ddd88[-i, ] # train: the rest obs
  test <- ddd88[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre4.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy4.3 <- ddd88[,"c"] # observed true outcome value
mse4.3 <- sum((obsy4.3-testy_pre4.3)^2) / nrow(ddd88)  ## MSE calculation
mse4.3 # 37.57847
rsq4.3 <- (cor(obsy4.3,testy_pre4.3))^2 ## R squared calculation
rsq4.3 # 0.4970028
MAE(testy_pre4.3, obsy4.3) ## MAE calculation
#4.854456
MAPE(testy_pre4.3, obsy4.3) #  mean absolute percentage error calculation
#0.1257003


##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy,testy_pre,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy2,testy_pre2,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy1.2,testy_pre1.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy2.2,testy_pre2.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy1.3,testy_pre1.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy2.3,testy_pre2.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# Note: From the plot, there are several obviously outlier around c=20, we can draw all the subjects for the sigmoidal fit and pick some subjects not successfully fit and leave out these subjects.

################ draw all the subjects to check their fit ###################
somePDFPath = "/Users/mengjiawei/Desktop/PKD/plot.pdf"
pdf(file=somePDFPath)  
par(mfrow = c(3,3))
for (i in idx)   
{   
  m3(ddd3[which(ddd3$z == i),]) 
} 
dev.off() 
#############################################

# Results: Leaving out some subjects (2,28,34,35,40,42,49,50,51,59,62,63,64,66,71,73,78,
# 79, 81,89,91,92,98,99,101,111,115,119,121,135,130,145,146,148,151,159,161,168,169,172,179,186,188,193,196,216,235,240) 
# that obviously can not successfully fit 

# note: will use these successful fit data for below analysis
######## remove the birthdate, dvdate etc which are used to calculate missing age ########
colnames(ddd3)[which(names(ddd3) == "x")] <- "age"
ddd3 <- subset(ddd2, select = -age)
ddd3 <- subset(ddd3, select = -y)  
ddd3 <- subset(ddd3, select = -w)
ddd3 <- subset(ddd3, select = -e)
ddd3 <- subset(ddd3, select = -h)

######## import estimated 'c' data ########
rawdata2 <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/OutputData/sigmoidalParameterEst.csv")

######## add c and pvalue for c to the clinical information data ########
ddd3$c<-NA
ddd3$cpval<-NA
for(i in (1:nrow(rawdata2))){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
}

ddd3[ddd3$z == 215,]
rawdata2[215,]
# Note: find out z=215 is empty

# remove the no.215
ddd3<-ddd3[-which((ddd$z == 215)),]

######add c and pvalue for c for each patient
for(i in (1:214)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
for(i in (216:242)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
ddd4<-ddd3[-which(is.na(ddd3$c)),] # remove the NA in exist c
ddd5 <- ddd4 %>% filter(cpval <= 0.05) # keep the c which p_value < 0.05
ddd5 <- ddd5  %>% filter(c < 60) # keep the c which c < 60

# delete ids and variables with a lot of NA
ddd5<-data.table(ddd5, key="z")
sort(colSums(is.na(ddd5))) 
sort(rowSums(is.na(ddd5)))
ddd5 <- ddd5[,-c("sldl","uprotc_ca","uprote_ca","esode_cc")] 
ddd5 <- ddd5[-which(rowSums(is.na(ddd5))>10),]  

################leave out patients that cannot fit well###############
ddd5<-ddd5[-which((ddd5$z %in% c(2,28,34,35,40,42,49,50,51,59,62,63,64,66,71,73,78,79,
                                 81,89,91,92,98,99,101,111,115,119,121,135,130,145,146,
                                 148,151,159,161,168,169,172,179,186,188,193,196,216,235,240))),]

#######################################################################
# imputation for variables with NA using MICE, variables without NA just keep
ddd5imp <- ddd5 %>% select_if(colSums(is.na(ddd5))>0)
ddd5keep <- ddd5 %>% select_if(colSums(is.na(ddd5))==0)
imp <- mice(ddd5imp,maxit=5)
ddd55 <- mice::complete(imp,5)
ddd555 <- cbind(ddd55,ddd5keep) 
sort(colSums(is.na(ddd555)))  # checked no missing values, this is the final data for use

#change the name of column for age
colnames(ddd555)[which(names(ddd555) == "x")] <- "age"


######################################## machine learning prediction ###################################
######################## prediction 1):using everything else + visit 1 biomarker #################
# select visit 1 biomarkers
ddd555
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd66 <- ddd6[,head(.SD,1),by=pkdid]  
# ddd66 is the visit 1 biomarkers and everything else
## checked no missing values, this is the final data for use
sort(rowSums(is.na(ddd66)))
sort(colSums(is.na(ddd66))) 
## convert to dataframe to run machine learning
final1 <- as.data.frame(ddd66) 

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre5<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre5[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy5 <- final1[,"c"] # observed true outcome value
rsq5 <- (cor(obsy5,testy_pre5))^2 ## R squared calculation
rsq5 # 0.6765279
MAE(testy_pre5, obsy5) ## MAE calculation
#3.364389
MAPE(testy_pre5, obsy5) #  mean absolute percentage error calculation
# 0.08020167

# Method = Random Forest to predict c
testy_pre5.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre5.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy5.2 <- final1[,"c"] # observed true outcome value
rsq5.2 <- (cor(obsy5.2,testy_pre5.2))^2 ## R squared calculation
rsq5.2 # 0.5657328
MAE(testy_pre5.2, obsy5.2) ## MAE calculation
#4.293558
MAPE(testy_pre5.2, obsy5.2) #  mean absolute percentage error calculation
#0.1065778

# Method = SVM to predict c
final1$genetype <- nnet::class.ind(final1$genetype)  # Genetype need to be used as ind
final1$genetype <- final1$genetype[,2:4]
final1$genetype
testy_pre5.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre5.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy5.3 <- final1[,"c"] # observed true outcome value
rsq5.3 <- (cor(obsy5.3,testy_pre5.3))^2 ## R squared calculation
rsq5.3 #0.53484
MAE(testy_pre5.3, obsy5.3) ## MAE calculation
# 4.42397
MAPE(testy_pre5.3, obsy5.3) #  mean absolute percentage error calculation
#0.109471


##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy3,testy_pre3,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy5,testy_pre5,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy3.2,testy_pre3.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy5.2,testy_pre5.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy3.3,testy_pre3.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy5.3,testy_pre5.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots




################### find the outlier in the plot #################
abs(obsy-testy_pre)
which.max(abs(obsy-testy_pre)) #21
which.max(abs(obsy1.2-testy_pre1.2)) #54

(abs(obsy-testy_pre) > 10) #7, 21,54,88

final1[7,] #pkdid = 108730
final1[21,] #pkdid = 120777
final1[54,] #pkdid = 213454
final1[88,] #pkdid = 300696

ddd555
ddd555[ddd555$pkdid==108730,] #z=9
ddd555[ddd555$pkdid==120777,] #z=23
ddd555[ddd555$pkdid==213454,] #z=93
ddd555[ddd555$pkdid==300696,] #z=156

# Note: For the absolute value of difference for the observed and predicted c larger than 10 is #7(pkdid = 108730, z=9), 21(pkdid = 120777, z==23), 54(pkdid = 213454, z=93), 88(pkdid = 300696, z=156) patient. 

####################### prediction 2):using everything else + visit 2 biomarker #################
#derive visit 2 biomarkers
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd77 <- ddd6[,.SD[2L],by=pkdid]  
# ddd77 is the visit 2 biomarkers dataset

ddd88<-as.data.frame(ddd77)
sort(rowSums(is.na(ddd88)))
sort(colSums(is.na(ddd88))) ## checked no missing values, this is the final data for use
# note: have two missing rows for biomarkers

ddd88 <- na.omit(ddd88) ## delete the missing values

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre6<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(ddd88$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- ddd88[-i, ] # train: the rest obs
  test <- ddd88[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre6[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy6 <- ddd88[,"c"] # observed true outcome value
rsq6 <- (cor(obsy6,testy_pre6))^2 ## R squared calculation
rsq6 #0.7318242
MAE(testy_pre6, obsy6) ## MAE calculation
#2.939634
MAPE(testy_pre6, obsy6) #  mean absolute percentage error calculation
#0.06942608

# Method = Random Forest to predict c
testy_pre6.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(ddd88$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- ddd88[-i, ] # train: the rest obs
  test <- ddd88[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre6.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy6.2 <- ddd88[,"c"] # observed true outcome value
rsq6.2 <- (cor(obsy6.2,testy_pre6.2))^2 ## R squared calculation
rsq6.2 #0.5461713
MAE(testy_pre6.2, obsy6.2) ## MAE calculation
#4.45372
MAPE(testy_pre6.2, obsy6.2) #  mean absolute percentage error calculation
#0.1106011


# Method = SVM to predict c
ddd88$genetype <- nnet::class.ind(ddd88$genetype)  # Genetype need to be used as ind
ddd88$genetype <- ddd88$genetype[,2:4]
ddd88$genetype
testy_pre6.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(ddd88$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- ddd88[-i, ] # train: the rest obs
  test <- ddd88[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre6.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy6.3 <- ddd88[,"c"] # observed true outcome value
rsq6.3 <- (cor(obsy6.3,testy_pre6.3))^2 ## R squared calculation
rsq6.3 #0.402714
MAE(testy_pre6.3, obsy6.3) ## MAE calculation
# 4.939207
MAPE(testy_pre6.3, obsy6.3) #  mean absolute percentage error calculation
#0.1217987


##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy,testy_pre,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy2,testy_pre2,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy1.2,testy_pre1.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy2.2,testy_pre2.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy1.3,testy_pre1.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy2.3,testy_pre2.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:120),ylim=range(0:120),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

########## Using clinical information (add health) predict inflection point parameter (c) #######
########## original data cleaning 
######## import data(after adding health data) ########
rawdata <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/rawdata4.csv")

######## delete the data which will not be used ########
raw1 <- rawdata[,-which(names(rawdata) %in% c("pkdid5","vis","visc",
                                              "race","ff","gene_drharris","il18","ngal","trunc_grp","brwgt"))]
# delete "ff", "il18","ngal" because they are all NA; gene_drharris is same as genetype; 

raw1$pkdid <- as.factor(raw1$pkdid)
d <- data.table(raw1, key="pkdid") # transform to data table
d<-d[-which(is.na(d$pkdid)),]

######## transform the date data to date type ########
d$dvdate <- as.Date(d$dvdate,format = "%m/%d/%Y")
d$basedate <- as.Date(d$basedate,format = "%m/%d/%Y")
d$birthdate <- as.Date(d$birthdate,format = "%m/%d/%Y")

######## calculate the visit time sequence for each patient ########
dd2<-getanID(data = d, id.vars = "pkdid")
d$yrs <- dd2$.id

######## calculating age ########
# x is the age
x <- rep(NA,length(d$age))

# y is the ckd_epi for each visit for each patient
y <- d$ckd_epi

# z is the patient group (number 1 is for first patient..)
z<- d$pkdid

# w is the visit date
w<- d$dvdate

# e is the birthdate for each patient
e<- d$birthdate

# h is the visit time sequence
h<- d$yrs

# conbine these to a data frame dd
dd<-cbind(x,y,z,w,e,h)
dd<-as.data.frame(dd)

# add this 'dd' data frame to original data frame 'd'
newd<- cbind(dd,d)
ddd<-as.data.frame(newd)
ddd$e<-as.Date(e)
ddd$w<-as.Date(w)

# fill missing birthdate using same birthdate from the same person
for (s in 1:max(ddd$z)){
  ddd$e[ddd$z == s][which(is.na(ddd$e[ddd$z == s]))] <- (ddd$e[ddd$z == s])[1]
}

# get each age data for patients using birthdate and dvdate
for (f in 1:max(ddd$z)){
  diff <- round(as.numeric(difftime(ddd$w[ddd$z == f], ddd$e[ddd$z==f]))/365.25,1)
  ddd$x[ddd$z == f] <- c(diff)
}

# remove the missing age row in data (because that one doesn't have dvdata)
ddd1<-ddd[-which(is.na(ddd$x)),]

# remove junk
for (f in 1:max(ddd1$z)){
  ddd1$x[ddd1$z == f][duplicated(ddd1$x[ddd1$z == f],fromLast=TRUE)] <- NA
}
ddd2<-ddd1[-which(is.na(ddd1$x)),]

# remove uncompleted age column from original clinical information data frame, keep the completed age
ddd3 <- subset(ddd2, select = -age)
ddd3 #ddd3 is the clinical information for each patient which has no missing age data

##################################### nlsLM function in minpack #####################################
######## write a function for plotting ########
m3 = function(data){
  x = data$x
  y = data$y
  
  # find the best starting point for c
  c.start = x[which(abs(y-60)==min(abs(y-60),na.rm=T))]
  
  # fit model
  fitmodel <- nlsLM(y ~(5+((110-5)/(1 + (x/c)^b))),start=c(b=0,c=c.start))
  b = round(summary(fitmodel)$parameters[1,1],2)
  c = round(summary(fitmodel)$parameters[2,1],1)
  plot(x, y, ylim=c(0,130), xlab="Age", ylab="eGFR", main=paste("SubjectID = ", eval(i)))
  legend("bottomleft", c( paste("slope = ", eval(b)),paste("Age_50% = ", eval(c)) ), bty="n")
  lines(x[!is.na(y)], fitted(fitmodel),col = 2, lwd = 2)  
}

## Notes:  a  small portion still do not fit well  due to flat eGFR over time.

######## fit on everyone and get c,b and p values ######## 
n = max(ddd3$z)
bcEst = matrix(NA,n, 5)
names(bcEst) =c("b","pval_b","c","pval_c", "cor_coef")

for(i in (1:n)){
  print(i) # monitor the progress
  data = ddd3[which(ddd3$z == i),]
  x = data$x
  y = data$y
  c.start = x[which(abs(y-60)==min(abs(y-60),na.rm=T))]
  
  # fit model (Errors are generated for some subjects.  Use tryCatch to skip them)
  fitmodel <- tryCatch(nlsLM(y ~(5+((110-5)/(1 + (x/c)^b))),start=c(b=0,c=c.start)),error = function(e) NULL)
  # get the results for the subjects that works for fitting the model
  if (!is.null(fitmodel)){
    # get the output
    bcEst[i,1] = summary(fitmodel)$parameters[1,1]
    bcEst[i,2] = summary(fitmodel)$parameters[1,4]
    bcEst[i,3] = summary(fitmodel)$parameters[2,1]
    bcEst[i,4] = summary(fitmodel)$parameters[2,4]
    bcEst[i,5] = cor(y[!is.na(y)], predict(fitmodel)) 
  }
} 

bcEst = matrix(NA,n, 6)
names(bcEst) =c("b","pval_b","c","pval_c", "cor_coef", "ESRD_age")

for(i in (1:n)){
  print(i) # monitor the progress
  data = ddd3[which(ddd3$z == i),]
  x = data$x
  y = data$y
  c.start = x[which(abs(y-60)==min(abs(y-60),na.rm=T))]
  # fit model (Errors are generated for some subjects.  Use tryCatch to skip them)
  fitmodel <- tryCatch(nlsLM(y ~(5+((110-5)/(1 + (x/c)^b))),start=c(b=0,c=c.start)),error = function(e) NULL)
  # get the results for the subjects that works for fitting the model
  if (!is.null(fitmodel)){
    # get the output
    bcEst[i,1] = summary(fitmodel)$parameters[1,1]
    bcEst[i,2] = summary(fitmodel)$parameters[1,4]
    bcEst[i,3] = summary(fitmodel)$parameters[2,1]
    bcEst[i,4] = summary(fitmodel)$parameters[2,4]
    bcEst[i,5] = cor(y[!is.na(y)], predict(fitmodel)) 
    b= summary(fitmodel)$parameters[1,1]
    c= summary(fitmodel)$parameters[2,1]
    bcEst[i,6] = c*((110-15)/(15-5))^(1/b)
  }
} 

# Note: For some reason, the output matrix has some extra stuff at the end.  It is likely related to the matrix defination.  Try to change it to data frame in the future.

# round the coef and pval outputs
bcEstimate = cbind(round(bcEst[1:242,1],2 ), round(bcEst[1:242,2],3 ),round(bcEst[1:242,3],1 ), round(bcEst[1:242,4],3 ),round(bcEst[1:242,5],3 ),round(bcEst[1:242,6], 1 ) ) 
colnames(bcEstimate) =c("b","pval_b","c","pval_c", "cor_coef", "ESRD_age")

# remove the extra stuff
BCestimate = bcEstimate[1:242,]
write.csv(BCestimate, file="/Users/mengjiawei/Documents/CRISP-PKD/DATA/OutputData/sigmoidalParameterEst2.csv") 


### Look at the R2 of the subjects with successful fit
idx = (which(BCestimate[,4]<0.05)) #176 (73%)
idx

### draw all the subjects that successfully fit
somePDFPath = "/Users/mengjiawei/Desktop/PKD/plot.pdf"
pdf(file=somePDFPath)  
par(mfrow = c(3,3))
for (i in idx)   
{   
  m3(ddd3[which(ddd3$z == i),]) 
} 
dev.off() 
#############################################
hist(BCestimate[idx,5])
hist(BCestimate[idx,5]^2) # 116 with R2 >0.7

################### combine estimated c from fit and clinical info
######## remove the birthdate, dvdate etc which are used to calculate missing age ########
colnames(ddd3)[which(names(ddd3) == "x")] <- "age"
ddd3 <- subset(ddd2, select = -age)
ddd3 <- subset(ddd3, select = -y)  
ddd3 <- subset(ddd3, select = -w)
ddd3 <- subset(ddd3, select = -e)
ddd3 <- subset(ddd3, select = -h)

######## import estimated 'c' data ########
rawdata2 <- read.csv("/Users/mengjiawei/Documents/CRISP-PKD/DATA/OutputData/sigmoidalParameterEst2.csv")

######## add c and pvalue for c to the clinical information data ########
ddd3$c<-NA
ddd3$cpval<-NA
for(i in (1:nrow(rawdata2))){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
}

ddd3[ddd3$z == 215,]
rawdata2[215,]
# Note: find out z=215 is empty

# remove the no.215
ddd3<-ddd3[-which((ddd$z == 215)),]

######add c and pvalue for c for each patient
for(i in (1:214)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
for(i in (216:242)){
  print(i) # monitor the progress
  ddd3[ddd3$z == i,]$c <- rawdata2[i,]$c
  ddd3[ddd3$z == i,]$cpval <- rawdata2[i,]$pval_c
}
ddd4<-ddd3[-which(is.na(ddd3$c)),] # remove the NA in exist c
ddd5 <- ddd4 %>% filter(cpval <= 0.05) # keep the c which p_value < 0.05
ddd5 <- ddd5  %>% filter(c < 60) # keep the c which c < 60

# delete ids and variables with a lot of NA
ddd5<-data.table(ddd5, key="z")
sort(colSums(is.na(ddd5))) 
sort(rowSums(is.na(ddd5)))
ddd5 <- ddd5[,-c("sldl","uprotc_ca","uprote_ca","esode_cc")] 
ddd5 <- ddd5[-which(rowSums(is.na(ddd5))>10),]  

################leave out patients that cannot fit well###############
ddd5<-ddd5[-which((ddd5$z %in% c(2,28,34,35,40,42,49,50,51,59,62,63,64,66,71,73,78,79,
                                 81,89,91,92,98,99,101,111,115,119,121,135,130,145,146,
                                 148,151,159,161,168,169,172,179,186,188,193,196,216,235,240))),]

#######################################################################
# imputation for variables with NA using MICE, variables without NA just keep
ddd5imp <- ddd5 %>% select_if(colSums(is.na(ddd5))>0)
ddd5imp
ddd5keep <- ddd5 %>% select_if(colSums(is.na(ddd5))==0)
imp <- mice(ddd5imp,maxit=5)
ddd55 <- mice::complete(imp,5)
ddd555 <- cbind(ddd55,ddd5keep) 
sort(colSums(is.na(ddd555)))  # checked no missing values, this is the final data for use

#change the name of column for age
colnames(ddd555)[which(names(ddd555) == "x")] <- "age"

######################################## machine learning prediction ###################################
######################## prediction 1):using everything else + visit 1 biomarker #################
# select visit 1 biomarkers
ddd555
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval)) #visc will all be the baeline:0
ddd6 <- data.table(ddd6, key="pkdid")
ddd66 <- ddd6[,head(.SD,1),by=pkdid]  
# ddd66 is the visit 1 biomarkers and everything else
sort(rowSums(is.na(ddd66)))
sort(colSums(is.na(ddd66))) ## checked no missing values, this is the final data for use
final1 <- as.data.frame(ddd66) ## convert to dataframe to run machine learning
dim(final1)

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre7<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre7[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy7 <- final1[,"c"] # observed true outcome value
rsq7 <- (cor(obsy7,testy_pre7))^2 ## R squared calculation
rsq7 #0.6723594
MAE(testy_pre7, obsy7) ## MAE calculation
#3.354925
MAPE(testy_pre7, obsy7) #  mean absolute percentage error calculation
#0.07958759


# Method = Random Forest to predict c
testy_pre7.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))]  # trainx: predictors for model building, "c .mean", "z .mean", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre7.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy7.2 <- final1[,"c"] # observed true outcome value
rsq7.2 <- (cor(obsy7.2,testy_pre7.2))^2 ## R squared calculation
rsq7.2 # 0.5690382
MAE(testy_pre7.2, obsy7.2) ## MAE calculation
# 4.245902
MAPE(testy_pre7.2, obsy7.2) #  mean absolute percentage error calculation
#0.1049742

# Method = SVM to predict c
final1$genetype <- nnet::class.ind(final1$genetype)  # Genetype need to be used as ind
final1$genetype <- final1$genetype[,2:4]
final1$genetype
final1$visc
testy_pre7.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final1$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final1[-i, ] # train: the rest obs
  test <- final1[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c","z","pkdid"))] # trainx: predictors for model building, "c", "z", "pkdid" should not be used
  trainy <- train[,c("c")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c","z","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre7.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy7.3 <- final1[,"c"] # observed true outcome value
rsq7.3 <- (cor(obsy7.3,testy_pre7.3))^2 ## R squared calculation
rsq7.3 #  0.4958824
MAE(testy_pre7.3, obsy7.3) ## MAE calculation
#4.531943
MAPE(testy_pre7.3, obsy7.3) #  mean absolute percentage error calculation
#0.1123968



##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy5,testy_pre5,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy7,testy_pre7,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy5.2,testy_pre5.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy7.2,testy_pre7.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy5.3,testy_pre5.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy7.3,testy_pre7.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


######################## prediction using all three year biomarkers #################
#derive 3 year data
ddd555
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,cpval,visc))
ddd6 <- data.table(ddd6, key="pkdid")
ddd66 <- ddd6[,head(.SD,3),by=pkdid]  #get 3 year data for all patients
head(ddd66)

#each patient should have at least 3 visit records
check <- as.data.frame(table(ddd66$yrs))
del <- as.matrix(check[check$Freq < 120,])
ddd66 <- ddd66[!(ddd66$yrs %in% del[,1]),]
check2 <- as.data.frame(table(ddd66$pkdid))
del2 <- as.matrix(check2[check2$Freq < 3,])
ddd66 <- ddd66[!(ddd66$pkdid %in% del2[,1]),]

## checked no missing values
sort(rowSums(is.na(ddd66)))
sort(colSums(is.na(ddd66)))

# change the long form to wide form
ddd66<- as.data.frame(ddd66)
ddd777<-reshape(ddd66,direction = "wide",idvar = "pkdid",timevar = "yrs")
ddd777_2 <- ddd777 %>% select_if(colSums(is.na(ddd777))<10)
final4 <- subset(ddd777_2, select = -c(`c.1`,`c.2`,`z.1`,
                                       `z.2`,`z.3`,`educ.1`,`educ.2`,`gender.1`,`gender.2`,`race4.1`,`race4.2`,`genetype.1`,`genetype.2`))
# check for missing value
sort(rowSums(is.na)(final4))
sort(colSums(is.na(final4))) 
# convert c to numeric
final4$`c.3`<-as.numeric(final4$`c.3`)
str(final4)
head(final4)
# convert to wide form 

ddd666<- dcast(melt(ddd66, id.vars = c("pkdid", "yrs")), pkdid ~ yrs + variable)

ddd666_2 <- ddd666 %>% select_if(colSums(is.na(ddd666))<10)
final3 <- subset(ddd666_2, select = -c(`1_c`,`2_c`,`1_z`,
                                       `2_z`,`3_z`,`1_educ`,`2_educ`,`1_gender`,`2_gender`,`1_race4`,`2_race4`,`1_genetype`,`2_genetype`))
# check for missing value
sort(rowSums(is.na)(final3))
sort(colSums(is.na(final3))) 
# convert c to numeric
final3$`3_c`<-as.numeric(final3$`3_c`)
str(final3)
head(final3)
dim(final3)

### machine learning predict###
# leave one out 
# Method = Lasso to predict c
testy_pre8<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final4$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final4[-i, ] # train: the rest obs
  test <- final4[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c.3","pkdid"))] # trainx: predictors for model building, "c.3", "pkdid" should not be used
  trainy <- train[,c("c.3")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c.3","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_lasso <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.glmnet") # Lasso regression model
  testy_pre8[i] <- sl_lasso$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy8 <- final4[,"c.3"] # observed true outcome value
rsq8 <- (cor(obsy8,testy_pre8))^2 ## R squared calculation
rsq8 #0.7113394
MAE(testy_pre8, obsy8) ## MAE calculation
#3.110319
MAPE(testy_pre8, obsy8) #  mean absolute percentage error calculation
#0.07428092


# Method = Random Forest to predict c
testy_pre8.2 <- vector() # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final4$pkdid)) # loop for all observations
{
  set.seed(123456)
  train <- final4[-i, ] # train: the rest obs
  test <- final4[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c.3","pkdid"))]  # trainx: predictors for model building, "c.3", "pkdid" should not be used
  trainy <- train[,c("c.3")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c.3","pkdid"))] # testx:  predictors for the left out obs prediction
  testy <- test[,c("c.3")]
  sl_rf <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.randomForest") # Random Forest model
  testy_pre8.2[i] <- sl_rf$SL.predict # the predicted outcome value of c, after the loop get a vector
}
obsy8.2 <- final4[,"c.3"] # observed true outcome value
rsq8.2 <- (cor(obsy8.2,testy_pre8.2))^2 ## R squared calculation
rsq8.2 # 0.6086483
MAE(testy_pre8.2, obsy8.2) ## MAE calculation
#3.974922
MAPE(testy_pre8.2, obsy8.2) #  mean absolute percentage error calculation
#0.09790185


# Method = SVM to predict c
final4$genetype.3 <- nnet::class.ind(final4$genetype.3)  # Genetype need to be used as ind
final4$genetype.3 <- final4$genetype.3[,2:4]
final4$genetype.3
testy_pre8.3<- vector()  # the vector of predicted 'c' value, initialled empty
for (i in 1:length(final4$pkdid))  # loop for all observations
{
  set.seed(123456) # set a seed for SuperLearner model
  train <- final4[-i, ] # train: the rest obs
  test <- final4[i, ] # test: the left out obs (obs_i)
  trainx <- train[,!(names(train) %in% c("c.3","pkdid"))] # trainx: predictors for model building, "c .mean", "z .mean", "pkdid" should not be used
  trainy <- train[,c("c.3")] # trainy: outcome for model building
  testx <- test[,!(names(test) %in% c("c.3","pkdid"))] # testx:  predictors for the left out obs prediction
  sl_svm <- SuperLearner(Y=trainy,X=trainx,newX=testx,SL.library="SL.svm")
  testy_pre8.3[i] <- sl_svm$SL.predict # the predicted outcome value of c, after the loop get a vector
}

obsy8.3 <- final4[,"c.3"]  # observed true outcome value
rsq8.3 <- (cor(obsy8.3,testy_pre8.3))^2 ## R squared calculation
rsq8.3 #  0.5189616
MAE(testy_pre8.3, obsy8.3) ## MAE calculation
# 4.48
MAPE(testy_pre8.3, obsy8.3) #  mean absolute percentage error calculation
#0.11


##plot together

par(mfrow=c(3,2))
op<- par(cex=0.5)
plot(obsy5,testy_pre5,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05, x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for lasso
plot(obsy8,testy_pre8,lty="solid",col="black",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("Lasso"),
       col=c("black"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for random forest
plot(obsy5.2,testy_pre5.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for random forest
plot(obsy8.2,testy_pre8.2,lty="solid",col="red",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("RandomForest"),
       col=c("red"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# for svm
plot(obsy5.3,testy_pre5.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots

# for svm
plot(obsy8.3,testy_pre8.3,lty="solid",col="green",cex=0.5,
     ylab="Predicted c",xlab="Fitted c",xlim=range(0:95),ylim=range(0:95),cex.axis=1.5, cex.lab=1.5)
legend("topleft",legend=c("SVM"),
       col=c("green"),lty=1,cex=2.0,bty="n",
       seg.len=0.05,x.intersp=0.01)
abline(1,1)  # adding a diagonal line to show dispersion of dots


# Note: Compared to only use visit 1 biomarker, the R2 for lasso and random forest are increased, the MAE and MAPE are decreased for using all three years biomarkers.


###### Plotted the non-zero median coefficients of the predictors in Lasso regression for prediction of c value ########3
## Std final1: use standardized value of continuous predictors, this does not affect the prediction result,
##              but will be further used for significant coefficient identification

library(dplyr)
library(data.table)
library(car)
library(caret)
library(mice)
library(randomForest)
library(ranger)
library(nnet)
library(e1071)
library(ggplot2)
library(reshape2)
library(MASS)
# Model: visit 1 

coef_lasso <- list()
library(SuperLearner)

final1s1 <- final1
ind <- sapply(final1s1,is.numeric)              # select numeric variables
final1s1[ind] <- lapply(final1s1[ind],scale)   # standardize numeric variables
final1s1$c <- final1$c         # c should not be scaled
# checked all variables are scaled except for non-numeric variables and c
final1s1


###  Select important variables from Lasso non-zero coefficients 
##      Use standardized predictors
##      Use loop to perform 10-fold cross validation, 9/10 for training and 1/10 for testing
##      SL internally generates 100 folds for cv. Select the median of 100 betas

# Model: baseline visit 
coef_lasso <- list()
for (i in 1:10)
{
  train_ind <- seq((i-1)*23+1,i*23,1)  # 10-fold external cv
  train <- final1s1[-train_ind,]
  x <- train[,!(names(train) %in% c("c","pkdid","z"))]
  y <- train[,c("c")]
  sl_lasso <- SuperLearner(Y=y,X=x,SL.library="SL.glmnet")
  coef_lasso_all <- as.matrix(sl_lasso$fitLibrary$SL.glmnet_All$object$glmnet.fit$beta)  # Extract internal cv coef matrix
  coef_lasso[[i]] <- as.matrix(apply(coef_lasso_all,1,median)) # select only the median coef for each variable among internal cv
}

coef1 <- do.call(cbind, coef_lasso) # the median Lasso beta coefficient matrix, n=10 as of external cv
coef1.1 <- as.data.frame(t(coef1[which(abs(rowSums(coef1))>1),])) # select significant variables by biggest sum beta values
coef1.1 <- melt(coef1.1)  # wide -> long format
coef1.1$rowid <- 1:10     # assign for each variable coef1 - coef10

# sort variables by max coef value
coef1.1p <- coef1.1 %>% group_by(variable) %>% 
  mutate(mx = max(value)) %>% 
  arrange(desc(mx), desc(value))  
coef1.1.1 <- dplyr::select(coef1.1p,-mx)
coef1.1.1$variable <- factor(coef1.1.1$variable,levels=rev(unique(coef1.1.1$variable)))
coef1.1.1

# assign group of Demo, Biomarker or Imaging
group <- ifelse(coef1.1.1$variable=='age','demographic',
                ifelse(coef1.1.1$variable=='gender','demographic',
                       ifelse(coef1.1.1$variable=='race4','demographic',
                              ifelse(coef1.1.1$variable=='educ','demographic',
                                     ifelse(coef1.1.1$variable=='visc','demographic',
                                            ifelse(coef1.1.1$variable=='tkv','imaging',
                                                   ifelse(coef1.1.1$variable=='agehtn','health',
                                                          ifelse(coef1.1.1$variable=='systole','health',
                                                                 ifelse(coef1.1.1$variable=='map','health',
                                                                        ifelse(coef1.1.1$variable=='diastol','health',
                                                                               'biomarker'))))))))))
coef1.1.1$group <- group

# Plot
Lasso_Coef1 <- ggplot(data=coef1.1.1, aes(x=variable,y=value,group=factor(rowid))) + 
  geom_jitter(aes(color = group)) + 
  geom_hline(yintercept=0) +
  ggtitle("Effect size of baseline visit health, demographic, biomarker and genetic predictors for c") +
  ylab("Beta Coefficients") +
  theme(plot.title=element_text(size=15,face="bold"),text=element_text(size=13),axis.text.x=element_text(angle=90))
ggsave("/Users/mengjiawei/Desktop/PKD/Lasso_Coef1o.png",scale=1,width=24,height=16,units="cm")


##### visit1
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd77 <- ddd6[,.SD[2L],by=pkdid]  
# ddd77 is the visit 1 biomarkers dataset

ddd88<-as.data.frame(ddd77)
sort(rowSums(is.na(ddd88)))
sort(colSums(is.na(ddd88))) ## checked no missing values, this is the final data for use
# note: have two missing rows for biomarkers

ddd88 <- na.omit(ddd88) ## delete the missing values

coef_lasso2 <- list()

for (i in 1:10)
{
  train_ind <- seq((i-1)*23+1,i*23,1)
  train <- ddd88[-train_ind,]
  x <- train[,!(names(train) %in% c("c","pkdid","z"))]
  y <- train[,c("c")]
  sl_lasso <- SuperLearner(Y=y,X=x,SL.library="SL.glmnet")
  coef_lasso_all <- as.matrix(sl_lasso$fitLibrary$SL.glmnet_All$object$glmnet.fit$beta)
  coef_lasso2[[i]] <- as.matrix(apply(coef_lasso_all,1,median))
}

coef2 <- do.call(cbind, coef_lasso2)
coef2.1 <- as.data.frame(t(coef2[which(abs(rowSums(coef2))>0.4),]))
coef2.1 <- melt(coef2.1) 
coef2.1$rowid <- 1:10

# sort variables by max coef value
coef2.1 <- coef2.1 %>% group_by(variable) %>% 
  mutate(mx = max(value)) %>% 
  arrange(desc(mx), desc(value)) 
coef2.1.1 <- dplyr::select(coef2.1,-mx)
coef2.1.1$variable <- factor(coef2.1.1$variable,levels=rev(unique(coef2.1.1$variable)))


# assign group of Demo, Biomarker or Imaging
group2 <- ifelse(coef2.1.1$variable=='age','demographic',
                 ifelse(coef2.1.1$variable=='gender','demographic',
                        ifelse(coef2.1.1$variable=='race4','demographic',
                               ifelse(coef2.1.1$variable=='educ','demographic',
                                      ifelse(coef2.1.1$variable=='visc','demographic',
                                             ifelse(coef2.1.1$variable=='tkv','imaging',
                                                    ifelse(coef2.1.1$variable=='agehtn','health',
                                                           ifelse(coef2.1.1$variable=='systole','health',
                                                                  ifelse(coef2.1.1$variable=='map','health',
                                                                         ifelse(coef2.1.1$variable=='diastol','health',
                                                                                'biomarker'))))))))))
coef2.1.1$group <- group2

# Plot
Lasso_Coef2 <- ggplot(coef2.1.1, aes(variable,value,group=factor(rowid))) + 
  geom_jitter(aes(color = group)) + 
  geom_hline(yintercept=0) +
  ggtitle("Effect size of visit 1 health, demographic, biomarker and genetic predictors for c") +
  ylab("Beta Coefficients") +
  theme(plot.title=element_text(size=15,face="bold"),text=element_text(size=13),axis.text.x=element_text(angle=90))
ggsave("/Users/mengjiawei/Desktop/PKD/Lasso_Coef2o.png",scale=1,width=24,height=16,units="cm")



#visit2
ddd3 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd3 <- data.table(ddd3, key="pkdid")
ddd33 <- ddd3[,.SD[3L],by=pkdid]  
# ddd33 is the visit 3 biomarkers dataset

ddd333<-as.data.frame(ddd33)
sort(rowSums(is.na(ddd333)))
sort(colSums(is.na(ddd333))) ## checked no missing values, this is the final data for use
# note: have two missing rows for biomarkers

ddd333 <- na.omit(ddd333) ## delete the missing values


coef_lasso3 <- list()

for (i in 1:10)
{
  train_ind <- seq((i-1)*23+1,i*23,1)
  train <- ddd333[-train_ind,]
  x <- train[,!(names(train) %in% c("c","pkdid","z"))]
  y <- train[,c("c")]
  sl_lasso <- SuperLearner(Y=y,X=x,SL.library="SL.glmnet")
  coef_lasso_all <- as.matrix(sl_lasso$fitLibrary$SL.glmnet_All$object$glmnet.fit$beta)
  coef_lasso3[[i]] <- as.matrix(apply(coef_lasso_all,1,median))
}

coef3 <- do.call(cbind, coef_lasso3)
coef3.1 <- as.data.frame(t(coef3[which(abs(rowSums(coef3))>0.4),]))
coef3.1 <- melt(coef3.1) 
coef3.1$rowid <- 1:10

# sort variables by max coef value
coef3.1 <- coef3.1 %>% group_by(variable) %>% 
  mutate(mx = max(value)) %>% 
  arrange(desc(mx), desc(value)) 
coef3.1.1 <- dplyr::select(coef3.1,-mx)
coef3.1.1$variable <- factor(coef3.1.1$variable,levels=rev(unique(coef3.1.1$variable)))


# assign group of Demo, Biomarker or Imaging
group3 <- ifelse(coef3.1.1$variable=='age','demographic',
                 ifelse(coef3.1.1$variable=='gender','demographic',
                        ifelse(coef3.1.1$variable=='race4','demographic',
                               ifelse(coef3.1.1$variable=='educ','demographic',
                                      ifelse(coef3.1.1$variable=='visc','demographic',
                                             ifelse(coef3.1.1$variable=='tkv','imaging',
                                                    ifelse(coef3.1.1$variable=='agehtn','health',
                                                           ifelse(coef3.1.1$variable=='systole','health',
                                                                  ifelse(coef3.1.1$variable=='map','health',
                                                                         ifelse(coef3.1.1$variable=='diastol','health',
                                                                                'biomarker'))))))))))
coef3.1.1$group <- group3

# Plot
Lasso_Coef3 <- ggplot(coef3.1.1, aes(variable,value,group=factor(rowid))) + 
  geom_jitter(aes(color = group)) + 
  geom_hline(yintercept=0) +
  ggtitle("Effect size of visit 2 health, demographic, biomarker and genetic predictors for c") +
  ylab("Beta Coefficients") +
  theme(plot.title=element_text(size=15,face="bold"),text=element_text(size=13),axis.text.x=element_text(angle=90))
ggsave("/Users/mengjiawei/Desktop/PKD/Lasso_Coef3o.png",scale=1,width=24,height=16,units="cm")


#visit3
ddd4 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,yrs,cpval))
ddd4 <- data.table(ddd4, key="pkdid")
ddd44 <- ddd4[,.SD[4L],by=pkdid]  
# ddd33 is the visit 3 biomarkers dataset

ddd444<-as.data.frame(ddd44)
sort(rowSums(is.na(ddd444)))
sort(colSums(is.na(ddd444))) ## checked no missing values, this is the final data for use
# note: have two missing rows for biomarkers

ddd444 <- na.omit(ddd444) ## delete the missing values

coef_lasso4 <- list()

for (i in 1:10)
{
  train_ind <- seq((i-1)*23+1,i*23,1)
  train <- ddd444[-train_ind,]
  x <- train[,!(names(train) %in% c("c","pkdid","z"))]
  y <- train[,c("c")]
  sl_lasso <- SuperLearner(Y=y,X=x,SL.library="SL.glmnet")
  coef_lasso_all <- as.matrix(sl_lasso$fitLibrary$SL.glmnet_All$object$glmnet.fit$beta)
  coef_lasso4[[i]] <- as.matrix(apply(coef_lasso_all,1,median))
}

coef4 <- do.call(cbind, coef_lasso4)
coef4.1 <- as.data.frame(t(coef4[which(abs(rowSums(coef4))>0.4),]))
coef4.1 <- melt(coef4.1) 
coef4.1$rowid <- 1:10

# sort variables by max coef value
coef4.1 <- coef4.1 %>% group_by(variable) %>% 
  mutate(mx = max(value)) %>% 
  arrange(desc(mx), desc(value)) 
coef4.1.1 <- dplyr::select(coef4.1,-mx)
coef4.1.1$variable <- factor(coef4.1.1$variable,levels=rev(unique(coef4.1.1$variable)))


# assign group of Demo, Biomarker or Imaging
group4 <- ifelse(coef4.1.1$variable=='age','demographic',
                 ifelse(coef4.1.1$variable=='gender','demographic',
                        ifelse(coef4.1.1$variable=='race4','demographic',
                               ifelse(coef4.1.1$variable=='educ','demographic',
                                      ifelse(coef4.1.1$variable=='visc','demographic',
                                             ifelse(coef4.1.1$variable=='tkv','imaging',
                                                    ifelse(coef4.1.1$variable=='agehtn','health',
                                                           ifelse(coef4.1.1$variable=='systole','health',
                                                                  ifelse(coef4.1.1$variable=='map','health',
                                                                         ifelse(coef4.1.1$variable=='diastol','health',
                                                                                'biomarker'))))))))))
coef4.1.1$group <- group4

# Plot
Lasso_Coef4 <- ggplot(coef4.1.1, aes(variable,value,group=factor(rowid))) + 
  geom_jitter(aes(color = group)) + 
  geom_hline(yintercept=0) +
  ggtitle("Effect size of visit 3 health, demographic, biomarker and genetic predictors for c") +
  ylab("Beta Coefficients") +
  theme(plot.title=element_text(size=15,face="bold"),text=element_text(size=13),axis.text.x=element_text(angle=90))
ggsave("/Users/mengjiawei/Desktop/PKD/Lasso_Coef4o.png",scale=1,width=24,height=16,units="cm")



######################## using all three year biomarkers #################
#derive 3 year data
ddd555
ddd6 <- subset(ddd555, select = -c(basedate,dvdate,birthdate,ckd_epi,cpval))
ddd6 <- data.table(ddd6, key="pkdid")
ddd66 <- ddd6[,head(.SD,3),by=pkdid]  #get 3 year data for all patients
head(ddd66)

#each patient should have at least 3 visit records
check <- as.data.frame(table(ddd66$yrs))
del <- as.matrix(check[check$Freq < 120,])
ddd66 <- ddd66[!(ddd66$yrs %in% del[,1]),]
check2 <- as.data.frame(table(ddd66$pkdid))
del2 <- as.matrix(check2[check2$Freq < 3,])
ddd66 <- ddd66[!(ddd66$pkdid %in% del2[,1]),]

## checked no missing values
sort(rowSums(is.na(ddd66)))
sort(colSums(is.na(ddd66)))
dim(ddd66)
length(unique(ddd66$pkdid))
#change the long form to wide form
ddd66<- as.data.frame(ddd66)

ddd777<-reshape(ddd66,direction = "wide",idvar = "pkdid",timevar = "yrs")

ddd777_2 <- ddd777 %>% select_if(colSums(is.na(ddd777))<10)
ddd777_2
final4 <- subset(ddd777_2, select = -c(`c.1`,`c.2`,`z.1`,
                                       `z.2`,`z.3`,`educ.1`,`educ.2`,`gender.1`,`gender.2`,`race4.1`,`race4.2`,`genetype.1`,`genetype.2`))
sort(rowSums(is.na)(final4))
sort(colSums(is.na(final4))) 

final4$`c.3`<-as.numeric(final4$`c.3`)
str(final4)
head(final4)

coef_lasso5 <- list()

for (i in 1:10)
{
  train_ind <- seq((i-1)*23+1,i*23,1)
  train <- final4[-train_ind,]
  x <- train[,!(names(train) %in% c("c.3","pkdid"))]
  y <- train[,c("c.3")]
  sl_lasso <- SuperLearner(Y=y,X=x,SL.library="SL.glmnet")
  coef_lasso_all <- as.matrix(sl_lasso$fitLibrary$SL.glmnet_All$object$glmnet.fit$beta)
  coef_lasso5[[i]] <- as.matrix(apply(coef_lasso_all,1,median))
}

coef5 <- do.call(cbind, coef_lasso5)
coef5.1 <- as.data.frame(t(coef5[which(abs(rowSums(coef5))>0.4),]))
coef5.1 <- melt(coef5.1) 
coef5.1$rowid <- 1:10
coef5.1
# sort variables by max coef value
coef5.1 <- coef5.1 %>% group_by(variable) %>% 
  mutate(mx = max(value)) %>% 
  arrange(desc(mx), desc(value)) 
coef5.1.1 <- dplyr::select(coef5.1,-mx)
coef5.1.1$variable <- factor(coef5.1.1$variable,levels=rev(unique(coef5.1.1$variable)))


# assign group of Demo, Biomarker or Imaging
group5 <- ifelse(coef5.1.1$variable=='age.1'&coef5.1.1$variable=='age.2'&coef5.1.1$variable=='age.3' ,'demographic',
                 ifelse(coef5.1.1$variable=='gender.3','demographic',
                        ifelse(coef5.1.1$variable=='race4.3','demographic',
                               ifelse(coef5.1.1$variable=='educ.3','demographic',
                                      ifelse(coef5.1.1$variable=='visc.1'&coef5.1.1$variable=='visc.2'&coef5.1.1$variable=='visc.3','demographic',
                                             ifelse(coef5.1.1$variable=='tkv.1'&coef5.1.1$variable=='tkv.2'&coef5.1.1$variable=='tkv.3','imaging',
                                                    ifelse(coef5.1.1$variable=='agehtn.1'&coef5.1.1$variable=='agehtn.2'&coef5.1.1$variable=='agehtn.3','health',
                                                           ifelse(coef5.1.1$variable=='systole.1'&coef5.1.1$variable=='systole.2'&coef5.1.1$variable=='systole.3','health',
                                                                  ifelse(coef5.1.1$variable=='map.1'&coef5.1.1$variable=='map.2'&coef5.1.1$variable=='map.3','health',
                                                                         ifelse(coef5.1.1$variable=='diastol.1'&coef5.1.1$variable=='diastol.2'&coef5.1.1$variable=='diastol.3','health',
                                                                                'biomarker'))))))))))
coef5.1.1$group <- group5

# Plot
Lasso_Coef5 <- ggplot(coef5.1.1, aes(variable,value,group=factor(rowid))) + 
  geom_jitter(aes(color = group)) + 
  geom_hline(yintercept=0) +
  ggtitle("Effect size of visit all three years health, demographic, biomarker and genetic predictors for c") +
  ylab("Beta Coefficients") +
  theme(plot.title=element_text(size=15,face="bold"),text=element_text(size=13),axis.text.x=element_text(angle=90))
ggsave("/Users/mengjiawei/Desktop/PKD/Lasso_Coef5o.png",scale=1,width=24,height=16,units="cm")






