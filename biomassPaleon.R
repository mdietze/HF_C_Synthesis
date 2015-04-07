###Code for calculating biomass estimates using ring-width data and field notes formatted for PalEON project.
###R code originally formatted for Goose Egg State Forest biomass estimates. Ring-width files employ a
###Site-Tree-Core identifier of SSSTTTC (nchar = 7) for each series (e.g., GE1001n = Goose Egg Plot 1, Tree 1, Core N).
###Field note formatting (field.data) may vary between sites.

# Load packages (not all required by potentially useful for data evaluation at user's convenience)
library(dplR)
library(sp)
library(foreign)
library(lattice)
library(ggplot2)
library(sciplot)
library(data.table)

#Set working directory to user's preference, read in field data*****
setwd("C:\\Users\\Daniel\\Documents\\PalEON\\GooseEgg\\")
field.data<-read.table("GooseEggAllPlots.csv", header=T,sep=",",skip=7,colClasses=c("character","integer","character","character","numeric","numeric","numeric","character","character","factor","character"))

#Read in data frame with coefficients (a & b) for alometric equations of the form a*DBH^b
#Set to user's preference****
bm.eqs<-read.table("C:\\Users\\Daniel\\Documents\\PalEON\\biomass_coeff.csv",header=T,sep=",")

#Function for current biomass
#dbh.data is a data.frame with dbh data
#Read in ring-width (.rwl) files from current working directory
files_to_read<-list.files(pattern=glob2rx("*.rwl"))
#Initialize output
core.means.all<-NULL
for (i in files_to_read){
  file_read<-read.rwl(i)
  assign(i,file_read)
  #Retrieve & "recode" series ID (SSSTTTC)
  cores<-data.frame("code"=names(file_read))
  #Manipulate string to return tree ID (SSSTTT; nchar=6)
  #Update to user's preference based on site used*****
  cores$recode<- substr(as.character(cores$code),1,6)
  core.means<-data.frame("year"=row.names(file_read))
  core.means<-NULL
  n<-1
  file_read$year<-row.names(file_read)
  for(j in unique(cores$recode)){
    print(j)
    print(n)
    tree1<-data.frame(file_read[,c(as.character(cores$code[cores$recode==j]),"year")])
    if(dim(tree1)[2]>1){
      tree2<-data.frame(rowMeans(subset(tree1, select=-c(year)),na.rm=T))
    }else{
      tree2=tree1
    }
    tree2<-na.omit(tree2)
    tree2$tree<-j
    tree2$year<-as.numeric(as.character(row.names(tree2)))
    tree2$rACUM<-0
    tree2$rACUM.inv<-0
    for(k in 1:dim(tree2)[1]){
      tree2$rACUM[k]<-sum(tree2[c(1:k),1])
      tree2$rACUM.inv[k-1]<-sum(tree2[c((k):(dim(tree2)[1])),1])
    }
    names(tree2)<-c("meanRW","recode","year","rACUM","rACUM.inv")
    if(n==1)(core.means<-tree2)else(core.means<-rbind(core.means,tree2))
    n=n+1
  }
  core.means.all<-rbind(core.means.all,core.means)
}
core.means.all$Site<-substr(as.character(core.means.all$recode),1,3)
#dim(unique(subset(core.means.all, select=tree)))

core.means.all$Tree.Number<-as.numeric(as.character(substr(as.character(core.means.all$recode),4,6)))

ab.all<-merge(core.means.all,field.data,by=c("Site","Tree.Number"))
bm.eqs.sp<-bm.eqs[bm.eqs$eq==2,]
ab.all<-merge(ab.all,bm.eqs.sp,by="Species")

#dbhest1==> dbh at the BEGINNING of the growing season
ab.all$dbhest1<-ab.all$DBH-0.2*(ab.all$rACUM.inv+ab.all$meanRW)
ab.all$dbhest1[ab.all$dbhest1<0]<-NA
#dbhest2==> dbh at the END of the growing season
ab.all$dbhest2<-ab.all$DBH-0.2*(ab.all$rACUM.inv)
ab.all$dbhest2[ab.all$dbhest2<0]<-NA
ab.all$AB<-ab.all$a*ab.all$dbhest2^ab.all$b

#annAB==>biomass acculumated during the year.
# that is, biomass at the END of the growing season
# minus the biomass at the BEGINNING of the growing season
ab.all$annAB<-ab.all$a*(ab.all$dbhest2^ab.all$b-ab.all$dbhest1^ab.all$b)

sum.fn<-function(x) sum(x, na.rm=TRUE)
count.fn<-function(x) length(unique(x, na.rm=TRUE))

# Function to calculate total biomass at each plot (non-spatially averaged)
bioAll.fn<-function(x){
  convHA <- 10000/(pi*(c(13,20,30)^2))
  nestDBH <- c(10,20,30)
  nestDist <- c(13,20,30)
  tree <- rep(sort(unique(x$Group.2)),3)
  sites <- rep(unique(x$Group.1),each=length(tree)/3)
  plotBio <- data.frame(sites,tree)
  for(i in 1:3){
    plotHec <- aggregate(x$AB[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],by=list(x$Group.1[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],x$Group.2[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]]),sum.fn)
    plotHec[,3] <- plotHec[,3]*convHA[i]
    names(plotHec) <- c("sites","tree",paste(c("AB.ha",i),collapse=""))
    plotBio <- merge(plotBio,plotHec,by=intersect(names(plotBio),names(plotHec)),all.x=TRUE)
  }
  return(plotBio)
}

totalBio.all <- aggregate(ab.all[,c("DBH","a","b","Distance")],by=list(ab.all$Site,ab.all$Tree.Number),mean)
totalBio.all$AB <- totalBio.all$a*totalBio.all$DBH^totalBio.all$b
totalBio.nest <- bioAll.fn(totalBio.all)
totalBio.nest[,3:5] <- totalBio.nest[,3:5]/1000
#mean(colSums(totalBio.nest[totalBio.nest$sites=="GE1",3:5],na.rm=TRUE),na.rm=TRUE)


#Function to calculate biomass/hectare for each subplot (nested plot)
bio.fn<-function(x){
  convHA <- 10000/(pi*(c(13,20,30)^2))
  nestDBH <- c(10,20,30)
  nestDist <- c(13,20,30)
  years <- rep(sort(unique(x$year)),3)
  sites <- rep(unique(x$Site),each=length(years)/3)
  plotBio <- data.frame(years,sites)
  for(i in 1:3){
    plotHec <- aggregate(x$annAB[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],by=list(x$year[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],x$Site[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]]),sum.fn)
    plotHec[,3] <- plotHec[,3]*convHA[i]
    names(plotHec) <- c('years',"sites",paste(c("annAB.ha",i),collapse=""))
    plotBio <- merge(plotBio,plotHec,by=intersect(names(plotBio),names(plotHec)),all.x=TRUE)
  }
  return(plotBio)
}

#Aboveground biomass per year per hectare for each nested plot
ab.ha.nest<-bio.fn(ab.all)
ab.ha.nest[,3:5] <- ab.ha.nest[,3:5]/1000

#Create function to calculate biomass/hectare for each species
bio.sp<-function(x){
  convHA <- 10000/(pi*(c(13,20,30)^2))
  nestDBH <- c(10,20,30)
  nestDist <- c(13,20,30)
  years <- rep(sort(unique(x$year)),length(unique(x$Species))*length(unique(x$Site)))
  sites <- rep(unique(x$Site),each=length(years)/3)
  species <- rep(rep(unique(x$Species),each=length(years)/(3*length(unique(x$Species)))),3)
  plotBio <- data.frame(years,sites,species)
  for(i in 1:3){
    plotHec <- aggregate(x$annAB[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],by=list(x$year[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],x$Site[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],x$Species[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]]),sum.fn)
    plotHec[,4] <- plotHec[,4]*convHA[i]
    names(plotHec) <- c('years',"sites","species",paste(c("annAB.ha",i),collapse=""))
    plotBio <- merge(plotBio,plotHec,by=intersect(names(plotBio),names(plotHec)),all.x=TRUE)
  }
  #plotBio <- data.frame(plotBio[,1:3],rowMeans(plotBio[4:6],na.rm=TRUE))
  #names(plotBio)[4] <- "annAB.ha"
  return(plotBio)
}

#Aboveground biomass per year per hectare for each species
ab.ha.spec<-bio.sp(ab.all)
ab.ha.spec[,4:6] <- ab.ha.spec[,4:6]/1000


#Create function to calculate biomass/hectare for each canopy position
bio.can<-function(x){
  convHA <- 10000/(pi*(c(13,20,30)^2))
  nestDBH <- c(10,20,30)
  nestDist <- c(13,20,30)
  years <- rep(sort(unique(x$year)),length(unique(x$Canopy[!is.na(x$Canopy)]))*length(unique(x$Site)))
  sites <- rep(unique(x$Site),each=length(years)/3)
  canopy <- rep(rep(unique(x$Canopy[!is.na(x$Canopy)]),each=length(years)/(3*length(unique(x$Canopy[!is.na(x$Canopy)])))),3)
  plotBio <- data.frame(years,sites,canopy)
  for(i in 1:3){
    plotHec <- aggregate(x$annAB[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],by=list(x$year[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],x$Site[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]],x$Canopy[x$DBH>=nestDBH[i] & x$Distance<=nestDist[i]]),sum.fn)
    plotHec[,4] <- plotHec[,4]*convHA[i]
    names(plotHec) <- c('years',"sites","canopy",paste(c("annAB.ha",i),collapse=""))
    plotBio <- merge(plotBio,plotHec,by=intersect(names(plotBio),names(plotHec)),all.x=TRUE)
  }
  #plotBio <- data.frame(plotBio[,1:3],rowMeans(plotBio[4:6],na.rm=TRUE))
  #names(plotBio)[4] <- "annAB.ha"
  return(plotBio)
}

#Aboveground biomass per year per hectare for each species
ab.ha.can<-bio.can(ab.all)
ab.ha.can[,4:6] <- ab.ha.can[,4:6]/1000


#Density
n.trees<-data.frame(tapply(X=ab.all$Tree.Number,INDEX=list(ab.all$year,ab.all$Site),count.fn))
n.trees[is.na(n.trees)]=0
#windows()
#plot(as.numeric(row.names(n.trees)),rowSums(n.trees),type="l",ylab="Sample Depth",xlab="Year")

#####END OF BIOMASS CODE



#################################
########## some Plots ###########
#################################

#Plot each nested plot's annual biomass/hectare
for(j in unique(ab.ha.nest$sites)){
  colPlot <- c("black","blue","red")
  windows(10,10)
  plot(ab.ha.nest$years[ab.ha.nest$sites==j & ab.ha.nest$years %in% c(1960:2013)],ab.ha.nest[ab.ha.nest$sites==j & ab.ha.nest$years %in% c(1960:2013),1+2],type="l",col=colPlot[1],ylab="Annual Aboveground Biomass (Mg/hectare/year)",xlab="Year",main=j,ylim=c(0,6))
  for(i in 2:3){
    lines(ab.ha.nest$years[ab.ha.nest$sites==j & ab.ha.nest$years %in% c(1960:2013)],ab.ha.nest[ab.ha.nest$sites==j & ab.ha.nest$years %in% c(1960:2013),i+2],col=colPlot[i])
  }
  legend("topleft",legend=c("13-m (>10cm)","20-m (>20cm)","30-m (>30cm)"),lwd=c(2,2,2),col=colPlot,bty="n")
}

#Plot site-level biomass
windows(10,10)
par(mar=c(5,4,4,5))
ab.ha.nest[is.na(ab.ha.nest)] <- 0
plot(aggregate(ab.ha.nest$years[ab.ha.nest$years %in% c(1960:2013)],list(ab.ha.nest$years[ab.ha.nest$years %in% c(1960:2013)]),unique)[,2],aggregate(rowMeans(ab.ha.nest[ab.ha.nest$years %in% c(1960:2013),3:5],na.rm=TRUE),list(ab.ha.nest$years[ab.ha.nest$years %in% c(1960:2013)]),mean)[,2],lwd=3,ylab="Annual Aboveground Biomass (Mg/hectare/year)",xlab="Year",main="Goose Egg",ylim=c(0,6),type="l")

