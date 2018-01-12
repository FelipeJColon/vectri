#
# plot to make timeseries of historical stations
#
require(ncdf)

# loop over population

dir="/scratch/tompkins/kericho/output"
varname="PRd" # variable to plot
nens=57       # number of ensemble members

# MUST HAVE EVEN NUMBER OF QUANTILES PROBABILITIES:
quanprob<-c(0.1,0.25,0.75,0.9) # quantiles 
nquanprob<-length(quanprob)

# station loop on the outside
for (rainsource in c("trmm","arc2","fews")){
#for (station in c("GULU","MASINDI","MBARARA","TORORO","KABALE")){
#for (pop in c(10,30,100,300,1000)){


#d1<-1961
#d0<-1926
#if (station=="TORORO") d0<-1929

#
# read default run first
#
pattern=paste(rainsource,"_run0",sep="")

file<-list.files(path=dir,pattern=pattern,full.names=T)

#
# date handling
#
ymd0<-substr(file,nchar(file)-19,nchar(file)-12)
ymd1<-substr(file,nchar(file)-10,nchar(file)-3)
d0 <- substr(ymd0,1,4)
d1 <- substr(ymd1,1,4)
nyear=strtoi(d1)-strtoi(d0)+1
mmdd0 <-  paste("-",substr(ymd0,5,6),"-",substr(ymd0,7,8),sep="")
basedate=paste(d0,mmdd0,sep="")
print (basedate)

#
# default run here
#
nc<-open.ncdf(file)
var<-get.var.ncdf(nc,varname)
t2m <- get.var.ncdf(nc,"t2m")
rain <-  get.var.ncdf(nc,"rain")
time<-get.var.ncdf(nc,"time")
date<-as.Date(time,origin=basedate)
ntime<-dim(time)
close.ncdf(nc)



#
# now read the ensemble members
#
varens<-matrix(,nrow=nens,ncol=ntime)
for (iens in 1:nens){
    pattern=paste(rainsource,"_run",iens,sep="")
    file<-list.files(path=dir,pattern=pattern,full.names=T)
    nc<-open.ncdf(file)
    varens[iens,]<-get.var.ncdf(nc,varname)
    close.ncdf(nc)
}

#
# now calculate the quartile limits 
#
quan<-matrix(,nrow=nquanprob,ncol=ntime)
for (itime in 1:ntime){
quan[,itime]<-quantile(varens[,itime],probs=quanprob)
}

#
# make the plot
#
outfile=paste(dir,"/../vectri_plot_",rainsource,".pdf",sep="")
pdf(outfile)

#
# plot default:
#
ymax=signif(max(quan+0.05),digits=1)
plot(date,var,xaxt="n",ylim=c(0,ymax),xlab="Date",
ylab="Malaria Parasite Ratio (PR)"
,type='n',lwd=3) #,cey.axis=0.15)

#
# loop over paired quantiles
#
for (iquan in 1:(nquanprob/2)){
cshade=paste("grey",80-(iquan-1)*30,sep="")
print (paste("color ",cshade))
uquan=nquanprob+1-iquan
polygon(c(date, rev(date)), c(quan[iquan,],rev(quan[uquan,])),
     col = cshade, border = NA)
}

lines(date,var,lwd=3)

ticks<-as.Date(paste(c(d0:d1),mmdd0,sep=""),origin=ticks[1])
labels <- paste(c(d0:d1),substr(mmdd0,2,3),sep="")
years<-c(d0:d1)
axis(1,ticks,labels=labels, las=2,cex.axis=0.7)
#atloc <- seq(1,nyear,by=1)
#print (atloc)
#axis(1,atloc,labels=F,cex.axis=0.8)
#text(atloc,par("usr")[3] - 0.2,labels=ticks,pos=1,xpd=T)

#
# climate plots
#

filet=paste(dir,"/../t2m_",rainsource,".pdf",sep="")
print (filet)

pdf(filet)
plot(date,t2m,xaxt="n",ylim=c(14,24),xlab="Date",
ylab="ERAI Temperature (C)"
,type='l',lwd=3) #,cey.axis=0.15)
axis(1,ticks,labels=labels, las=2,cex.axis=0.7)
print (rain)

pdf(paste(dir,"/../",rainsource,".pdf",sep=""))
plot(date,rain,xaxt="n",ylim=c(0,100),xlab="Date",
ylab="Precipitation (mm/day)"
,type='l',lwd=3) #,cey.axis=0.15)
axis(1,ticks,labels=labels, las=2,cex.axis=0.7)


} # loop
#}
warnings() 

