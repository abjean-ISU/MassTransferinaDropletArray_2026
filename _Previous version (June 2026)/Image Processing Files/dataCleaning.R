rm(list=ls())

## CHANGE FILES ##
#import the relevant high speed image csvs

folder = 'C:\\Users\\anbarron\\High Speed Camera Files\\09242024\\Droplet Data Checks'

files = list.files(path=folder,full.names=T)
fileNames = list.files(path=folder)

dropletData = lapply(files,read.csv)

#loop for variables 
designation=character()
frameNo = integer()
xDiameter = numeric()
yDiameter = numeric()

for ( i in 1:length(dropletData)){
  df = data.frame(dropletData[i])
  df = df[order(df$Droplet.No,df$Frame.No),]
  des.list=matrix(data=NA,nrow=sum(df$Droplet.No!=0),ncol=1)
  frameNo.list=matrix(data=NA,nrow=sum(df$Droplet.No!=0),ncol=1)
  xDiameter.list=matrix(data=NA,nrow=sum(df$Droplet.No!=0),ncol=1)
  yDiameter.list=matrix(data=NA,nrow=sum(df$Droplet.No!=0),ncol=1)
  counter = 1
  for (j in 1:length(df$Frame.No)){
      des.list[counter] = paste(fileNames[i],as.character(df$Droplet.No[j]),sep='.')
      frameNo.list[counter] = df$Frame.No[j]
      xDiameter.list[counter] = df$X.Diameter..mm.[j]
      yDiameter.list[counter] = df$Y.Diameter..mm.[j]
      counter = counter+1
  }
  designation = c(designation,des.list)
  frameNo = c(frameNo,frameNo.list)
  xDiameter=c(xDiameter,xDiameter.list)
  yDiameter=c(yDiameter,yDiameter.list)
}

compositeData=data.frame(Droplet.Designation=designation,FrameNo=frameNo,xDiameter.mm=xDiameter,yDiameter.mm=yDiameter)

#loop for droplet no 
dropletNo = 1
counter = 1

for (i in 2:length(compositeData$FrameNo)){
  if (isTRUE(compositeData$Droplet.Designation[i]==compositeData$Droplet.Designation[i-1])){
    dropletNo = c(dropletNo,counter)
  } else {
    counter = counter+1
    dropletNo = c(dropletNo,counter)
  }
}

compositeData$Droplet.No = dropletNo

write.csv(compositeData,'C:\\Users\\anbarron\\High Speed Camera Files\\09242024\\AllDropletData.csv')
