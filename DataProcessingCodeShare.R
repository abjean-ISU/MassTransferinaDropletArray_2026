#Liquid Holdup/Experimental results data processing
rm(list=ls())

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


#Importing Mass Transfer Experimental Data 
ExperimentalData.Filepath = 'Experimental Data/MassTransfer Experimental Data 10222025.xlsx'
ExperimentalSheets = excel_sheets(ExperimentalData.Filepath)

ExperimentalData.df = lapply(ExperimentalSheets,function(sheet){
  read_excel(ExperimentalData.Filepath,sheet=sheet)
})
names(ExperimentalData.df) = ExperimentalSheets

CleanData.df = ExperimentalData.df


#Identifying outliers from each trial of the data set using Modified Z Score 
for (i in 1:length(ExperimentalData.df)){
  Medians = (tapply(ExperimentalData.df[[i]]$`CO2 (g/L)`,ExperimentalData.df[[i]]$`Trial`,FUN=median))
  
  medDiff = vector()
  meds = vector()
  for (j in 1:length(ExperimentalData.df[[i]]$Trial)){
    TrialNo = ExperimentalData.df[[i]]$Trial[[j]]
    meds[j] = Medians[TrialNo]
    medDiff[j] = abs(ExperimentalData.df[[i]]$`CO2 (g/L)`[[j]]-Medians[TrialNo])
  }
  
  ExperimentalData.df[[i]]$median = meds
  ExperimentalData.df[[i]]$medianDiff = medDiff
  mads = tapply(ExperimentalData.df[[i]]$medianDiff,ExperimentalData.df[[i]]$`Trial`,FUN=median)
  
  mad = vector()
  ZScore = vector()
  for (j in 1:length(ExperimentalData.df[[i]]$Trial)){
    TrialNo = ExperimentalData.df[[i]]$Trial[[j]]
    mad[j] = mads[TrialNo]
  }
  
  ExperimentalData.df[[i]]$mad = mad
  ExperimentalData.df[[i]]$mad[ExperimentalData.df[[i]]$mad==0] = 0.005*1.96   #floor at half a display step
  ExperimentalData.df[[i]]$ZScore = 0.6745*(ExperimentalData.df[[i]]$`CO2 (g/L)`-ExperimentalData.df[[i]]$median)/ExperimentalData.df[[i]]$mad
  ExperimentalData.df[[i]] = replace(ExperimentalData.df[[i]],is.na(ExperimentalData.df[[i]]),0)
  
  ExperimentalData.df[[i]]$outlier = ExperimentalData.df[[i]]$ZScore<(-3.5)|ExperimentalData.df[[i]]$ZScore>(3.5)
  
}

#Creating data frame of all trials
CO2_ExperimentalResults.df = data.frame(FlowRate.LPM = rep(c(3,4.5),each=9),
                                        Height.m = rep(c(rep(c(0.53,0.68,0.98),each=3)),2),
                                        Trial = rep(c(1,2,3),6))

meanCO2.gpL = vector()
sdCO2.gpL = vector()

meanT.C = vector()
sdT.C = vector()


#Determine mean and error in each trial at each height 
for (i in 1:length(CleanData.df)){
  mean.new.CO2 = tapply(CleanData.df[[i]]$`CO2 (g/L)`,CleanData.df[[i]]$`Trial`,FUN=mean)
  sd.new.CO2 = tapply(CleanData.df[[i]]$`CO2 (g/L)`,CleanData.df[[i]]$`Trial`,FUN=sd)

  mean.new.T = tapply(CleanData.df[[i]]$`Liquid Temperature (C)`,CleanData.df[[i]]$`Trial`,FUN=mean)
  sd.new.T = tapply(CleanData.df[[i]]$`Liquid Temperature (C)`,CleanData.df[[i]]$`Trial`,FUN=sd)
  
  meanCO2.gpL = c(meanCO2.gpL,mean.new.CO2)
  sdCO2.gpL = c(sdCO2.gpL,sd.new.CO2)
  
  meanT.C = c(meanT.C,mean.new.T)
  sdT.C = c(sdT.C,sd.new.T)
}

#Instrument uncertainties 
CO2.resolution = 0.005*1.96/sqrt(3) #display reads to 0.01 vol/vol
CO2.accuracy = 0.05/sqrt(3) #Carbo 5100 standard accuracy (+/- 0.05 g/L)
T.resolution = 0.005/sqrt (3) #display reads to 0.01 C

CO2_ExperimentalResults.df$CO2.gpL = meanCO2.gpL
CO2_ExperimentalResults.df$CO2.SD.gpL = sqrt(sdCO2.gpL^2+CO2.resolution^2) #scatter within run
CO2_ExperimentalResults.df$CO2.u.gpL = sqrt(sdCO2.gpL^2+CO2.resolution^2+CO2.accuracy^2) #total uncertainty for z-score
CO2_ExperimentalResults.df$Temp.C = meanT.C
CO2_ExperimentalResults.df$Temp.SD.C = sqrt(sdT.C^2+T.resolution^2)

#Determining Modified Z Scores   
CO2_ExperimentalResults.df$median = rep(as.vector(tapply(CO2_ExperimentalResults.df$CO2.gpL,
                                                     list(CO2_ExperimentalResults.df$Height.m,CO2_ExperimentalResults.df$FlowRate.LPM),
                                                     FUN=median)),each=3)
CO2_ExperimentalResults.df$medianDiff = abs(CO2_ExperimentalResults.df$CO2.gpL-CO2_ExperimentalResults.df$median)
CO2_ExperimentalResults.df$mad = rep(as.vector(tapply(CO2_ExperimentalResults.df$medianDiff,
                                                          list(CO2_ExperimentalResults.df$Height.m,CO2_ExperimentalResults.df$FlowRate.LPM),
                                                          FUN=median, constant=1)),each=3)

CO2_ExperimentalResults.df$ZScore = 0.6745*(CO2_ExperimentalResults.df$CO2.gpL-CO2_ExperimentalResults.df$median)/CO2_ExperimentalResults.df$mad

CO2_ExperimentalResults.df$outlier = CO2_ExperimentalResults.df$ZScore<(-3.5)|CO2_ExperimentalResults.df$ZScore>(3.5)


#Creating data subsets
ExperimentalResults_4.5LPM = CO2_ExperimentalResults.df[CO2_ExperimentalResults.df$FlowRate.LPM == 4.5,]
ExperimentalResults_3.0LPM = CO2_ExperimentalResults.df[CO2_ExperimentalResults.df$FlowRate.LPM == 3.0,]


#Importing matlab modeling results 
High.3LPM.T1.ModelResults = read.csv('Model Results/3LPM CO2 High 10272026/T1/MT_Models.csv')
High.3LPM.T2.ModelResults = read.csv('Model Results/3LPM CO2 High 10272026/T2/MT_Models.csv')
High.3LPM.T3.ModelResults = read.csv('Model Results/3LPM CO2 High 10272026/T3/MT_Models.csv')

Med.3LPM.T1.ModelResults = read.csv('Model Results/3LPM CO2 Med 10272026/T1/MT_Models.csv')
Med.3LPM.T2.ModelResults = read.csv('Model Results/3LPM CO2 Med 10272026/T2/MT_Models.csv')
Med.3LPM.T3.ModelResults = read.csv('Model Results/3LPM CO2 Med 10272026/T3/MT_Models.csv')

Low.3LPM.T1.ModelResults = read.csv('Model Results/3LPM CO2 Low 10272026/T1/MT_Models.csv')
Low.3LPM.T2.ModelResults = read.csv('Model Results/3LPM CO2 Low 10272026/T2/MT_Models.csv')
Low.3LPM.T3.ModelResults = read.csv('Model Results/3LPM CO2 Low 10272026/T3/MT_Models.csv')


High.4.5LPM.T1.ModelResults = read.csv('Model Results/4.5LPM CO2 High 10272026/T1/MT_Models.csv')
High.4.5LPM.T2.ModelResults = read.csv('Model Results/4.5LPM CO2 High 10272026/T2/MT_Models.csv')
High.4.5LPM.T3.ModelResults = read.csv('Model Results/4.5LPM CO2 High 10272026/T3/MT_Models.csv')

Med.4.5LPM.T1.ModelResults = read.csv('Model Results/4.5LPM CO2 Med 10272026/T1/MT_Models.csv')
Med.4.5LPM.T2.ModelResults = read.csv('Model Results/4.5LPM CO2 Med 10272026/T2/MT_Models.csv')
Med.4.5LPM.T3.ModelResults = read.csv('Model Results/4.5LPM CO2 Med 10272026/T3/MT_Models.csv')

Low.4.5LPM.T1.ModelResults = read.csv('Model Results/4.5LPM CO2 Low 10272026/T1/MT_Models.csv')
Low.4.5LPM.T2.ModelResults = read.csv('Model Results/4.5LPM CO2 Low 10272026/T2/MT_Models.csv')
Low.4.5LPM.T3.ModelResults = read.csv('Model Results/4.5LPM CO2 Low 10272026/T3/MT_Models.csv')

#Finding the minimum and maximum temperatures for each flowrate 
ID_minT_3.0 = which.min(ExperimentalResults_3.0LPM$Temp.C)
MinT_3.0 = ExperimentalResults_3.0LPM$Temp.C[ID_minT_3.0]
if (ExperimentalResults_3.0LPM$Height.m[ID_minT_3.0]==0.98){
  if(ExperimentalResults_3.0LPM$Trial[ID_minT_3.0]==1){
    MinT3.0_dataset = High.3LPM.T1.ModelResults
  }else if (ExperimentalResults_3.0LPM$Trial[ID_minT_3.0]==2){
    MinT3.0_dataset = High.3LPM.T2.ModelResults
  }else{
    MinT3.0_dataset = High.3LPM.T3.ModelResults
  }
}else if (ExperimentalResults_3.0LPM$Height.m[ID_minT_3.0]==0.68){
  if(ExperimentalResults_3.0LPM$Trial[ID_minT_3.0]==1){
    MinT3.0_dataset = Med.3LPM.T1.ModelResults
  }else if (ExperimentalResults_3.0LPM$Trial[ID_minT_3.0]==2){
    MinT3.0_dataset = Med.3LPM.T2.ModelResults
  }else{
    MinT3.0_dataset = Med.3LPM.T3.ModelResults
  }
}else {
  if(ExperimentalResults_3.0LPM$Trial[ID_minT_3.0]==1){
    MinT3.0_dataset = Low.3LPM.T1.ModelResults
  }else if (ExperimentalResults_3.0LPM$Trial[ID_minT_3.0]==2){
    MinT3.0_dataset = Low.3LPM.T2.ModelResults
  }else{
    MinT3.0_dataset = Low.3LPM.T3.ModelResults
  }
}

ID_maxT_3.0 = which.max(ExperimentalResults_3.0LPM$Temp.C)
MaxT_3.0 = ExperimentalResults_3.0LPM$Temp.C[ID_maxT_3.0]
if (ExperimentalResults_3.0LPM$Height.m[ID_maxT_3.0]==0.98){
  if(ExperimentalResults_3.0LPM$Trial[ID_maxT_3.0]==1){
    MaxT3.0_dataset = High.3LPM.T1.ModelResults
  }else if (ExperimentalResults_3.0LPM$Trial[ID_maxT_3.0]==2){
    MaxT3.0_dataset = High.3LPM.T2.ModelResults
  }else{
    MaxT3.0_dataset = High.3LPM.T3.ModelResults
  }
}else if (ExperimentalResults_3.0LPM$Height.m[ID_maxT_3.0]==0.68){
  if(ExperimentalResults_3.0LPM$Trial[ID_maxT_3.0]==1){
    MaxT3.0_dataset = Med.3LPM.T1.ModelResults
  }else if (ExperimentalResults_3.0LPM$Trial[ID_maxT_3.0]==2){
    MaxT3.0_dataset = Med.3LPM.T2.ModelResults
  }else{
    MaxT3.0_dataset = Med.3LPM.T3.ModelResults
  }
}else {
  if(ExperimentalResults_3.0LPM$Trial[ID_maxT_3.0]==1){
    MaxT3.0_dataset = Low.3LPM.T1.ModelResults
  }else if (ExperimentalResults_3.0LPM$Trial[ID_maxT_3.0]==2){
    MaxT3.0_dataset = Low.3LPM.T2.ModelResults
  }else{
    MaxT3.0_dataset = Low.3LPM.T3.ModelResults
  }
}


ID_minT_4.5 = which.min(ExperimentalResults_4.5LPM$Temp.C)
MinT_4.5 = ExperimentalResults_4.5LPM$Temp.C[ID_minT_4.5]
if (ExperimentalResults_4.5LPM$Height.m[ID_minT_4.5]==0.98){
  if(ExperimentalResults_4.5LPM$Trial[ID_minT_4.5]==1){
    MinT4.5_dataset = High.4.5LPM.T1.ModelResults
  }else if (ExperimentalResults_4.5LPM$Trial[ID_minT_4.5]==2){
    MinT4.5_dataset = High.4.5LPM.T2.ModelResults
  }else{
    MinT4.5_dataset = High.4.5LPM.T3.ModelResults
  }
}else if (ExperimentalResults_4.5LPM$Height.m[ID_minT_4.5]==0.68){
  if(ExperimentalResults_4.5LPM$Trial[ID_minT_4.5]==1){
    MinT4.5_dataset = Med.4.5LPM.T1.ModelResults
  }else if (ExperimentalResults_4.5LPM$Trial[ID_minT_4.5]==2){
    MinT4.5_dataset = Med.4.5LPM.T2.ModelResults
  }else{
    MinT4.5_dataset = Med.4.5LPM.T3.ModelResults
  }
}else {
  if(ExperimentalResults_4.5LPM$Trial[ID_minT_4.5]==1){
    MinT4.5_dataset = Low.4.5LPM.T1.ModelResults
  }else if (ExperimentalResults_4.5LPM$Trial[ID_minT_4.5]==2){
    MinT4.5_dataset = Low.4.5LPM.T2.ModelResults
  }else{
    MinT4.5_dataset = Low.4.5LPM.T3.ModelResults
  }
}

ID_maxT_4.5 = which.max(ExperimentalResults_4.5LPM$Temp.C)
MaxT_4.5 = ExperimentalResults_4.5LPM$Temp.C[ID_maxT_4.5]
if (ExperimentalResults_4.5LPM$Height.m[ID_maxT_4.5]==0.98){
  if(ExperimentalResults_4.5LPM$Trial[ID_maxT_4.5]==1){
    MaxT4.5_dataset = High.4.5LPM.T1.ModelResults
  }else if (ExperimentalResults_4.5LPM$Trial[ID_maxT_4.5]==2){
    MaxT4.5_dataset = High.4.5LPM.T2.ModelResults
  }else{
    MaxT4.5_dataset = High.4.5LPM.T3.ModelResults
  }
}else if (ExperimentalResults_4.5LPM$Height.m[ID_maxT_4.5]==0.68){
  if(ExperimentalResults_4.5LPM$Trial[ID_maxT_4.5]==1){
    MaxT4.5_dataset = Med.4.5LPM.T1.ModelResults
  }else if (ExperimentalResults_4.5LPM$Trial[ID_maxT_4.5]==2){
    MaxT4.5_dataset = Med.4.5LPM.T2.ModelResults
  }else{
    MaxT4.5_dataset = Med.4.5LPM.T3.ModelResults
  }
}else {
  if(ExperimentalResults_4.5LPM$Trial[ID_maxT_4.5]==1){
    MaxT4.5_dataset = Low.4.5LPM.T1.ModelResults
  }else if (ExperimentalResults_4.5LPM$Trial[ID_maxT_4.5]==2){
    MaxT4.5_dataset = Low.4.5LPM.T2.ModelResults
  }else{
    MaxT4.5_dataset = Low.4.5LPM.T3.ModelResults
  }
}

#Creating composite plots of all data and models  

#Creating legend aesthetics
legend_df = data.frame(x=0,y=0,label = factor(c('Ruckenstein',
                                                'Angelo',
                                                'Hsu',
                                                'Amokrane',
                                                'Experimental Data'),
                                              levels = c('Ruckenstein','Angelo','Hsu','Amokrane','Experimental Data')))
legend_color = c("Ruckenstein" = safe_colorblind_palette[11],
                 "Angelo"      = safe_colorblind_palette[2],
                 "Hsu"         = safe_colorblind_palette[4],
                 "Amokrane"    = safe_colorblind_palette[3],
                 "Experimental Data"  = "black")
legend_shape = c("Ruckenstein" = 22,"Angelo" = 22,"Hsu" = 22,
                 "Amokrane" = 22,"Experimental Data"= 16)


#3.0 LPM
compositePlot.3.0LPM = ggplot(MaxT3.0_dataset,aes(x=Height.m))+
  geom_line(aes(y=MinT3.0_dataset$Ruckenstein),color='black',linewidth=1)+
  geom_line(aes(y=MaxT3.0_dataset$Ruckenstein),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT3.0_dataset$Ruckenstein,
                  ymax=MinT3.0_dataset$Ruckenstein),fill=safe_colorblind_palette[11],alpha=0.8)+
  geom_line(aes(y=MinT3.0_dataset$Hsu),color='black',linewidth=1)+
  geom_line(aes(y=MaxT3.0_dataset$Hsu),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT3.0_dataset$Hsu,
                  ymax=MinT3.0_dataset$Hsu),fill=safe_colorblind_palette[4],alpha=0.8)+
  geom_line(aes(y=MinT3.0_dataset$Angelo),color='black',linewidth=1)+
  geom_line(aes(y=MaxT3.0_dataset$Angelo),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT3.0_dataset$Angelo,
                  ymax=MinT3.0_dataset$Angelo),fill=safe_colorblind_palette[2],alpha=0.8)+
  geom_line(aes(y=MinT3.0_dataset$Amokrane),color='black',linewidth=1)+
  geom_line(aes(y=MaxT3.0_dataset$Amokrane),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT3.0_dataset$Amokrane,
                  ymax=MinT3.0_dataset$Amokrane),fill=safe_colorblind_palette[3],alpha=0.8)+
  geom_point(data=ExperimentalResults_3.0LPM,aes(x=Height.m,y=CO2.gpL),color='black',size=3)+
  geom_errorbar(data=ExperimentalResults_3.0LPM,aes(x=Height.m,ymin=(CO2.gpL - CO2.u.gpL),
                                                      ymax = (CO2.gpL + CO2.u.gpL)),
                width=0)+
  theme_bw()+
  xlab('\nDistance from Droplet Manifold (m)')+
  ylab(bquote(paste(CO[2],' Concentration (g/L)')))+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0))+
  ylim(0,2.5)+
  theme(text = element_text(family = "sans"))+
  theme(aspect.ratio=.5)+
  theme(text=element_text(size=20))

compositePlot.3.0LPM = compositePlot.3.0LPM+
  geom_point(data=legend_df,aes(x=x,y=y,color=label,shape=label,fill=label),size=0)+
  scale_color_manual(values=legend_color,name=NULL)+
  scale_fill_manual(values=legend_color,name=NULL)+
  scale_shape_manual(values=legend_shape,name=NULL)+
  guides(color = guide_legend(override.aes = list(size=5)),
         fill = guide_legend(override.aes = list(size=5)),
         shape = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=20))+
  theme(axis.text = element_text(size = 20))+
  theme(legend.position='bottom')
  
print(compositePlot.3.0LPM)  

#4.5 LPM
compositePlot.4.5LPM = ggplot(MaxT4.5_dataset,aes(x=Height.m))+
  geom_line(aes(y=MinT4.5_dataset$Ruckenstein),color='black',linewidth=1)+
  geom_line(aes(y=MaxT4.5_dataset$Ruckenstein),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT4.5_dataset$Ruckenstein,
                  ymax=MinT4.5_dataset$Ruckenstein),fill=safe_colorblind_palette[11],alpha=0.8)+
  geom_line(aes(y=MinT4.5_dataset$Hsu),color='black',linewidth=1)+
  geom_line(aes(y=MaxT4.5_dataset$Hsu),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT4.5_dataset$Hsu,
                  ymax=MinT4.5_dataset$Hsu),fill=safe_colorblind_palette[4],alpha=0.8)+
  geom_line(aes(y=MinT4.5_dataset$Angelo),color='black',linewidth=1)+
  geom_line(aes(y=MaxT4.5_dataset$Angelo),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT4.5_dataset$Angelo,
                  ymax=MinT4.5_dataset$Angelo),fill=safe_colorblind_palette[2],alpha=0.8)+
  geom_line(aes(y=MinT4.5_dataset$Amokrane),color='black',linewidth=1)+
  geom_line(aes(y=MaxT4.5_dataset$Amokrane),color='black',linewidth=1)+
  geom_ribbon(aes(ymin=MaxT4.5_dataset$Amokrane,
                  ymax=MinT4.5_dataset$Amokrane),fill=safe_colorblind_palette[3],alpha=0.8)+
  geom_point(data=ExperimentalResults_4.5LPM,aes(x=Height.m,y=CO2.gpL),color='black',size=3)+
  geom_errorbar(data=ExperimentalResults_4.5LPM,aes(x=Height.m,ymin=(CO2.gpL - CO2.u.gpL),
                                                    ymax = (CO2.gpL + CO2.u.gpL)), width=0)+
  theme_bw()+
  xlab('\nDistance from Droplet Manifold (m)')+
  ylab(bquote(paste(CO[2],' Concentration (g/L)')))+
  scale_x_continuous(limits = c(0,1), expand = c(0, 0))+
  ylim(0,2.5)+
  theme(text = element_text(family = "sans"))+
  theme(aspect.ratio=.5)+
  theme(text=element_text(size=20))

compositePlot.4.5LPM = compositePlot.4.5LPM+
  geom_point(data=legend_df,aes(x=x,y=y,color=label,shape=label,fill=label),size=0)+
  scale_color_manual(values=legend_color,name=NULL)+
  scale_fill_manual(values=legend_color,name=NULL)+
  scale_shape_manual(values=legend_shape,name=NULL)+
  guides(color = guide_legend(override.aes = list(size=5)),
         fill = guide_legend(override.aes = list(size=5)),
         shape = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=20))+
  theme(axis.text = element_text(size = 20))+
  theme(legend.position='bottom')

print(compositePlot.4.5LPM)  

#Model uncertainty from runUncertainty.m (one row per trial)
ModelU = read.csv('Model Results/ModelUncertainty_10272026.csv')
ModelU = ModelU[order(ModelU$FlowRate.LPM,ModelU$Height.m,ModelU$Trial),]   #same row order as CO2_ExperimentalResults.df

#Determining model results for each experiment 
CO2_ExperimentalResults.df$Amokrane = c(Low.3LPM.T1.ModelResults$Amokrane[Low.3LPM.T1.ModelResults$Height.m==0.53],
                        Low.3LPM.T2.ModelResults$Amokrane[Low.3LPM.T2.ModelResults$Height.m==0.53],
                        Low.3LPM.T3.ModelResults$Amokrane[Low.3LPM.T3.ModelResults$Height.m==0.53],
                        Med.3LPM.T1.ModelResults$Amokrane[Med.3LPM.T1.ModelResults$Height.m==0.68],
                        Med.3LPM.T2.ModelResults$Amokrane[Med.3LPM.T2.ModelResults$Height.m==0.68],
                        Med.3LPM.T3.ModelResults$Amokrane[Med.3LPM.T3.ModelResults$Height.m==0.68],
                        High.3LPM.T1.ModelResults$Amokrane[High.3LPM.T1.ModelResults$Height.m==0.98],
                        High.3LPM.T2.ModelResults$Amokrane[High.3LPM.T2.ModelResults$Height.m==0.98],
                        High.3LPM.T3.ModelResults$Amokrane[High.3LPM.T3.ModelResults$Height.m==0.98],
                        Low.4.5LPM.T1.ModelResults$Amokrane[Low.4.5LPM.T1.ModelResults$Height.m==0.53],
                        Low.4.5LPM.T2.ModelResults$Amokrane[Low.4.5LPM.T2.ModelResults$Height.m==0.53],
                        Low.4.5LPM.T3.ModelResults$Amokrane[Low.4.5LPM.T3.ModelResults$Height.m==0.53],
                        Med.4.5LPM.T1.ModelResults$Amokrane[Med.4.5LPM.T1.ModelResults$Height.m==0.68],
                        Med.4.5LPM.T2.ModelResults$Amokrane[Med.4.5LPM.T2.ModelResults$Height.m==0.68],
                        Med.4.5LPM.T3.ModelResults$Amokrane[Med.4.5LPM.T3.ModelResults$Height.m==0.68],
                        High.4.5LPM.T1.ModelResults$Amokrane[High.4.5LPM.T1.ModelResults$Height.m==0.98],
                        High.4.5LPM.T2.ModelResults$Amokrane[High.4.5LPM.T2.ModelResults$Height.m==0.98],
                        High.4.5LPM.T3.ModelResults$Amokrane[High.4.5LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Amokrane.Residual = (CO2_ExperimentalResults.df$Amokrane-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Amokrane.Percenterror = abs(CO2_ExperimentalResults.df$Amokrane.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Amokrane.u.gpL = ModelU$Amokrane.u.gpL
CO2_ExperimentalResults.df$Amokrane.ZScore = CO2_ExperimentalResults.df$Amokrane.Residual/
  sqrt(CO2_ExperimentalResults.df$CO2.u.gpL^2+CO2_ExperimentalResults.df$Amokrane.u.gpL^2)

CO2_ExperimentalResults.df$Hsu = c(Low.3LPM.T1.ModelResults$Hsu[Low.3LPM.T1.ModelResults$Height.m==0.53],
                        Low.3LPM.T2.ModelResults$Hsu[Low.3LPM.T2.ModelResults$Height.m==0.53],
                        Low.3LPM.T3.ModelResults$Hsu[Low.3LPM.T3.ModelResults$Height.m==0.53],
                        Med.3LPM.T1.ModelResults$Hsu[Med.3LPM.T1.ModelResults$Height.m==0.68],
                        Med.3LPM.T2.ModelResults$Hsu[Med.3LPM.T2.ModelResults$Height.m==0.68],
                        Med.3LPM.T3.ModelResults$Hsu[Med.3LPM.T3.ModelResults$Height.m==0.68],
                        High.3LPM.T1.ModelResults$Hsu[High.3LPM.T1.ModelResults$Height.m==0.98],
                        High.3LPM.T2.ModelResults$Hsu[High.3LPM.T2.ModelResults$Height.m==0.98],
                        High.3LPM.T3.ModelResults$Hsu[High.3LPM.T3.ModelResults$Height.m==0.98],
                        Low.4.5LPM.T1.ModelResults$Hsu[Low.4.5LPM.T1.ModelResults$Height.m==0.53],
                        Low.4.5LPM.T2.ModelResults$Hsu[Low.4.5LPM.T2.ModelResults$Height.m==0.53],
                        Low.4.5LPM.T3.ModelResults$Hsu[Low.4.5LPM.T3.ModelResults$Height.m==0.53],
                        Med.4.5LPM.T1.ModelResults$Hsu[Med.4.5LPM.T1.ModelResults$Height.m==0.68],
                        Med.4.5LPM.T2.ModelResults$Hsu[Med.4.5LPM.T2.ModelResults$Height.m==0.68],
                        Med.4.5LPM.T3.ModelResults$Hsu[Med.4.5LPM.T3.ModelResults$Height.m==0.68],
                        High.4.5LPM.T1.ModelResults$Hsu[High.4.5LPM.T1.ModelResults$Height.m==0.98],
                        High.4.5LPM.T2.ModelResults$Hsu[High.4.5LPM.T2.ModelResults$Height.m==0.98],
                        High.4.5LPM.T3.ModelResults$Hsu[High.4.5LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Hsu.Residual = (CO2_ExperimentalResults.df$Hsu-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Hsu.Percenterror = abs(CO2_ExperimentalResults.df$Hsu.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Hsu.u.gpL = ModelU$Hsu.u.gpL
CO2_ExperimentalResults.df$Hsu.ZScore = CO2_ExperimentalResults.df$Hsu.Residual/
  sqrt(CO2_ExperimentalResults.df$CO2.u.gpL^2+CO2_ExperimentalResults.df$Hsu.u.gpL^2)

CO2_ExperimentalResults.df$Angelo = c(Low.3LPM.T1.ModelResults$Angelo[Low.3LPM.T1.ModelResults$Height.m==0.53],
                   Low.3LPM.T2.ModelResults$Angelo[Low.3LPM.T2.ModelResults$Height.m==0.53],
                   Low.3LPM.T3.ModelResults$Angelo[Low.3LPM.T3.ModelResults$Height.m==0.53],
                   Med.3LPM.T1.ModelResults$Angelo[Med.3LPM.T1.ModelResults$Height.m==0.68],
                   Med.3LPM.T2.ModelResults$Angelo[Med.3LPM.T2.ModelResults$Height.m==0.68],
                   Med.3LPM.T3.ModelResults$Angelo[Med.3LPM.T3.ModelResults$Height.m==0.68],
                   High.3LPM.T1.ModelResults$Angelo[High.3LPM.T1.ModelResults$Height.m==0.98],
                   High.3LPM.T2.ModelResults$Angelo[High.3LPM.T2.ModelResults$Height.m==0.98],
                   High.3LPM.T3.ModelResults$Angelo[High.3LPM.T3.ModelResults$Height.m==0.98],
                   Low.4.5LPM.T1.ModelResults$Angelo[Low.4.5LPM.T1.ModelResults$Height.m==0.53],
                   Low.4.5LPM.T2.ModelResults$Angelo[Low.4.5LPM.T2.ModelResults$Height.m==0.53],
                   Low.4.5LPM.T3.ModelResults$Angelo[Low.4.5LPM.T3.ModelResults$Height.m==0.53],
                   Med.4.5LPM.T1.ModelResults$Angelo[Med.4.5LPM.T1.ModelResults$Height.m==0.68],
                   Med.4.5LPM.T2.ModelResults$Angelo[Med.4.5LPM.T2.ModelResults$Height.m==0.68],
                   Med.4.5LPM.T3.ModelResults$Angelo[Med.4.5LPM.T3.ModelResults$Height.m==0.68],
                   High.4.5LPM.T1.ModelResults$Angelo[High.4.5LPM.T1.ModelResults$Height.m==0.98],
                   High.4.5LPM.T2.ModelResults$Angelo[High.4.5LPM.T2.ModelResults$Height.m==0.98],
                   High.4.5LPM.T3.ModelResults$Angelo[High.4.5LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Angelo.Residual = (CO2_ExperimentalResults.df$Angelo-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Angelo.Percenterror = abs(CO2_ExperimentalResults.df$Angelo.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Angelo.u.gpL = ModelU$Angelo.u.gpL
CO2_ExperimentalResults.df$Angelo.ZScore = CO2_ExperimentalResults.df$Angelo.Residual/
  sqrt(CO2_ExperimentalResults.df$CO2.u.gpL^2+CO2_ExperimentalResults.df$Angelo.u.gpL^2)

CO2_ExperimentalResults.df$Ruckenstein = c(Low.3LPM.T1.ModelResults$Ruckenstein[Low.3LPM.T1.ModelResults$Height.m==0.53],
                      Low.3LPM.T2.ModelResults$Ruckenstein[Low.3LPM.T2.ModelResults$Height.m==0.53],
                      Low.3LPM.T3.ModelResults$Ruckenstein[Low.3LPM.T3.ModelResults$Height.m==0.53],
                      Med.3LPM.T1.ModelResults$Ruckenstein[Med.3LPM.T1.ModelResults$Height.m==0.68],
                      Med.3LPM.T2.ModelResults$Ruckenstein[Med.3LPM.T2.ModelResults$Height.m==0.68],
                      Med.3LPM.T3.ModelResults$Ruckenstein[Med.3LPM.T3.ModelResults$Height.m==0.68],
                      High.3LPM.T1.ModelResults$Ruckenstein[High.3LPM.T1.ModelResults$Height.m==0.98],
                      High.3LPM.T2.ModelResults$Ruckenstein[High.3LPM.T2.ModelResults$Height.m==0.98],
                      High.3LPM.T3.ModelResults$Ruckenstein[High.3LPM.T3.ModelResults$Height.m==0.98],
                      Low.4.5LPM.T1.ModelResults$Ruckenstein[Low.4.5LPM.T1.ModelResults$Height.m==0.53],
                      Low.4.5LPM.T2.ModelResults$Ruckenstein[Low.4.5LPM.T2.ModelResults$Height.m==0.53],
                      Low.4.5LPM.T3.ModelResults$Ruckenstein[Low.4.5LPM.T3.ModelResults$Height.m==0.53],
                      Med.4.5LPM.T1.ModelResults$Ruckenstein[Med.4.5LPM.T1.ModelResults$Height.m==0.68],
                      Med.4.5LPM.T2.ModelResults$Ruckenstein[Med.4.5LPM.T2.ModelResults$Height.m==0.68],
                      Med.4.5LPM.T3.ModelResults$Ruckenstein[Med.4.5LPM.T3.ModelResults$Height.m==0.68],
                      High.4.5LPM.T1.ModelResults$Ruckenstein[High.4.5LPM.T1.ModelResults$Height.m==0.98],
                      High.4.5LPM.T2.ModelResults$Ruckenstein[High.4.5LPM.T2.ModelResults$Height.m==0.98],
                      High.4.5LPM.T3.ModelResults$Ruckenstein[High.4.5LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Ruckenstein.Residual = (CO2_ExperimentalResults.df$Ruckenstein-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Ruckenstein.Percenterror = abs(CO2_ExperimentalResults.df$Ruckenstein.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Ruckenstein.u.gpL = ModelU$Ruckenstein.u.gpL
CO2_ExperimentalResults.df$Ruckenstein.ZScore = CO2_ExperimentalResults.df$Ruckenstein.Residual/
  sqrt(CO2_ExperimentalResults.df$CO2.u.gpL^2+CO2_ExperimentalResults.df$Ruckenstein.u.gpL^2)

#creating a table to illustrate model accuracy 
ModelAccuracy.df = data.frame('Model'=c('Amokrane','Hsu','Angelo','Ruckenstein'))

ModelAccuracy.df$Bias = c(1/length(CO2_ExperimentalResults.df$Amokrane)*sum(CO2_ExperimentalResults.df$Amokrane.ZScore),
                          1/length(CO2_ExperimentalResults.df$Hsu)*sum(CO2_ExperimentalResults.df$Hsu.ZScore),
                          1/length(CO2_ExperimentalResults.df$Angelo)*sum(CO2_ExperimentalResults.df$Angelo.ZScore),
                          1/length(CO2_ExperimentalResults.df$Ruckenstein)*sum(CO2_ExperimentalResults.df$Ruckenstein.ZScore))

ModelAccuracy.df$meanZ = c(1/length(CO2_ExperimentalResults.df$Amokrane)*sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)),
                          1/length(CO2_ExperimentalResults.df$Hsu)*sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)),
                          1/length(CO2_ExperimentalResults.df$Angelo)*sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)),
                          1/length(CO2_ExperimentalResults.df$Ruckenstein)*sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)))

ModelAccuracy.df$medianz = c(median(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)),
                             median(abs(CO2_ExperimentalResults.df$Hsu.ZScore)),
                             median(abs(CO2_ExperimentalResults.df$Angelo.ZScore)),
                             median(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)))

ModelAccuracy.df$meanResidal = c(1/length(CO2_ExperimentalResults.df$Amokrane)*sum(abs(CO2_ExperimentalResults.df$Amokrane.Residual)),
                                 1/length(CO2_ExperimentalResults.df$Hsu)*sum(abs(CO2_ExperimentalResults.df$Hsu.Residual)),
                                 1/length(CO2_ExperimentalResults.df$Angelo)*sum(abs(CO2_ExperimentalResults.df$Angelo.Residual)),
                                 1/length(CO2_ExperimentalResults.df$Ruckenstein)*sum(abs(CO2_ExperimentalResults.df$Ruckenstein.Residual)))

ModelAccuracy.df$medianResidual = c(median(abs(CO2_ExperimentalResults.df$Amokrane.Residual)),
                             median(abs(CO2_ExperimentalResults.df$Hsu.Residual)),
                             median(abs(CO2_ExperimentalResults.df$Angelo.Residual)),
                             median(abs(CO2_ExperimentalResults.df$Ruckenstein.Residual)))

ModelAccuracy.df$meanPercenterror = c(1/length(CO2_ExperimentalResults.df$Amokrane)*sum(abs(CO2_ExperimentalResults.df$Amokrane.Percenterror)),
                                 1/length(CO2_ExperimentalResults.df$Hsu)*sum(abs(CO2_ExperimentalResults.df$Hsu.Percenterror)),
                                 1/length(CO2_ExperimentalResults.df$Angelo)*sum(abs(CO2_ExperimentalResults.df$Angelo.Percenterror)),
                                 1/length(CO2_ExperimentalResults.df$Ruckenstein)*sum(abs(CO2_ExperimentalResults.df$Ruckenstein.Percenterror)))

ModelAccuracy.df$medianPercenterror = c(median(abs(CO2_ExperimentalResults.df$Amokrane.Percenterror)),
                                    median(abs(CO2_ExperimentalResults.df$Hsu.Percenterror)),
                                    median(abs(CO2_ExperimentalResults.df$Angelo.Percenterror)),
                                    median(abs(CO2_ExperimentalResults.df$Ruckenstein.Percenterror)))


#Creating table to highlight z score outliers 
Outliers.table = data.frame(Z.Score.Range = c('<=1','<=3','<=10','>10'))

Outliers.table$Amokrane = c(sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)<=1),
                                    sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)>1&abs(CO2_ExperimentalResults.df$Amokrane.ZScore)<=3),
                                    sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)>3&abs(CO2_ExperimentalResults.df$Amokrane.ZScore)<=10),
                            sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)>10))

Outliers.table$Hsu = c(sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)<=1),
                                    sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)>1&abs(CO2_ExperimentalResults.df$Hsu.ZScore)<=3),
                                    sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)>3&abs(CO2_ExperimentalResults.df$Hsu.ZScore)<=10),
                                    sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)>10))

Outliers.table$Angelo = c(sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)<=1),
                                  sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)>1&abs(CO2_ExperimentalResults.df$Angelo.ZScore)<=3),
                                  sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)>3&abs(CO2_ExperimentalResults.df$Angelo.ZScore)<=10),
                                    sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)>10))

Outliers.table$Ruckenstein = c(sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)<=1),
                                       sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)>1&abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)<=3),
                                       sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)>3&abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)<=10),
                                    sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)>10))

Outliers.long <- Outliers.table %>%
  pivot_longer(cols = -Z.Score.Range, 
               names_to = "Model", 
               values_to = "Count") %>%
  mutate(Z.Score.Range = factor(Z.Score.Range, levels = c('>10','<=10','<=3','<=1')))

Outliers.long <- Outliers.long %>%
  group_by(Model) %>%
  mutate(Percent = Count / sum(Count) * 100)

#creating a stacked bar plot of z score outliers 
zScore.breakdown.plot = ggplot(Outliers.long, aes(x=Model,y=Percent,fill=Z.Score.Range))+
  theme_bw()+
  geom_bar(stat='identity')+
  scale_fill_manual(values = c(safe_colorblind_palette[7],safe_colorblind_palette[10],safe_colorblind_palette[12],safe_colorblind_palette[5]),
                    guide = guide_legend(reverse=T)) +
  ylab("Percent of Predictions \nWithin |Z-Score| Range") +
  xlab("Model") +
  labs(fill = "|Z-Score|") +
  theme(text = element_text(family = "sans"))+
  theme(aspect.ratio=.5)+
  theme(text=element_text(size=20))+
  theme(legend.text=element_text(size=20))+
  theme(axis.text = element_text(size = 20))+
  theme(legend.position = 'bottom')
  
print(zScore.breakdown.plot)  

#Creating a data frame for boxplot
accuracy.df = data.frame(zScore = c(CO2_ExperimentalResults.df$Amokrane.ZScore,CO2_ExperimentalResults.df$Hsu.ZScore,
                                        CO2_ExperimentalResults.df$Angelo.ZScore,CO2_ExperimentalResults.df$Ruckenstein.ZScore),
                         residual = c(CO2_ExperimentalResults.df$Amokrane.Percenterror,CO2_ExperimentalResults.df$Hsu.Percenterror,
                                          CO2_ExperimentalResults.df$Angelo.Percenterror,CO2_ExperimentalResults.df$Ruckenstein.Percenterror),
                         SD = rep(CO2_ExperimentalResults.df$CO2.SD.gpL, 4),
                             Model = factor(rep(c('Amokrane','Hsu','Angelo','Ruckenstein'),
                                                each = length(CO2_ExperimentalResults.df$Amokrane.ZScore))))
accuracy.df$Model <- factor(accuracy.df$Model,
                                levels = c('Amokrane','Hsu','Angelo','Ruckenstein'))

Model.Boxplot = ggplot(data=accuracy.df,aes(x=Model,y=zScore,fill=Model))+
  theme_bw()+
  geom_boxplot(outlier.shape=NA)+
  geom_hline(yintercept=c(-3,3),linetype='dashed',color='black',linewidth=0.8)+
  geom_jitter(width = 0.07, size = 1.5,alpha=0.5)+
  scale_fill_manual(values = c(safe_colorblind_palette[3],safe_colorblind_palette[4],safe_colorblind_palette[2],safe_colorblind_palette[11]))+
  theme(legend.position='none')+
  xlab(bquote(paste('\n ',k[L],' Model')))+
  ylab('Z-Score of Model Predictions')+
  theme(text = element_text(family = "sans"))+
  theme(legend.text=element_text(size=20))+
  theme(axis.text = element_text(size = 20))+
  theme(aspect.ratio=.5)+
  theme(text=element_text(size=20))

print(Model.Boxplot)

Model.Scatterplot = ggplot(data=accuracy.df,aes(x=abs(zScore),y=abs(residual),color=Model,shape=Model,size=SD))+
  theme_bw()+
  geom_point(alpha=0.7)+
  ylab('Percent Error between\nModel and Experimental Values')+
  xlab('|Z-Score| of Model Predictions')+
  scale_color_manual(values = c(safe_colorblind_palette[3],safe_colorblind_palette[4],safe_colorblind_palette[2],safe_colorblind_palette[11]))+
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  scale_size_continuous(range = c(2,8)) + 
  guides(color = guide_legend(nrow=1,title=NULL,override.aes = list(size=4)),
         shape = guide_legend(nrow=1, title=NULL,override.aes=list(size=4)),
         size = guide_legend(nrow=1, title='Standard Deviation'))+
  theme(legend.position = 'bottom', legend.box='vertical',legend.box.just='center')+
  scale_x_log10()+
  geom_vline(xintercept=3,linetype='dashed',color='black',linewidth=0.8)+
  annotate("text", y = max(abs(accuracy.df$residual))*0.8, x = 1, label = "|Z-score| = 3", vjust = -0.5,
           size = 8, family = 'sans')+
  theme(text = element_text(family = "sans"))+
  theme(legend.text=element_text(size=20))+
  theme(axis.text = element_text(size = 20))+
  theme(aspect.ratio=.5)+
  theme(text=element_text(size=20))
  
print(Model.Scatterplot)


#Importing Liquid Holdup 
LiquidHoldup.RawData = read.csv('Experimental Data/LH_Experimental.csv')
n.holdup = 10                #holdup replicates per condition
u.Flowrate.LPM = 0.5/sqrt (3)  #flowmeter graduations of 1 L/min

LiquidHoldup.Values = data.frame(FlowRate.LPM = rep(c(3,4.5),each=3),
                                Height.m = rep(c(0.53,0.68,0.98),2))


LiquidHoldup.Values$VL.L = as.vector(tapply(LiquidHoldup.RawData$LiquidHoldupVolume.L,
                             list(LiquidHoldup.RawData$Height.m,LiquidHoldup.RawData$Flowrate.LPM),FUN=mean))

LiquidHoldup.Values$VL.L.SD = as.vector(tapply(LiquidHoldup.RawData$LiquidHoldupVolume.L,
                                            list(LiquidHoldup.RawData$Height.m,LiquidHoldup.RawData$Flowrate.LPM),FUN=sd))
 

#Experimental residence time and its uncertainty (standard error of mean holdup + flowmeter)
LiquidHoldup.Values$tau.s = LiquidHoldup.Values$VL.L/(LiquidHoldup.Values$FlowRate.LPM/60)
LiquidHoldup.Values$tau.u = LiquidHoldup.Values$tau.s*sqrt((LiquidHoldup.Values$VL.L.SD/sqrt(n.holdup)/LiquidHoldup.Values$VL.L)^2+
                                                             (u.Flowrate.LPM/LiquidHoldup.Values$FlowRate.LPM)^2)

LiquidHoldup.3LPM.Model = read.csv('Model Results/3LPM_LH_t_model_10272026.csv')
LiquidHoldup.4.5LPM.Model =read.csv('Model Results/4.5LPM_LH_t_model_10272026.csv')

#LH Plotting 
LH.Plot.df = data.frame(Flowrate = rep(c('3.0 L/min Experimental','4.5 L/min Experimental','3.0 L/min Model','4.5 L/min Model' ),each=3),
                        Height = rep(c(0.53,0.68,0.98),4),
                        VL.L = c(LiquidHoldup.Values$VL.L,LiquidHoldup.3LPM.Model$volume.holdup.L,LiquidHoldup.4.5LPM.Model$volume.holdup.L),
                        VL.L.SD=c(LiquidHoldup.Values$VL.L.SD,0,0,0,0,0,0),
                        Type = rep(c('Experimental','Model'),each=6))

#Generating powerlaw fits for height vs. holdup volume 
fit_powerlaw = function(df){
  fit = nls(VL.L~a*Height^b,data=df,start=list(a=0.05,b=1))
  a = coef(fit)[1]
  b = coef(fit)[2]
  
  r2 = 1-sum(residuals(fit)^2)/sum((df$VL.L-mean(df$VL.L))^2)
  
  df_sorted <- df[order(df$Height), ]
  x_pos <- mean(df_sorted$Height)
  y_pos <- tail(df_sorted$VL.L, 1)
  
  data.frame(Flowrate = df$Flowrate[1],label = c(
    paste0("y == ", format(a, digits=3), " * x^", format(b, digits=3)), 
    paste0("R^2 == ", format(r2, digits=3))),x=c(x_pos,x_pos),y=c(y_pos,y_pos-0.004))
}

label.df <- LH.Plot.df %>%
  group_by(Flowrate) %>%
  do(fit_powerlaw(.)) 

#Plotting liquid holdup and models 
LiquidHoldup_plot = ggplot(LH.Plot.df, aes(x=Height,y=VL.L,group=Flowrate,color=Flowrate,shape=Flowrate,fill=Flowrate))+
  theme_bw()+
  geom_point(size=4)+
  geom_errorbar(data=subset(LH.Plot.df,Type=='Experimental'),aes(ymin=VL.L-VL.L.SD,ymax=VL.L+VL.L.SD,color=Flowrate),width=0.01,linewidth=0.9)+
  geom_smooth(data=subset(LH.Plot.df,Type=='Experimental'),method='nls',formula = y~a*x^b,
              method.args=list(start=list(a=0.05,b=1)),se=F,aes(color=Flowrate),linetype='dashed',linewidth=1)+
  geom_smooth(data=subset(LH.Plot.df,Type=='Model'),method='nls',formula = y~a*x^b,
              method.args=list(start=list(a=0.05,b=1)),se=F,aes(color=Flowrate),linetype='dotted',linewidth=1)+
  theme(text = element_text(family = "sans"))+
  theme(aspect.ratio=.5)+
  theme(text=element_text(size=20))+
  theme(legend.position='bottom')+
  theme(legend.title = element_blank())+
  ylab('Liquid Holdup Volume (L)\n')+
  xlab('\nHeight from Droplet Manifold (m)')+
  theme(legend.text=element_text(size=20))+
  theme(axis.text = element_text(size = 20))+
  scale_shape_manual(values=c('3.0 L/min Experimental'= 16,'4.5 L/min Experimental'=17,'3.0 L/min Model'=15,'4.5 L/min Model'=23))+
  scale_color_manual(values=c('3.0 L/min Experimental'= '#88CCEE','4.5 L/min Experimental'="#882255",'3.0 L/min Model'="#999933",'4.5 L/min Model'="#888888"))+
  scale_fill_manual(values=c('3.0 L/min Experimental'= '#88CCEE','4.5 L/min Experimental'="#882255",'3.0 L/min Model'="#999933",'4.5 L/min Model'="#888888"))+
  geom_text(data=label.df, aes(x=rep(c(0.9,0.9,0.9,0.9),each=2),y=c(0.057,0.052,0.013,0.007,0.076,0.071,0.036,0.031),label=label),color='black',parse=T,hjust =0, vjust=0, size=5,family='sans')+
  ylim(0,0.1)

print(LiquidHoldup_plot)

#Computations for kLa 

#Instrument accuracy as standard uncertainties (+/- limit divided by sqrt(3))
u.CO2.gpL = 0.05/sqrt(3)     #Carbo 5100 accuracy +/- 0.05 g/L
u.Temp.C = 0.2/sqrt(3)       #thermometer accuracy estimate; same value as the model
u.Pressure.atm = 0.2/14.7    #day-to-day barometric variation, 0.2 psi; same value as the model
u.Flowrate.LPM = 0.5/sqrt (3)  #flowmeter graduations of 1 L/min
n.holdup = 10                #holdup replicates per condition 

#Saturation concentration (g/L)
Csat.fun = function(Temp.K,Pressure.atm){
  H.atmm3.mol = (1/(exp(-159.854+8741.68/Temp.K+21.6694*log(Temp.K)-0.0011026*Temp.K)))/55342
  H.atmL.g = H.atmm3.mol*1000/44
  return(Pressure.atm/H.atmL.g)
}

#Experimental Computations 
kLa_Experimental.df = CleanData.df 
for (i in 1:length(kLa_Experimental.df)){
  kLa_Experimental.df[[i]]$Temp.K = kLa_Experimental.df[[i]]$`Liquid Temperature (C)`+273.15
  kLa_Experimental.df[[i]]$Pressure.atm = (kLa_Experimental.df[[i]]$`CO2 Pressure (psig)`+14.2)/14.7
  kLa_Experimental.df[[i]]$Flowrate.Lps = kLa_Experimental.df[[i]]$`Flowrate (L/min)`/60
  kLa_Experimental.df[[i]]$Csat = Csat.fun(kLa_Experimental.df[[i]]$Temp.K,kLa_Experimental.df[[i]]$Pressure.atm)
  
  FlowRate = kLa_Experimental.df[[i]]$`Flowrate (L/min)`[1]
  Height = kLa_Experimental.df[[i]]$Height.m[1]
  kLa_Experimental.df[[i]]$VL.L = rep(LiquidHoldup.Values$VL.L[LiquidHoldup.Values$FlowRate.LPM==FlowRate & LiquidHoldup.Values$Height.m==Height],length(kLa_Experimental.df[[i]]$Trial))
  kLa_Experimental.df[[i]]$VL.L.SD = rep(LiquidHoldup.Values$VL.L.SD[LiquidHoldup.Values$FlowRate.LPM==FlowRate & LiquidHoldup.Values$Height.m==Height],length(kLa_Experimental.df[[i]]$Trial))
  kLa_Experimental.df[[i]]$VL.L.u = kLa_Experimental.df[[i]]$VL.L.SD/sqrt(n.holdup)
  
  kLa_Experimental.df[[i]]$kLa.s = kLa_Experimental.df[[i]]$Flowrate.Lps/kLa_Experimental.df[[i]]$VL.L*
    log((kLa_Experimental.df[[i]]$Csat)/(kLa_Experimental.df[[i]]$Csat-kLa_Experimental.df[[i]]$`CO2 (g/L)`))
  
  #How much kLa changes per unit change in each input
  kLa_Experimental.df[[i]]$dkLa.dC = kLa_Experimental.df[[i]]$Flowrate.Lps/kLa_Experimental.df[[i]]$VL.L/
    (kLa_Experimental.df[[i]]$Csat-kLa_Experimental.df[[i]]$`CO2 (g/L)`)
  kLa_Experimental.df[[i]]$dkLa.dCsat = -kLa_Experimental.df[[i]]$Flowrate.Lps/kLa_Experimental.df[[i]]$VL.L*
    kLa_Experimental.df[[i]]$`CO2 (g/L)`/(kLa_Experimental.df[[i]]$Csat*(kLa_Experimental.df[[i]]$Csat-kLa_Experimental.df[[i]]$`CO2 (g/L)`))
  kLa_Experimental.df[[i]]$dCsat.dT = (Csat.fun(kLa_Experimental.df[[i]]$Temp.K+0.01,kLa_Experimental.df[[i]]$Pressure.atm)-
                                         Csat.fun(kLa_Experimental.df[[i]]$Temp.K-0.01,kLa_Experimental.df[[i]]$Pressure.atm))/0.02
  kLa_Experimental.df[[i]]$dCsat.dP = kLa_Experimental.df[[i]]$Csat/kLa_Experimental.df[[i]]$Pressure.atm
  
  #Systematic uncertainty: the same for every reading, so averaging does NOT reduce it
  kLa_Experimental.df[[i]]$kLa.s.u.sys = sqrt((kLa_Experimental.df[[i]]$kLa.s*kLa_Experimental.df[[i]]$VL.L.u/kLa_Experimental.df[[i]]$VL.L)^2+
                                                (kLa_Experimental.df[[i]]$dkLa.dC*u.CO2.gpL)^2+
                                                (kLa_Experimental.df[[i]]$dkLa.dCsat*kLa_Experimental.df[[i]]$dCsat.dT*u.Temp.C)^2+
                                                (kLa_Experimental.df[[i]]$dkLa.dCsat*kLa_Experimental.df[[i]]$dCsat.dP*u.Pressure.atm)^2+
                                                (kLa_Experimental.df[[i]]$kLa.s*u.Flowrate.LPM/kLa_Experimental.df[[i]]$`Flowrate (L/min)`)^2)
}

#Creating data frame of all trials
kLa_Results.df = data.frame(FlowRate.LPM = rep(c(3,4.5),each=9),
                            Height.m = rep(c(rep(c(0.53,0.68,0.98),each=3)),2),
                            Trial = rep(c(1,2,3),6))

meankLa.s = vector()
sdkLa.s = vector()
syskLa.s = vector()
meanP.atm = vector()

#Determine mean and error in each trial at each height 
for (i in 1:length(kLa_Experimental.df)){
  mean.new.kLa = tapply(kLa_Experimental.df[[i]]$kLa.s,kLa_Experimental.df[[i]]$Trial,FUN=mean)
  sys.new.kLa = tapply(kLa_Experimental.df[[i]]$kLa.s.u.sys,kLa_Experimental.df[[i]]$Trial,FUN=mean)
  rand.new.kLa = tapply(kLa_Experimental.df[[i]]$kLa.s,kLa_Experimental.df[[i]]$Trial,FUN=sd)
  error.new.kLa = sqrt(sys.new.kLa^2+rand.new.kLa^2)
  mean.new.P = tapply(kLa_Experimental.df[[i]]$Pressure.atm,kLa_Experimental.df[[i]]$Trial,FUN=mean)
  
  meankLa.s = c(meankLa.s,mean.new.kLa)
  sdkLa.s = c(sdkLa.s,error.new.kLa)
  syskLa.s = c(syskLa.s,sys.new.kLa)
  meanP.atm = c(meanP.atm,mean.new.P)
}

kLa_Results.df$kLa = meankLa.s
kLa_Results.df$kLa.sd = sdkLa.s
kLa_Results.df$kLa.sys = syskLa.s
kLa_Results.df$Pressure.atm = meanP.atm


#kLa Modeling Results 
kLa_Model.df = data.frame(Temp.K=CO2_ExperimentalResults.df$Temp.C+273.15,
                          CO2.Angelo = CO2_ExperimentalResults.df$Angelo,
                          vL.L = c(rep(LiquidHoldup.3LPM.Model$volume.holdup.L,each=3),rep(LiquidHoldup.4.5LPM.Model$volume.holdup.L,each=3)),
                          flowrate.lps = rep(c(3.0/60,4.5/60),each=9),
                          height.m= rep(c(rep(c(0.53,0.68,0.98),each=3)),2))

kLa_Model.df$CSat = Csat.fun(kLa_Model.df$Temp.K,kLa_Results.df$Pressure.atm)   #measured trial pressure, not fixed 5 psig

kLa_Model.df$kLa.s = kLa_Model.df$flowrate.lps/kLa_Model.df$vL.L*log(kLa_Model.df$CSat/(kLa_Model.df$CSat-kLa_Model.df$CO2.Angelo))

kLa_Results.df$kLa.model = kLa_Model.df$kLa.s
kLa_Results.df$kLa.model.u = ModelU$Angelo.kLa.u   #same row order (flow, height, trial)
kLa_Results.df$Temp.C = kLa_Model.df$Temp.K-273.15

#Composite kLa results 
kLa.overall.df = data.frame(Flowrate = rep(c(3.0,4.5),each=3),
                            Height = rep(c(0.53,0.68,0.98),2))

kLa.overall.df$Temp = as.vector(tapply(kLa_Results.df$Temp.C,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))
kLa.overall.df$Temp.sd = as.vector(tapply(kLa_Results.df$Temp.C,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=sd))
kLa.overall.df$model.kLa = as.vector(tapply(kLa_Results.df$kLa.model,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))
kLa.overall.df$model.kLa.sd = as.vector(tapply(kLa_Results.df$kLa.model,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=sd))
#model kLa uncertainty from the inputs (shared between trials, so averaged, not divided by sqrt(3))
kLa.overall.df$model.kLa.u = as.vector(tapply(kLa_Results.df$kLa.model.u,
                                              list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))

kLa.overall.df$Experimental.kLa = as.vector(tapply(kLa_Results.df$kLa,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))
kLa.random = as.vector(tapply(kLa_Results.df$kLa,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=sd))/sqrt(3)   #wobble between trials, shrinks with averaging
kLa.sys = as.vector(tapply(kLa_Results.df$kLa.sys,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))       #shared error, does not shrink
kLa.overall.df$Experimental.kLa.sd = sqrt(kLa.random^2+kLa.sys^2)

Table1.df = CO2_ExperimentalResults.df[,!(names(CO2_ExperimentalResults.df) %in% c('median','medianDiff','mad','ZScore','outlier'))]
kLa.overall.df$model.kLa.sd = NULL   #spread across trial temperatures, not an uncertainty; model.kLa.u is the uncertainty

#### kLa and absorption rate per unit column volume (R2-1) -- replaces DataAnalysis_ColumnBasis_28Oct.R ####
#Paste after the model kLa block and BEFORE writexl::write_xlsx(), so the new columns export with Table4_kLa.
#Column = cylinder enclosing the manifold (manifold diameter + 0.1 m); height = fall height of each condition.
#V_C is a defined reference volume, so it is not assigned an uncertainty.
D.col.m = 9*0.0254+0.1                 #0.329 m
A.col.m2 = pi/4*D.col.m^2

## Experimental: same propagation as kLa, with the holdup term removed (holdup cancels on this basis)
for (i in 1:length(kLa_Experimental.df)){
  d = kLa_Experimental.df[[i]]
  V.col.L = A.col.m2*d$Height.m*1000
  d$kLa.col.s = d$Flowrate.Lps/V.col.L*log(d$Csat/(d$Csat-d$`CO2 (g/L)`))
  s = d$VL.L/V.col.L                   #kLa-block derivatives are per V_L; rescale to V_C
  d$kLa.col.u.sys = sqrt((s*d$dkLa.dC*u.CO2.gpL)^2+
                           (s*d$dkLa.dCsat*d$dCsat.dT*u.Temp.C)^2+
                           (s*d$dkLa.dCsat*d$dCsat.dP*u.Pressure.atm)^2+
                           (d$kLa.col.s*u.Flowrate.LPM/d$`Flowrate (L/min)`)^2)
  d$Rate.gm3s = d$Flowrate.Lps*d$`CO2 (g/L)`/V.col.L*1000   #g CO2 per s per m3 of column (inlet C = 0)
  kLa_Experimental.df[[i]] = d
}

col.mean = c(); col.sys = c(); rate.mean = c()
for (i in 1:length(kLa_Experimental.df)){
  d = kLa_Experimental.df[[i]]
  col.mean  = c(col.mean,  tapply(d$kLa.col.s,d$Trial,FUN=mean))
  col.sys   = c(col.sys,   tapply(d$kLa.col.u.sys,d$Trial,FUN=mean))
  rate.mean = c(rate.mean, tapply(d$Rate.gm3s,d$Trial,FUN=mean))
}
kLa_Results.df$kLa.col = col.mean
kLa_Results.df$kLa.col.sys = col.sys
kLa_Results.df$Rate.gm3s = rate.mean

grp = list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM)
kLa.overall.df$Experimental.kLa.col = as.vector(tapply(kLa_Results.df$kLa.col,grp,FUN=mean))
col.random = as.vector(tapply(kLa_Results.df$kLa.col,grp,FUN=sd))/sqrt(3)      #between trials, shrinks with averaging
col.sys    = as.vector(tapply(kLa_Results.df$kLa.col.sys,grp,FUN=mean))         #shared, does not shrink
kLa.overall.df$Experimental.kLa.col.u = sqrt(col.random^2+col.sys^2)
kLa.overall.df$Experimental.Rate.gm3s = as.vector(tapply(kLa_Results.df$Rate.gm3s,grp,FUN=mean))

## Model: one-at-a-time propagation using the signed changes in C from runUncertainty.m (ModelU)
#CHECK names(ModelU) and rename the five Angelo.d* columns below if they differ
stopifnot(all(ModelU$FlowRate.LPM==CO2_ExperimentalResults.df$FlowRate.LPM),
          all(ModelU$Height.m==CO2_ExperimentalResults.df$Height.m),
          all(ModelU$Trial==CO2_ExperimentalResults.df$Trial))
kLa_Model.df$V.col.L = A.col.m2*kLa_Model.df$height.m*1000
C  = kLa_Model.df$CO2.Angelo
Cs = kLa_Model.df$CSat
P  = kLa_Results.df$Pressure.atm
QV = kLa_Model.df$flowrate.lps/kLa_Model.df$V.col.L
kLa_Model.df$kLa.col.s = QV*log(Cs/(Cs-C))
kLa_Model.df$Rate.gm3s = QV*C*1000
dfdC  = QV/(Cs-C)
dfdCs = -QV*C/(Cs*(Cs-C))
dCs.dT = (Csat.fun(kLa_Model.df$Temp.K+0.01,P)-Csat.fun(kLa_Model.df$Temp.K-0.01,P))/0.02
d.T   = dfdC*ModelU$Angelo.dT + dfdCs*dCs.dT*u.Temp.C
d.P   = dfdC*ModelU$Angelo.dP + dfdCs*Cs/P*u.Pressure.atm
d.Q   = dfdC*ModelU$Angelo.dQ + kLa_Model.df$kLa.col.s*u.Flowrate.LPM/(kLa_Model.df$flowrate.lps*60)
d.DAB = dfdC*ModelU$Angelo.dDAB
d.do  = dfdC*ModelU$Angelo.dd_o
#Henry's constant scales C and C* by the same fraction, so ln(C*/(C*-C)) is unchanged: no term
kLa_Model.df$kLa.col.u = sqrt(d.T^2+d.P^2+d.Q^2+d.DAB^2+d.do^2)

grp.m = list(kLa_Model.df$height.m,kLa_Model.df$flowrate.lps)
kLa.overall.df$model.kLa.col   = as.vector(tapply(kLa_Model.df$kLa.col.s,grp.m,FUN=mean))
kLa.overall.df$model.kLa.col.u = as.vector(tapply(kLa_Model.df$kLa.col.u,grp.m,FUN=mean))   #shared inputs: not /sqrt(3)
kLa.overall.df$model.Rate.gm3s = as.vector(tapply(kLa_Model.df$Rate.gm3s,grp.m,FUN=mean))
kLa.overall.df$V.col.L = A.col.m2*kLa.overall.df$Height*1000

print(kLa.overall.df[,c('Flowrate','Height','V.col.L','Experimental.kLa.col','Experimental.kLa.col.u',
                        'model.kLa.col','model.kLa.col.u','Experimental.Rate.gm3s','model.Rate.gm3s')],digits=3)
#Expected (6 Nov estimate from Tables 3-4): exp kLa_col x1e-4 = 8.0+/-0.9, 6.7+/-0.7, 5.1+/-0.7, 13.5+/-1.5, 11.5+/-1.3, 9.1+/-1.0;
#model ~7.1, 5.9, 5.3, 10.0, 8.7, 7.8 (x1e-4) with u of roughly 3%.

Table5_ColumnBasis = data.frame(
  Flowrate.LPM        = kLa.overall.df$Flowrate,
  Height.m            = kLa.overall.df$Height,
  V.col.L             = kLa.overall.df$V.col.L,
  Exp.kLa.col.e4      = kLa.overall.df$Experimental.kLa.col*1e4,
  Exp.kLa.col.u.e4    = kLa.overall.df$Experimental.kLa.col.u*1e4,
  Model.kLa.col.e4    = kLa.overall.df$model.kLa.col*1e4,
  Model.kLa.col.u.e4  = kLa.overall.df$model.kLa.col.u*1e4,
  Exp.Rate.gm3s       = kLa.overall.df$Experimental.Rate.gm3s,
  Model.Rate.gm3s     = kLa.overall.df$model.Rate.gm3s)

writexl::write_xlsx(list(Table1_CO2 = Table1.df,
                         Table2_ModelAccuracy = ModelAccuracy.df,
                         Table3_Holdup = LiquidHoldup.Values,
                         Holdup_Model_3.0LPM = LiquidHoldup.3LPM.Model,
                         Holdup_Model_4.5LPM = LiquidHoldup.4.5LPM.Model,
                         Table4_kLa = kLa.overall.df,
                         Table5_ColumnBasis = Table5_ColumnBasis),
                    'ManuscriptTables.xlsx')