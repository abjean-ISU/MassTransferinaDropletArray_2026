#Liquid Holdup/Experimental results data processing
rm(list=ls())

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


#Importing Mass Transfer Experimental Data 
ExperimentalData.Filepath = 'MassTransfer Experimental Data 10222025.xlsx'
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

CO2_ExperimentalResults.df$CO2.gpL = meanCO2.gpL
CO2_ExperimentalResults.df$CO2.SD.gpL = sdCO2.gpL+0.005
CO2_ExperimentalResults.df$Temp.C = meanT.C
CO2_ExperimentalResults.df$Temp.SD.C = sdT.C+0.005 

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
High.3LPM.T1.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 High 10222025\\T1\\MT_Models.csv')
High.3LPM.T2.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 High 10222025\\T2\\MT_Models.csv')
High.3LPM.T3.ModelResults = read.csv('3LPM CO2 High 10222025\\T3\\MT_Models.csv')

Med.3LPM.T1.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 Med 10222025\\T1\\MT_Models.csv')
Med.3LPM.T2.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 Med 10222025\\T2\\MT_Models.csv')
Med.3LPM.T3.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 Med 10222025\\T3\\MT_Models.csv')

Low.3LPM.T1.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 Low 10222025\\T1\\MT_Models.csv')
Low.3LPM.T2.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 Low 10222025\\T2\\MT_Models.csv')
Low.3LPM.T3.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\3LPM CO2 Low 10222025\\T3\\MT_Models.csv')


High.4.5LPM.T1.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 High 10222025\\T1\\MT_Models.csv')
High.4.5LPM.T2.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 High 10222025\\T2\\MT_Models.csv')
High.4.5LPM.T3.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 High 10222025\\T3\\MT_Models.csv')

Med.4.5LPM.T1.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 Med 10222025\\T1\\MT_Models.csv')
Med.4.5LPM.T2.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 Med 10222025\\T2\\MT_Models.csv')
Med.4.5LPM.T3.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 Med 10222025\\T3\\MT_Models.csv')

Low.4.5LPM.T1.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 Low 10222025\\T1\\MT_Models.csv')
Low.4.5LPM.T2.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 Low 10222025\\T2\\MT_Models.csv')
Low.4.5LPM.T3.ModelResults = read.csv('Mass Transfer & Liquid Holdup Models\\4.5LPM CO2 Low 10222025\\T3\\MT_Models.csv')

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
    MaxT3.0_dataset = Low.4.5LPM.T2.ModelResults
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
  geom_errorbar(data=ExperimentalResults_3.0LPM,aes(x=Height.m,ymin=(CO2.gpL - CO2.SD.gpL),
                                                      ymax = (CO2.gpL + CO2.SD.gpL)),
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
  geom_errorbar(data=ExperimentalResults_4.5LPM,aes(x=Height.m,ymin=(CO2.gpL - CO2.SD.gpL),
                                                    ymax = (CO2.gpL + CO2.SD.gpL)), width=0)+
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

#Determining model results for each experiment 
CO2_ExperimentalResults.df$Amokrane = c(Low.3LPM.T1.ModelResults$Amokrane[Low.3LPM.T1.ModelResults$Height.m==0.53],
                        Low.3LPM.T2.ModelResults$Amokrane[Low.3LPM.T2.ModelResults$Height.m==0.53],
                        Low.3LPM.T3.ModelResults$Amokrane[Low.3LPM.T3.ModelResults$Height.m==0.53],
                        Med.3LPM.T1.ModelResults$Amokrane[Med.3LPM.T1.ModelResults$Height.m==0.68],
                        Med.3LPM.T2.ModelResults$Amokrane[Med.3LPM.T2.ModelResults$Height.m==0.68],
                        Med.3LPM.T3.ModelResults$Amokrane[Med.3LPM.T3.ModelResults$Height.m==0.68],
                        High.3LPM.T1.ModelResults$Amokrane[High.3LPM.T1.ModelResults$Height.m==0.98],
                        High.3LPM.T2.ModelResults$Amokrane[High.3LPM.T2.ModelResults$Height.m==0.98],
                        High.3LPM.T3.ModelResults$Amokrane[High.3LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Amokrane.Residual = (CO2_ExperimentalResults.df$Amokrane-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Amokrane.Percenterror = abs(CO2_ExperimentalResults.df$Amokrane.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Amokrane.ZScore = (CO2_ExperimentalResults.df$Amokrane-CO2_ExperimentalResults.df$CO2.gpL)/CO2_ExperimentalResults.df$CO2.SD.gpL

CO2_ExperimentalResults.df$Hsu = c(Low.3LPM.T1.ModelResults$Hsu[Low.3LPM.T1.ModelResults$Height.m==0.53],
                        Low.3LPM.T2.ModelResults$Hsu[Low.3LPM.T2.ModelResults$Height.m==0.53],
                        Low.3LPM.T3.ModelResults$Hsu[Low.3LPM.T3.ModelResults$Height.m==0.53],
                        Med.3LPM.T1.ModelResults$Hsu[Med.3LPM.T1.ModelResults$Height.m==0.68],
                        Med.3LPM.T2.ModelResults$Hsu[Med.3LPM.T2.ModelResults$Height.m==0.68],
                        Med.3LPM.T3.ModelResults$Hsu[Med.3LPM.T3.ModelResults$Height.m==0.68],
                        High.3LPM.T1.ModelResults$Hsu[High.3LPM.T1.ModelResults$Height.m==0.98],
                        High.3LPM.T2.ModelResults$Hsu[High.3LPM.T2.ModelResults$Height.m==0.98],
                        High.3LPM.T3.ModelResults$Hsu[High.3LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Hsu.Residual = (CO2_ExperimentalResults.df$Hsu-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Hsu.Percenterror = abs(CO2_ExperimentalResults.df$Hsu.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Hsu.ZScore = (CO2_ExperimentalResults.df$Hsu-CO2_ExperimentalResults.df$CO2.gpL)/CO2_ExperimentalResults.df$CO2.SD.gpL

CO2_ExperimentalResults.df$Angelo = c(Low.3LPM.T1.ModelResults$Angelo[Low.3LPM.T1.ModelResults$Height.m==0.53],
                   Low.3LPM.T2.ModelResults$Angelo[Low.3LPM.T2.ModelResults$Height.m==0.53],
                   Low.3LPM.T3.ModelResults$Angelo[Low.3LPM.T3.ModelResults$Height.m==0.53],
                   Med.3LPM.T1.ModelResults$Angelo[Med.3LPM.T1.ModelResults$Height.m==0.68],
                   Med.3LPM.T2.ModelResults$Angelo[Med.3LPM.T2.ModelResults$Height.m==0.68],
                   Med.3LPM.T3.ModelResults$Angelo[Med.3LPM.T3.ModelResults$Height.m==0.68],
                   High.3LPM.T1.ModelResults$Angelo[High.3LPM.T1.ModelResults$Height.m==0.98],
                   High.3LPM.T2.ModelResults$Angelo[High.3LPM.T2.ModelResults$Height.m==0.98],
                   High.3LPM.T3.ModelResults$Angelo[High.3LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Angelo.Residual = (CO2_ExperimentalResults.df$Angelo-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Angelo.Percenterror = abs(CO2_ExperimentalResults.df$Angelo.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Angelo.ZScore = (CO2_ExperimentalResults.df$Angelo-CO2_ExperimentalResults.df$CO2.gpL)/CO2_ExperimentalResults.df$CO2.SD.gpL

CO2_ExperimentalResults.df$Ruckenstein = c(Low.3LPM.T1.ModelResults$Ruckenstein[Low.3LPM.T1.ModelResults$Height.m==0.53],
                      Low.3LPM.T2.ModelResults$Ruckenstein[Low.3LPM.T2.ModelResults$Height.m==0.53],
                      Low.3LPM.T3.ModelResults$Ruckenstein[Low.3LPM.T3.ModelResults$Height.m==0.53],
                      Med.3LPM.T1.ModelResults$Ruckenstein[Med.3LPM.T1.ModelResults$Height.m==0.68],
                      Med.3LPM.T2.ModelResults$Ruckenstein[Med.3LPM.T2.ModelResults$Height.m==0.68],
                      Med.3LPM.T3.ModelResults$Ruckenstein[Med.3LPM.T3.ModelResults$Height.m==0.68],
                      High.3LPM.T1.ModelResults$Ruckenstein[High.3LPM.T1.ModelResults$Height.m==0.98],
                      High.3LPM.T2.ModelResults$Ruckenstein[High.3LPM.T2.ModelResults$Height.m==0.98],
                      High.3LPM.T3.ModelResults$Ruckenstein[High.3LPM.T3.ModelResults$Height.m==0.98])

CO2_ExperimentalResults.df$Ruckenstein.Residual = (CO2_ExperimentalResults.df$Ruckenstein-CO2_ExperimentalResults.df$CO2.gpL)
CO2_ExperimentalResults.df$Ruckenstein.Percenterror = abs(CO2_ExperimentalResults.df$Ruckenstein.Residual)/(CO2_ExperimentalResults.df$CO2.gpL)*100
CO2_ExperimentalResults.df$Ruckenstein.ZScore = (CO2_ExperimentalResults.df$Ruckenstein-CO2_ExperimentalResults.df$CO2.gpL)/CO2_ExperimentalResults.df$CO2.SD.gpL

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

Outliers.table$Amokrane = c(sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)<1),
                                    sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)>1&abs(CO2_ExperimentalResults.df$Amokrane.ZScore)<3),
                                    sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)>3&abs(CO2_ExperimentalResults.df$Amokrane.ZScore)<10),
                            sum(abs(CO2_ExperimentalResults.df$Amokrane.ZScore)>10))

Outliers.table$Hsu = c(sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)<1),
                                    sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)>1&abs(CO2_ExperimentalResults.df$Hsu.ZScore)<3),
                                    sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)>3&abs(CO2_ExperimentalResults.df$Hsu.ZScore)<10),
                                    sum(abs(CO2_ExperimentalResults.df$Hsu.ZScore)>10))

Outliers.table$Angelo = c(sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)<1),
                                  sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)>1&abs(CO2_ExperimentalResults.df$Angelo.ZScore)<3),
                                  sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)>3&abs(CO2_ExperimentalResults.df$Angelo.ZScore)<10),
                                    sum(abs(CO2_ExperimentalResults.df$Angelo.ZScore)>10))

Outliers.table$Ruckenstein = c(sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)<1),
                                       sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)>1&abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)<3),
                                       sum(abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)>3&abs(CO2_ExperimentalResults.df$Ruckenstein.ZScore)<10),
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
                    guide = guide_legend(revers=T)) +
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
  geom_boxplot(outlier.color = 'red',outlier.size =3)+
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
LiquidHoldup.RawData = read.csv('LH_Experimental.csv')

LiquidHoldup.Values = data.frame(FlowRate.LPM = rep(c(3,4.5),each=3),
                                Height.m = rep(c(0.53,0.68,0.98),2))


LiquidHoldup.Values$VL.L = as.vector(tapply(LiquidHoldup.RawData$LiquidHoldupVolume.L,
                             list(LiquidHoldup.RawData$Height.m,LiquidHoldup.RawData$Flowrate.LPM),FUN=mean))

LiquidHoldup.Values$VL.L.SD = as.vector(tapply(LiquidHoldup.RawData$LiquidHoldupVolume.L,
                                            list(LiquidHoldup.RawData$Height.m,LiquidHoldup.RawData$Flowrate.LPM),FUN=sd))

LiquidHoldup.3LPM.Model = read.csv('3LPM_LH_t_model.csv')
LiquidHoldup.4.5LPM.Model =read.csv('4.5LPM_LH_t_model.csv')

#LH Plotting 
LH.Plot.df = data.frame(Flowrate = rep(c('3.0 L/min Experimental','4.5 L/min Experimental','3.0 L/min Model','4.5 L/min Model' ),each=3),
                        Height = rep(c(0.53,0.68,0.98),4),
                        VL.L = c(LiquidHoldup.Values$VL.L,LiquidHoldup.3LPM.Model$volume.holdup.L,LiquidHoldup.4.5LPM.Model$volume.holdup.L),
                        VL.L.SD=c(LiquidHoldup.Values$VL.L.SD,0,0,0,0,0,0),
                        Type = rep(c('Experimental','Model'),each=6))

#Generating powerlaw fits for height vs. holdup volume 
fit_powerlaw = function(df){
  fit = lm(log(VL.L)~log(Height),data=df)
  a= exp(coef(fit)[1])
  b=coef(fit)[2]
  
  r2 = summary(fit)$r.squared 
  
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
  geom_smooth(data=subset(LH.Plot.df,Type=='Experimental'),method='lm',formula = y~log(x),se=F,aes(color=Flowrate),linetype='dashed',linewidth=1)+
  geom_smooth(data=subset(LH.Plot.df,Type=='Model'),method='lm',formula = y~log(x),se=F,aes(color=Flowrate),linetype='dotted',linewidth=1)+
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

#Experimental Computations 
kLa_Experimental.df = CleanData.df 
for (i in 1:length(kLa_Experimental.df)){
  kLa_Experimental.df[[i]]$Temp.K = kLa_Experimental.df[[i]]$`Liquid Temperature (C)`+273.15
  kLa_Experimental.df[[i]]$Pressure.atm = (kLa_Experimental.df[[i]]$`CO2 Pressure (psig)`+14.7)/14.696
  kLa_Experimental.df[[i]]$Flowrate.Lps = kLa_Experimental.df[[i]]$`Flowrate (L/min)`/60
  kLa_Experimental.df[[i]]$H.atmm3.mol = (1/(exp(-159.854+8741.68/kLa_Experimental.df[[i]]$Temp.K+
                                                21.6694*log(kLa_Experimental.df[[i]]$Temp.K)-0.0011026*kLa_Experimental.df[[i]]$Temp.K)))/55342
  kLa_Experimental.df[[i]]$H.atmL.g = kLa_Experimental.df[[i]]$H.atmm3.mol*1000/44
  kLa_Experimental.df[[i]]$Csat = kLa_Experimental.df[[i]]$Pressure.atm/kLa_Experimental.df[[i]]$H.atmL.g
  
  FlowRate = kLa_Experimental.df[[i]]$`Flowrate (L/min)`[1]
  Height = kLa_Experimental.df[[i]]$Height.m[1]
  kLa_Experimental.df[[i]]$VL.L = rep(LiquidHoldup.Values$VL.L[LiquidHoldup.Values$FlowRate.LPM==FlowRate & LiquidHoldup.Values$Height.m==Height],length(kLa_Experimental.df[[i]]$Trial))
  kLa_Experimental.df[[i]]$VL.L.SD = rep(LiquidHoldup.Values$VL.L.SD[LiquidHoldup.Values$FlowRate.LPM==FlowRate & LiquidHoldup.Values$Height.m==Height],length(kLa_Experimental.df[[i]]$Trial))
  
  
  kLa_Experimental.df[[i]]$kLa.s = kLa_Experimental.df[[i]]$Flowrate.Lps/kLa_Experimental.df[[i]]$VL.L*
    log((kLa_Experimental.df[[i]]$Csat)/(kLa_Experimental.df[[i]]$Csat-kLa_Experimental.df[[i]]$`CO2 (g/L)`))
  kLa_Experimental.df[[i]]$kLa.s.SD = sqrt((-kLa_Experimental.df[[i]]$Flowrate.Lps/kLa_Experimental.df[[i]]$VL.L^2*
                                              log((kLa_Experimental.df[[i]]$Csat)/(kLa_Experimental.df[[i]]$Csat-kLa_Experimental.df[[i]]$`CO2 (g/L)`))*
                                             kLa_Experimental.df[[i]]$VL.L.SD)^2)
  kLa_Experimental.df[[i]]$kLa.s.weight = 1/kLa_Experimental.df[[i]]$kLa.s.SD^2
}

#Creating data frame of all trials
kLa_Results.df = data.frame(FlowRate.LPM = rep(c(3,4.5),each=9),
                                        Height.m = rep(c(rep(c(0.53,0.68,0.98),each=3)),2),
                                        Trial = rep(c(1,2,3),6))

meankLa.s = vector()
sdkLa.s = vector()


#Determine mean and error in each trial at each height 
for (i in 1:length(kLa_Experimental.df)){
  mean.new.kLa = tapply(kLa_Experimental.df[[i]]$kLa.s,kLa_Experimental.df[[i]]$Trial,function(x,w){
    weighted.mean(x,w[match(x,kLa_Experimental.df[[i]]$kLa.s)])},w=kLa_Experimental.df[[i]]$kLa.s.weight)
  error.new.kLa = sqrt(1/tapply(kLa_Experimental.df[[i]]$kLa.s.weight,kLa_Experimental.df[[i]]$Trial,FUN=sum))
  
  meankLa.s = c(meankLa.s,mean.new.kLa)
  sdkLa.s = c(sdkLa.s,error.new.kLa)
}

kLa_Results.df$kLa = meankLa.s
kLa_Results.df$kLa.sd = sdkLa.s
kLa_Results.df$kLa.weight = 1/kLa_Results.df$kLa.sd^2


#kLa Modeling Results 
kLa_Model.df = data.frame(Temp.K=CO2_ExperimentalResults.df$Temp.C+273.15,
                          CO2.Angelo = CO2_ExperimentalResults.df$Angelo,
                          vL.L = c(rep(LiquidHoldup.3LPM.Model$volume.holdup.L,each=3),rep(LiquidHoldup.4.5LPM.Model$volume.holdup.L,each=3)),
                          flowrate.lps = rep(c(3.0/60,4.5/60),each=9),
                          height.m= rep(c(rep(c(0.53,0.68,0.98),each=3)),2))

kLa_Model.df$H = (1/(exp(-159.854+8741.68/kLa_Model.df$Temp.K+
                           21.6694*log(kLa_Model.df$Temp.K)-0.0011026*kLa_Model.df$Temp.K)))/55342
kLa_Model.df$H = kLa_Model.df$H*1000/44

kLa_Model.df$CSat = ((5+14.7)/14.696)/kLa_Model.df$H

kLa_Model.df$kLa.s = kLa_Model.df$flowrate.lps/kLa_Model.df$vL.L*log(kLa_Model.df$CSat/(kLa_Model.df$CSat-kLa_Model.df$CO2.Angelo))


kLa_Results.df$kLa.model = kLa_Model.df$kLa.s
kLa_Results.df$Temp.C = kLa_Model.df$Temp.K-273.15

#Composite kLa results 
kLa.overall.df = data.frame(Flowrate = rep(c(3.0,4.5),each=3),
                            Height = rep(c(0.53,0.68,0.98),2))

kLa.overall.df$Temp = as.vector(tapply(kLa_Results.df$Temp.C,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))
kLa.overall.df$Temp.sd = as.vector(tapply(kLa_Results.df$Temp.C,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=sd))
kLa.overall.df$model.kLa = as.vector(tapply(kLa_Results.df$kLa.model,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=mean))
kLa.overall.df$model.kLa.sd = as.vector(tapply(kLa_Results.df$kLa.model,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=sd))


kLa.overall.df$Experimental.kLa = as.vector(tapply(kLa_Results.df$kLa,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),function(x,w){
  weighted.mean(x,w[match(x,kLa_Results.df$kLa)])},w=kLa_Results.df$kLa.weight))
kLa.overall.df$Experimental.kLa.sd = as.vector(sqrt(1/tapply(kLa_Results.df$kLa.weight,list(kLa_Results.df$Height.m,kLa_Results.df$FlowRate.LPM),FUN=sum)))

