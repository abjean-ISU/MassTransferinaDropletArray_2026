#VelocityAnalysis.R
#Droplet velocity and acceleration from the tracked centroids in the Frames.csv files written by StackAnalysis.R.
#Answers R1-3 (drops are far from terminal velocity), R2-3 (velocity model), R2-4 (can slow drops explain
#the holdup gap?) and gives the drafting hint for R1-8.
#
#For every tracked drop: position vs time is fitted with a straight line (mean speed) and a quadratic
#(speed at the track midpoint and acceleration). The model acceleration depends only on drop size and
#current speed, a = g - (9/2) mu_g/(r^2 rho_l) v (1 + 0.158 Re^(2/3)) (same drag law as the MATLAB model),
#so measured and model acceleration can be compared without knowing where the imaging windows were.
#The measured speeds are then compared with the model speed at the estimated window distances.

rm(list=ls())
library(ggplot2)

################################################################################
## Settings
################################################################################
#one row per data set: Frames.csv from StackAnalysis.R and its flow rate
dataSets = data.frame(file = c('../../Experimental Data/Stack Analysis/09232024/09232024_Frames.csv',
                               '../../Experimental Data/Stack Analysis/09242024/09242024_Frames.csv'),
                      FlowRate.LPM = c(3.0,4.5))
outFolder = getwd()   #write outputs to the current working directory

fps = 3400                 #camera frame rate (frames/s)
trimEnds = 1               #frames dropped at each end of a track (centroid bias when a drop enters/leaves the field)
minFrames = 15             #shortest track (after trimming) used for the quadratic fit
maxResid.px = 1            #tracks whose centroids scatter more than this about the fit are dropped (merged/broken tracks)
maxSE.a = 3                #tracks with acceleration standard error above this (m/s2) are dropped

d.orifice.m = 0.8e-3       #orifice diameter (m)
N0 = 200                   #number of orifices
g = 9.8                    #m/s2
rho.l = 997                #water density (kg/m3)
rho.g = 1.184              #air density (kg/m3); imaging-rig drops fell through room air at ~25 C, 1 atm
mu.g = 1.86e-5             #air viscosity (Pa s), same value as the model's air runs

#estimated distance of the imaging windows below the manifold (author, from the setup photo and backlight mounting)
windowDist.m = list(top = 0.15, mid = c(0.35,0.375,0.40))
tauHeight.m = 0.98         #fall height for the residence-time sensitivity check

dir.create(outFolder,showWarnings=F)
outFile = function(name) file.path(outFolder,name)

################################################################################
## Model functions (drag law of the MATLAB model)
################################################################################
a.model = function(v,r) g-(9/2)*(mu.g/(r^2*rho.l))*v*(1+0.158*(2*rho.g*v*r/mu.g)^(2/3))
v.terminal = function(r) uniroot(function(v) a.model(v,r),c(1e-6,50))$root

#march the model from the orifice (speed v0) to distance x.end; returns speed and time at x.end
modelTrajectory = function(r,v0,x.end,dx=1e-3){
  x = 0; v = v0; t = 0
  while (x < x.end-1e-12){ t = t+dx/v; v = v+a.model(v,r)/v*dx; x = x+dx }
  c(v=v,t=t)
}

#distance a drop of radius r leaving at v0 must fall to reach speed v.target (NA if it never does within 5 m)
impliedDistance = function(v.target,r,v0,dx=1e-3,xmax=5){
  if (is.na(v.target)) return(NA)
  if (v.target<=v0) return(0)
  x = 0; v = v0
  while (v<v.target & x<xmax){ v = v+a.model(v,r)/v*dx; x = x+dx }
  if (v<v.target) NA else x
}

################################################################################
## Per-drop fits
################################################################################
fitDrop = function(d){
  d = d[order(d$FrameNo),]
  if (trimEnds>0 & nrow(d)>2*trimEnds) d = d[(trimEnds+1):(nrow(d)-trimEnds),]
  n = nrow(d)
  out = data.frame(Droplet.Designation = d$Droplet.Designation[1], StackID = d$StackID[1],
                   Height = d$Height[1], Zone = d$Zone[1], Depth = d$Depth[1],
                   nFrames = n, frameGaps = sum(diff(d$FrameNo)>1),
                   Diameter.mm = mean(d$Diameter.mm), scale.pxmm = median(d$yCentroid.px/d$yCentroid.mm),
                   v.lin = NA, v.mid = NA, a = NA, a.se = NA, resid.px = NA)
  if (n < minFrames) return(out)
  t = (d$FrameNo-d$FrameNo[1])/fps             #time (s); fitting on frame number handles missing frames
  y = d$yCentroid.mm/1000                       #position (m); image y increases downward
  quad = lm(y~t+I(t^2))
  cq = summary(quad)$coefficients
  out$v.lin = coef(lm(y~t))[2]                  #mean speed over the track (m/s)
  out$v.mid = cq[2,1]+2*cq[3,1]*mean(t)         #speed at the track's mean time (m/s)
  out$a = 2*cq[3,1]                             #acceleration (m/s2)
  out$a.se = 2*cq[3,2]                          #standard error of the acceleration (m/s2)
  out$resid.px = summary(quad)$sigma*1000*out$scale.pxmm   #centroid scatter about the fit (px)
  out
}

med.iqr = function(x) sprintf('%.2f (IQR %.2f to %.2f)',median(x,na.rm=T),quantile(x,0.25,na.rm=T),quantile(x,0.75,na.rm=T))

Summary.list = list(); Height.list = list(); Bins.list = list(); Window.list = list(); Drops.list = list()

for (k in 1:nrow(dataSets)){
  Q.LPM = dataSets$FlowRate.LPM[k]
  v0 = Q.LPM/1000/60/(N0*pi*(d.orifice.m/2)^2)  #orifice exit velocity (m/s)
  Frames.df = read.csv(dataSets$file[k],stringsAsFactors=F)
  Frames.df = Frames.df[!is.na(Frames.df$yCentroid.mm),]

  Drops.df = do.call(rbind,lapply(split(Frames.df,Frames.df$Droplet.Designation),fitDrop))
  rownames(Drops.df) = NULL
  if (median(Drops.df$v.lin,na.rm=T)<0) stop('Velocities are negative: image y increases upward; flip the sign of y.')

  fitted = !is.na(Drops.df$a)
  Drops.df$used = fitted & Drops.df$resid.px<maxResid.px & Drops.df$a.se<maxSE.a
  D = Drops.df[Drops.df$used,]
  D$r.m = D$Diameter.mm/2000
  D$a.model = mapply(a.model,D$v.mid,D$r.m)
  D$a.resid = D$a-D$a.model
  D$v.terminal = sapply(D$r.m,v.terminal)
  D$v.frac.terminal = D$v.mid/D$v.terminal
  D$x.implied.m = mapply(impliedDistance,D$v.mid,D$r.m,MoreArgs=list(v0=v0))

  #measured speed / model speed at the estimated window distance(s), one row per height and distance
  Window.list[[k]] = do.call(rbind,lapply(names(windowDist.m),function(h){
    s = D[D$Height==h,]
    do.call(rbind,lapply(windowDist.m[[h]],function(xw){
      v.mod = sapply(s$r.m,function(r) modelTrajectory(r,v0,xw)['v'])
      data.frame(FlowRate.LPM = Q.LPM, Height = h, WindowDist.m = xw, n = nrow(s),
                 v.measured = round(median(s$v.mid),3), v.model = round(median(v.mod),3),
                 ratio.measured.over.model = round(median(s$v.mid/v.mod),3))
    }))
  }))

  #residence-time sensitivity: model time to fall tauHeight.m for the median drop
  tau.med = modelTrajectory(median(D$r.m),v0,tauHeight.m)['t']

  Summary.list[[k]] = data.frame(
    FlowRate.LPM = Q.LPM, v0.ms = round(v0,3),
    Drops = nrow(Drops.df),
    Fitted = sum(fitted),
    Used = nrow(D),
    TracksWithGaps.pct = round(100*mean(Drops.df$frameGaps>0),1),
    MedianFrames = median(D$nFrames),
    MedianDuration.ms = round(median((D$nFrames-1)/fps*1000),1),
    CentroidScatter.px = round(median(D$resid.px),3),
    SE.a.median = round(median(D$a.se),3),
    Resolved.pct = round(100*mean(abs(D$a)>2*D$a.se),1),
    a.measured = med.iqr(D$a),
    a.model = med.iqr(D$a.model),
    a.measured.minus.model = med.iqr(D$a.resid),
    a.resid.weighted = sprintf('%.3f +/- %.3f',sum(D$a.resid/D$a.se^2)/sum(1/D$a.se^2),1/sqrt(sum(1/D$a.se^2))),
    v.mid = med.iqr(D$v.mid),
    v.over.terminal = round(median(D$v.frac.terminal),3),
    FasterThanModelAllows = sum(is.na(D$x.implied.m)),
    tau.model.median.drop.s = round(tau.med,3),
    tau.if.10pct.slower.s = round(tau.med/0.9,3))

  #by imaging height
  Height.list[[k]] = do.call(rbind,lapply(split(D,D$Height),function(s){
    data.frame(FlowRate.LPM = Q.LPM, Height = s$Height[1], n = nrow(s),
               v.mid = round(median(s$v.mid),3), a = round(median(s$a),3), a.model = round(median(s$a.model),3),
               x.implied.m = round(median(s$x.implied.m,na.rm=T),3))
  }))

  #by diameter class
  D$d.bin = cut(D$Diameter.mm,c(0,1.2,1.5,1.8,2.1,2.5,3.5,4.5,10))
  Bins.list[[k]] = do.call(rbind,lapply(split(D,D$d.bin,drop=T),function(s)
    data.frame(FlowRate.LPM = Q.LPM, d.bin = as.character(s$d.bin[1]), n = nrow(s),
               v.mid = round(median(s$v.mid),3), a = round(median(s$a),3), a.model = round(median(s$a.model),3),
               a.wmean = round(sum(s$a/s$a.se^2)/sum(1/s$a.se^2),3), a.wse = round(1/sqrt(sum(1/s$a.se^2)),3),
               x.implied.m = round(median(s$x.implied.m,na.rm=T),3))))

  D$FlowRate.LPM = Q.LPM
  Drops.list[[k]] = D

  #plots
  dGrid = seq(min(D$Diameter.mm),max(D$Diameter.mm),length.out=100)
  Term.df = data.frame(Diameter.mm=dGrid,v.terminal=sapply(dGrid/2000,v.terminal))
  p.speed = ggplot(D,aes(x=Diameter.mm,y=v.mid,color=Height))+theme_bw()+
    geom_point(alpha=0.3)+
    geom_line(data=Term.df,aes(x=Diameter.mm,y=v.terminal),color='black',linewidth=1,inherit.aes=F)+
    geom_hline(yintercept=v0,linetype='dashed')+
    labs(x='Droplet Diameter (mm)',y='Droplet Velocity (m/s)',title=paste(Q.LPM,'L/min: speed vs terminal velocity (line)'))
  ggsave(outFile(paste0(Q.LPM,'LPM_Speed_vs_Diameter.png')),p.speed,width=8,height=5,dpi=300)

  dShow = quantile(D$Diameter.mm,c(0.1,0.5,0.9))
  vGrid = seq(0,max(D$v.mid)*1.1,length.out=100)
  Curve.df = do.call(rbind,lapply(dShow,function(d) data.frame(v=vGrid,a=a.model(vGrid,d/2000),
                                                              d.label=sprintf('model, d = %.2f mm',d))))
  p.accel = ggplot(D,aes(x=v.mid,y=a))+theme_bw()+
    geom_errorbar(aes(ymin=a-a.se,ymax=a+a.se),alpha=0.15,width=0)+
    geom_point(aes(color=Diameter.mm),alpha=0.5)+
    geom_line(data=Curve.df,aes(x=v,y=a,linetype=d.label),linewidth=1)+
    labs(x='Droplet Velocity (m/s)',y=bquote(Acceleration~(m/s^2)),color='d (mm)',linetype=NULL,
         title=paste(Q.LPM,'L/min: measured acceleration vs model drag law'))+
    theme(legend.position='bottom',legend.box='vertical')
  ggsave(outFile(paste0(Q.LPM,'LPM_Accel_vs_Speed.png')),p.accel,width=8,height=6,dpi=300)

  p.dist = ggplot(D,aes(x=x.implied.m,fill=Height))+theme_bw()+
    geom_histogram(bins=40,position='identity',alpha=0.5)+
    geom_vline(xintercept=unlist(windowDist.m),linetype='dashed')+
    labs(x='Model-Implied Fall Distance (m)',y='Drops',title=paste(Q.LPM,'L/min (dashed: estimated window distances)'))
  ggsave(outFile(paste0(Q.LPM,'LPM_ImpliedFallDistance.png')),p.dist,width=8,height=5,dpi=300)
}

################################################################################
## Results
################################################################################
Summary.df = do.call(rbind,Summary.list)
Height.df = do.call(rbind,Height.list)
Window.df = do.call(rbind,Window.list)
Bins.df = do.call(rbind,Bins.list)
rownames(Bins.df) = NULL

SummaryT = t(Summary.df); colnames(SummaryT) = paste0(Summary.df$FlowRate.LPM,' L/min')
cat('\n==== Summary (one column per flow rate) ====\n'); print(SummaryT,quote=F)
cat('\n==== By imaging height ====\n'); print(Height.df,row.names=F)
cat('\n==== Measured vs model speed at the estimated window distances ====\n'); print(Window.df,row.names=F)
cat('\n==== By diameter class ====\n'); print(Bins.df,row.names=F)

write.csv(SummaryT,outFile('Velocity_Summary.csv'))
write.csv(Height.df,outFile('Velocity_ByHeight.csv'),row.names=F)
write.csv(Window.df,outFile('Velocity_WindowComparison.csv'),row.names=F)
write.csv(Bins.df,outFile('Velocity_ByDiameter.csv'),row.names=F)
allDrops = do.call(rbind,lapply(Drops.list,function(x){x$d.bin=as.character(x$d.bin); x}))
write.csv(allDrops,outFile('Velocity_PerDrop.csv'),row.names=F)
cat('\nOutputs written to',outFolder,'\n')
