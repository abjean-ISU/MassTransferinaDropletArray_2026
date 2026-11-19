#AspectRatioAnalysis.R
#Drop shape oscillation from the tracked drops in the Frames.csv files written by StackAnalysis.R (R2-7).
#Questions answered:
#  1. How large is the shape deformation (aspect ratio) of the drops in the imaging windows, and is it
#     above the pixel noise floor?
#  2. Is there a periodic shape oscillation, and is its frequency close to the Rayleigh n = 2 mode
#     f2 = (1/2pi) sqrt(8 sigma/(rho r^3)) that the oscillating-drop correlations (Angelo, Hsu) assume?
#  3. How often do drops come near the circularity filter limit (aspect ratio 2.21, SupportingCalculations.R §2)?
#  4. Does the deformation change between the two imaging heights (damping)?
#Aspect ratio AR = xDiameter/yDiameter from the bounding box (AR > 1: wider than tall). The analysis uses
#log(AR), which treats wider and taller shapes symmetrically.
#Frequency is estimated two ways: a Lomb-Scargle periodogram of each long track (handles missing frames),
#and the autocorrelation of log(AR) pooled over all tracks in a diameter class (first minimum = half period).

rm(list=ls())
library(ggplot2)

################################################################################
## Settings
################################################################################
dataSets = data.frame(file = c('../../Experimental Data/Stack Analysis/09232024/09232024_Frames.csv',
                               '../../Experimental Data/Stack Analysis/09242024/09242024_Frames.csv'),
                      FlowRate.LPM = c(3.0,4.5))
outFolder = getwd()   #write outputs to the current working directory

fps = 3400                 #camera frame rate (frames/s)
trimEnds = 2               #frames dropped at each end of a track (drop partly outside the field distorts the box)
minFrames = 10             #shortest track (after trimming) used for the amplitude statistics
minFramesFreq = 30         #shortest track for a per-drop frequency (about 1.5 periods of a 1.6 mm drop)
maxLag = 40                #autocorrelation lags (frames)
fGrid = seq(30,1200,by=2)  #periodogram frequencies (Hz); Nyquist is fps/2 = 1700 Hz
dBreaks = c(0,1.3,1.7,2.1,2.5,3.5,6) #diameter classes (mm) for the pooled autocorrelation
AR.filter = 2.21           #aspect ratio at the circularity limit 0.8 (SupportingCalculations.R §2)

T.C = 22                   #water temperature for surface tension (imaging rig, room temperature)
rho.l = 997                #kg/m3
mu.l = 0.95e-3             #water viscosity at 22 C (Pa s)
sigma = function(T.C){T = T.C+273.15; 235.8e-3*(1-T/647.15)^1.256*(1-0.625*(1-T/647.15))} #N/m, as in the model
f.rayleigh = function(r) sqrt(8*sigma(T.C)/(rho.l*r^3))/(2*pi)   #n = 2 mode (Hz), r in m
tau.damp = function(r) rho.l*r^2/(5*mu.l)                         #Lamb viscous damping time of n = 2 (s)

#imaging-rig drag law and window distances (same as VelocityAnalysis.R) for the decay estimate
g = 9.8; rho.g = 1.184; mu.g = 1.86e-5; d.orifice.m = 0.8e-3; N0 = 200
windowDist.m = c(top = 0.15, mid = 0.375)

dir.create(outFolder,showWarnings=F)
outFile = function(name) file.path(outFolder,name)

################################################################################
## Helper functions
################################################################################
#Lomb-Scargle normalized periodogram (Press et al., Numerical Recipes 13.8)
lombScargle = function(t,y,f){
  y = y-mean(y); s2 = var(y)
  sapply(f,function(ff){
    w = 2*pi*ff
    tau = atan2(sum(sin(2*w*t)),sum(cos(2*w*t)))/(2*w)
    cc = cos(w*(t-tau)); ss = sin(w*(t-tau))
    (sum(y*cc)^2/sum(cc^2)+sum(y*ss)^2/sum(ss^2))/(2*s2)
  })
}

#time for a drop of radius r leaving the orifice at v0 to reach distance x.end (air, imaging rig)
a.model = function(v,r) g-(9/2)*(mu.g/(r^2*rho.l))*v*(1+0.158*(2*rho.g*v*r/mu.g)^(2/3))
timeTo = function(r,v0,x.end,dx=1e-3){ x = 0; v = v0; t = 0
  while (x < x.end-1e-12){ t = t+dx/v; v = v+a.model(v,r)/v*dx; x = x+dx }; t }

med.iqr = function(x) sprintf('%.3f (IQR %.3f to %.3f)',median(x,na.rm=T),quantile(x,0.25,na.rm=T),quantile(x,0.75,na.rm=T))

################################################################################
## Per-drop statistics
################################################################################
fitDrop = function(d){
  d = d[order(d$FrameNo),]
  if (trimEnds>0 & nrow(d)>2*trimEnds) d = d[(trimEnds+1):(nrow(d)-trimEnds),]
  n = nrow(d)
  out = data.frame(Droplet.Designation = d$Droplet.Designation[1], StackID = d$StackID[1],
                   Height = d$Height[1], nFrames = n, Diameter.mm = mean(d$Diameter.mm),
                   AR.mean = NA, AR.max = NA, logAR.sd = NA, noise.sd = NA, sd.over.noise = NA,
                   f.peak = NA, power.peak = NA, f.rayleigh = NA, f.ratio = NA)
  if (n < minFrames) return(out)
  t = (d$FrameNo-d$FrameNo[1])/fps
  y = log(d$AR)
  yd = residuals(lm(y~t))                       #remove slow drift (e.g. slight focus change)
  out$AR.mean = exp(mean(y))
  out$AR.max = exp(max(abs(y)))                 #largest deformation either way, as an aspect ratio >= 1
  out$logAR.sd = sd(yd)
  out$noise.sd = median(d$noise.logAR)
  out$sd.over.noise = out$logAR.sd/out$noise.sd
  r = out$Diameter.mm/2000
  out$f.rayleigh = f.rayleigh(r)
  if (n >= minFramesFreq && sd(yd) > 0){
    P = lombScargle(t,yd,fGrid)
    out$f.peak = fGrid[which.max(P)]
    out$power.peak = max(P)
    out$f.ratio = out$f.peak/out$f.rayleigh
  }
  out
}

Summary.list = list(); Drops.list = list(); ACF.list = list(); ACFmin.list = list(); Height.list = list()

for (k in 1:nrow(dataSets)){
  Q.LPM = dataSets$FlowRate.LPM[k]
  v0 = Q.LPM/1000/60/(N0*pi*(d.orifice.m/2)^2)
  Frames.df = read.csv(dataSets$file[k],stringsAsFactors=F)
  Frames.df = Frames.df[!is.na(Frames.df$xDiameter.mm) & !is.na(Frames.df$yDiameter.mm),]
  if (is.null(Frames.df$AR)) Frames.df$AR = Frames.df$xDiameter.mm/Frames.df$yDiameter.mm
  if (!is.null(Frames.df$AspectRatio)) Frames.df$AR = Frames.df$AspectRatio

  #pixel noise floor of log(AR) for each frame: bounding-box width and height each uncertain by
  #+/-0.5 px (rectangular, SD 1/sqrt(12) px) -> SD(log AR) = sqrt(1/dx^2+1/dy^2)/sqrt(12)
  if (all(c('xDiameter.px','yDiameter.px') %in% names(Frames.df)) && !all(is.na(Frames.df$xDiameter.px))){
    dx = Frames.df$xDiameter.px; dy = Frames.df$yDiameter.px
  } else {
    sc = Frames.df$yCentroid.px/Frames.df$yCentroid.mm          #px/mm
    dx = Frames.df$xDiameter.mm*sc; dy = Frames.df$yDiameter.mm*sc
  }
  Frames.df$noise.logAR = sqrt(1/dx^2+1/dy^2)/sqrt(12)

  Drops.df = do.call(rbind,lapply(split(Frames.df,Frames.df$Droplet.Designation),fitDrop))
  rownames(Drops.df) = NULL
  D = Drops.df[!is.na(Drops.df$AR.mean),]
  D$FlowRate.LPM = Q.LPM
  D$d.bin = cut(D$Diameter.mm,dBreaks)
  Drops.list[[k]] = D

  #frames used for the frame-level statistics: the trimmed frames of the drops in D
  usedFrames = do.call(rbind,lapply(split(Frames.df,Frames.df$Droplet.Designation),function(d){
    d = d[order(d$FrameNo),]
    if (trimEnds>0 & nrow(d)>2*trimEnds) d = d[(trimEnds+1):(nrow(d)-trimEnds),]
    d}))
  usedFrames = usedFrames[usedFrames$Droplet.Designation %in% D$Droplet.Designation,]
  ARsym = exp(abs(log(usedFrames$AR)))          #deformation as an aspect ratio >= 1

  #expected decay of the n = 2 oscillation between breakup and each window (median drop)
  r.med = median(D$Diameter.mm)/2000
  decay = sapply(windowDist.m,function(x) exp(-timeTo(r.med,v0,x)/tau.damp(r.med)))

  Long = D[!is.na(D$f.peak),]
  Summary.list[[k]] = data.frame(
    FlowRate.LPM = Q.LPM,
    Drops = nrow(Drops.df), Used = nrow(D), Frames = nrow(usedFrames),
    MedianFrames = median(D$nFrames),
    AR.mean.median = med.iqr(D$AR.mean),
    Frames.AR.over.1.1.pct = round(100*mean(ARsym>1.1),2),
    Frames.AR.over.1.2.pct = round(100*mean(ARsym>1.2),2),
    Frames.AR.over.1.5.pct = round(100*mean(ARsym>1.5),2),
    Frames.AR.over.2.0.pct = round(100*mean(ARsym>2.0),3),
    Frames.AR.over.filter.pct = round(100*mean(ARsym>AR.filter),3),
    logAR.sd.within.track = med.iqr(D$logAR.sd),
    noise.floor.sd = med.iqr(D$noise.sd),
    Tracks.sd.above.2x.noise.pct = round(100*mean(D$sd.over.noise>2,na.rm=T),1),
    Tracks.for.frequency = nrow(Long),
    f.peak.Hz = if (nrow(Long)>0) med.iqr(Long$f.peak) else NA,
    f.rayleigh.Hz = if (nrow(Long)>0) med.iqr(Long$f.rayleigh) else NA,
    f.peak.over.rayleigh = if (nrow(Long)>0) med.iqr(Long$f.ratio) else NA,
    Tracks.within.20pct.of.rayleigh.pct = if (nrow(Long)>0) round(100*mean(abs(Long$f.ratio-1)<0.2),1) else NA,
    damping.time.median.drop.s = round(tau.damp(r.med),3),
    amplitude.left.at.top.window = round(decay['top'],3),
    amplitude.left.at.mid.window = round(decay['mid'],3))

  #by imaging height
  Height.list[[k]] = do.call(rbind,lapply(split(D,D$Height),function(s)
    data.frame(FlowRate.LPM = Q.LPM, Height = s$Height[1], n = nrow(s),
               d.median.mm = round(median(s$Diameter.mm),3),
               AR.mean.median = round(median(s$AR.mean),4),
               logAR.sd.median = round(median(s$logAR.sd),4),
               noise.sd.median = round(median(s$noise.sd),4))))

  #pooled autocorrelation of detrended log(AR) by diameter class (pairs of frames in the same track)
  for (b in levels(droplevels(D$d.bin))){
    ids = D$Droplet.Designation[D$d.bin==b]
    num = numeric(maxLag+1); cnt = numeric(maxLag+1)
    for (id in ids){
      d = Frames.df[Frames.df$Droplet.Designation==id,]
      d = d[order(d$FrameNo),]
      if (trimEnds>0 & nrow(d)>2*trimEnds) d = d[(trimEnds+1):(nrow(d)-trimEnds),]
      t = d$FrameNo-d$FrameNo[1]
      y = residuals(lm(log(d$AR)~t))
      full = rep(NA,max(t)+1); full[t+1] = y
      for (L in 0:min(maxLag,length(full)-1)){
        p = full[1:(length(full)-L)]*full[(1+L):length(full)]
        num[L+1] = num[L+1]+sum(p,na.rm=T); cnt[L+1] = cnt[L+1]+sum(!is.na(p))
      }
    }
    acf = (num/cnt)/(num[1]/cnt[1])
    dmed = median(D$Diameter.mm[D$d.bin==b])
    fR = f.rayleigh(dmed/2000)
    ACF.list[[length(ACF.list)+1]] = data.frame(FlowRate.LPM = Q.LPM, d.bin = b, n.drops = length(ids),
                                                 lag.frames = 0:maxLag, lag.ms = 1000*(0:maxLag)/fps,
                                                 acf = acf, pairs = cnt)
    #first local minimum after lag 1 (half a period if the oscillation is periodic)
    iMin = which(diff(sign(diff(acf)))>0)+1
    iMin = iMin[iMin>2]
    lagMin = if (length(iMin)>0) iMin[1]-1 else NA
    ACFmin.list[[length(ACFmin.list)+1]] = data.frame(FlowRate.LPM = Q.LPM, d.bin = b, n.drops = length(ids),
      d.median.mm = round(dmed,3), f.rayleigh.Hz = round(fR),
      expected.half.period.frames = round(fps/(2*fR),1),
      acf.first.min.frames = lagMin,
      acf.at.first.min = if (is.na(lagMin)) NA else round(acf[lagMin+1],3),
      f.from.acf.Hz = if (is.na(lagMin)) NA else round(fps/(2*lagMin)),
      acf.lag1 = round(acf[2],3))
  }

  #plots: longest tracks
  show = head(D[order(-D$nFrames),'Droplet.Designation'],6)
  S = Frames.df[Frames.df$Droplet.Designation %in% show,]
  S = do.call(rbind,lapply(split(S,S$Droplet.Designation),function(d){d$t.ms = 1000*(d$FrameNo-min(d$FrameNo))/fps
    d$label = sprintf('%s (d = %.2f mm, f2 = %.0f Hz)',d$Droplet.Designation[1],mean(d$Diameter.mm),f.rayleigh(mean(d$Diameter.mm)/2000)); d}))
  p.ts = ggplot(S,aes(x=t.ms,y=AR))+theme_bw()+geom_line()+geom_point(size=0.8)+
    geom_hline(yintercept=1,linetype='dashed')+facet_wrap(~label,ncol=2,scales='free_x')+
    labs(x='Time (ms)',y='Aspect Ratio (x/y)',title=paste(Q.LPM,'L/min: longest tracks'))
  ggsave(outFile(paste0(Q.LPM,'LPM_AR_TimeSeries.png')),p.ts,width=10,height=8,dpi=300)

  p.amp = ggplot(D,aes(x=Diameter.mm,y=logAR.sd))+theme_bw()+
    geom_point(aes(color=Height),alpha=0.4)+geom_point(aes(y=noise.sd),shape=3,alpha=0.3)+
    labs(x='Droplet Diameter (mm)',y='SD of log(AR) within track',
         title=paste(Q.LPM,'L/min: shape fluctuation (points) vs pixel noise floor (+)'))
  ggsave(outFile(paste0(Q.LPM,'LPM_AR_Amplitude_vs_Diameter.png')),p.amp,width=8,height=5,dpi=300)

  if (nrow(Long)>0){
    dGrid = seq(min(Long$Diameter.mm),max(Long$Diameter.mm),length.out=100)
    p.f = ggplot(Long,aes(x=Diameter.mm,y=f.peak))+theme_bw()+geom_point(aes(size=power.peak),alpha=0.4)+
      geom_line(data=data.frame(Diameter.mm=dGrid,f=f.rayleigh(dGrid/2000)),aes(x=Diameter.mm,y=f),inherit.aes=F)+
      labs(x='Droplet Diameter (mm)',y='Periodogram Peak (Hz)',size='Peak power',
           title=paste(Q.LPM,'L/min: dominant shape frequency vs Rayleigh n = 2 (line)'))
    ggsave(outFile(paste0(Q.LPM,'LPM_AR_Frequency_vs_Diameter.png')),p.f,width=8,height=5,dpi=300)
  }
}

################################################################################
## Results
################################################################################
Summary.df = do.call(rbind,Summary.list)
Height.df = do.call(rbind,Height.list)
ACF.df = do.call(rbind,ACF.list)
ACFmin.df = do.call(rbind,ACFmin.list)
allDrops = do.call(rbind,lapply(Drops.list,function(x){x$d.bin=as.character(x$d.bin); x}))

SummaryT = t(Summary.df); colnames(SummaryT) = paste0(Summary.df$FlowRate.LPM,' L/min')
cat('\n==== Summary (one column per flow rate) ====\n'); print(SummaryT,quote=F)
cat('\n==== By imaging height ====\n'); print(Height.df,row.names=F)
cat('\n==== Pooled autocorrelation: first minimum vs expected half period ====\n'); print(ACFmin.df,row.names=F)

p.acf = ggplot(ACF.df,aes(x=lag.ms,y=acf,color=d.bin))+theme_bw()+geom_line()+geom_hline(yintercept=0)+
  facet_wrap(~FlowRate.LPM,labeller=label_both)+
  labs(x='Lag (ms)',y='Autocorrelation of log(AR)',color='d (mm)',title='Pooled shape autocorrelation by diameter class')
ggsave(outFile('AR_Autocorrelation.png'),p.acf,width=10,height=5,dpi=300)

write.csv(SummaryT,outFile('AR_Summary.csv'))
write.csv(Height.df,outFile('AR_ByHeight.csv'),row.names=F)
write.csv(ACFmin.df,outFile('AR_AutocorrelationMinimum.csv'),row.names=F)
write.csv(ACF.df,outFile('AR_Autocorrelation.csv'),row.names=F)
write.csv(allDrops,outFile('AR_PerDrop.csv'),row.names=F)
cat('\nOutputs written to',outFolder,'\n')
