rm(list=ls())
library(ggplot2)

################################################################################
# Per-stack droplet size analysis for ONE imaging data set (one date / flow rate)
#
# Run once per data set; each run reads only its own folders and writes to its
# own output folder.
#
# Input: the Droplet Data Checks folder. Each CSV there is the manually checked
# copy of one stack's dropData.csv, and its file name is the stack ID (same name
# as the image-stack folder, e.g. 09232024_backcenter_.8gpm_mid.csv). Bad
# droplets are removed by deleting their rows, so every row is an accepted droplet.
#
# Serves: R1-7 / R2-9 (per-location distributions), R1-6 (front/back
# comparison), Section 3C Step 1b (spatial breakdown, dripping fraction four
# ways). Also writes an AllDropletData file for the MATLAB model.
################################################################################

## CHANGE FILES ##
dataName = '09242024'   #label used in output file names
FlowRate.LPM = 4.5
checkFolder = '../Experimental Data/Droplet Data Checks'
#folder holding this date's image-stack folders (each has the .tiff frames)
stackRoot = 'path/to/09242024 image stacks'   #raw .tiff image-stack folders (not included in Data Sharing)
outFolder = '../Experimental Data/Stack Analysis/09242024'
locationFile = file.path(outFolder,paste0(dataName,'_StackLocations.csv')) #filled in by hand (see Section 3)

drip.cutoff.mm = 3.0   #dripping/jetting cutoff (valley between modes); check against the histogram
binwidth.mm = 0.1      #histogram bin width

dir.create(outFolder,showWarnings=F,recursive=T)
outFile = function(name) file.path(outFolder,paste0(dataName,'_',name))

################################################################################
## 1. Read every stack (one row per droplet per frame)
################################################################################
#columns needed from dropData.csv (names as read.csv converts them)
needCols = c('Frame.No','Droplet.No','X.Diameter..mm.','Y.Diameter..mm.')
optCols = c('X.Centroid..Pixels.','Y.Centroid..Pixels.','X.Centroid..mm.','Y.Centroid..mm.',
            'X.Diameter..Pixels.','Y.Diameter..Pixels.','Gradient','Darkness')

#reads one dropData-style csv; a stack with no accepted droplets (empty file,
#header only, or blank rows) returns a table with the right columns and 0 rows
readDrops = function(file){
  empty = data.frame(matrix(numeric(0),nrow=0,ncol=length(c(needCols,optCols)),
                            dimnames=list(NULL,c(needCols,optCols))))
  if (file.size(file)==0) return(empty)
  df = tryCatch(read.csv(file),error=function(e) empty)
  if (nrow(df)==0) return(empty)
  missing = setdiff(needCols,names(df))
  if (length(missing)>0) stop(paste(basename(file),'is missing column(s):',paste(missing,collapse=', ')))
  df = df[!is.na(df$Frame.No),] #drops blank rows left after deleting in Excel
  for (col in setdiff(optCols,names(df))) df[[col]] = NA #older files may lack centroids/gradient
  df[,c(needCols,optCols)]
}

files = list.files(path=checkFolder,pattern='\\.csv$',full.names=T)
if (length(files)==0) stop(paste('No csv files found in',checkFolder))
#stack ID = image-stack folder name: drop '.csv', the 'dropData.' prefix and the '_CHECK' suffix
stackIDof = function(f){
  id = sub('\\.csv$','',basename(f),ignore.case=T)
  id = sub('^dropData[._]','',id,ignore.case=T) #prefix may be 'dropData.' or 'dropData_'
  id = sub('_CHECK$','',id,ignore.case=T)
  sub('_0\\.','_.',id)
}
allStacks = sort(stackIDof(files)) #every stack, including ones with no droplets

frames.list = list()
for (file in files){
  df = readDrops(file)
  if (nrow(df)==0) next
  df$StackID = stackIDof(file)
  frames.list[[length(frames.list)+1]] = df
}
if (length(frames.list)==0) stop('Every stack file is empty')
Frames.df = do.call(rbind,frames.list)

emptyStacks = setdiff(allStacks,unique(Frames.df$StackID))
if (length(emptyStacks)>0) cat('Stacks with no accepted droplets (kept in the summary with n = 0):',
                               paste(emptyStacks,collapse=', '),'\n')

#shorter column names
names(Frames.df)[match(c(needCols,optCols),names(Frames.df))] =
  c('FrameNo','DropletNo','xDiameter.mm','yDiameter.mm','xCentroid.px','yCentroid.px',
    'xCentroid.mm','yCentroid.mm','xDiameter.px','yDiameter.px','Gradient','Darkness')
Frames.df = Frames.df[order(Frames.df$StackID,Frames.df$DropletNo,Frames.df$FrameNo),]

#number of frames in each stack = number of .tiff files in its image folder
stackFrames = data.frame(StackID = allStacks)
stackFrames$nFramesStack = sapply(file.path(stackRoot,stackFrames$StackID),function(p){
  if (!dir.exists(p)) return(NA)
  length(list.files(path=p,pattern='\\.tiff?$',ignore.case=T))
})
noFolder = stackFrames$StackID[is.na(stackFrames$nFramesStack)]
if (length(noFolder)>0) warning(paste('Image folder not found for',length(noFolder),'stack(s); frame counts left blank:',
                                      paste(noFolder,collapse=', ')))
#a tracked frame number above the tiff count means the wrong folder was matched
lastTracked = tapply(Frames.df$FrameNo,Frames.df$StackID,FUN=max)
badCount = stackFrames$StackID[which(lastTracked[stackFrames$StackID] > stackFrames$nFramesStack)]

if (length(badCount)>0) warning(paste('Tracked frame numbers exceed the tiff count for:',paste(badCount,collapse=', ')))

#equivalent diameter per frame, same definition as the MATLAB model: (dx^2*dy)^(1/3)
Frames.df$Diameter.mm = (Frames.df$xDiameter.mm^2*Frames.df$yDiameter.mm)^(1/3)
Frames.df$AspectRatio = Frames.df$xDiameter.mm/Frames.df$yDiameter.mm
Frames.df$Droplet.Designation = paste(Frames.df$StackID,Frames.df$DropletNo,sep='.')

################################################################################
## 2. One row per unique droplet
################################################################################

byDrop = split(Frames.df,Frames.df$Droplet.Designation)
Droplets.df = do.call(rbind,lapply(byDrop,function(d){
  data.frame(Droplet.Designation = d$Droplet.Designation[1],
             StackID = d$StackID[1],
             DropletNo = d$DropletNo[1],
             nFrames = nrow(d),
             FirstFrame = min(d$FrameNo),
             LastFrame = max(d$FrameNo),
             Diameter.mm = mean(d$Diameter.mm),
             Diameter.CV.track = ifelse(nrow(d)>1,sd(d$Diameter.mm)/mean(d$Diameter.mm),NA),
             AspectRatio = mean(d$AspectRatio))
}))
rownames(Droplets.df) = NULL
Droplets.df$Volume.mm3 = pi/6*Droplets.df$Diameter.mm^3
Droplets.df$Dripping = Droplets.df$Diameter.mm > drip.cutoff.mm

#quality checks on tracks
cat('\n----',dataName,'(',FlowRate.LPM,'L/min ) track checks ----\n')
cat('Stacks:',nrow(stackFrames),'  Unique droplets:',nrow(Droplets.df),'\n')
cat('Droplets seen in only one frame:',sum(Droplets.df$nFrames==1),'\n')
#a large diameter spread within one track can mean two drops were merged in the manual check
suspect = Droplets.df[!is.na(Droplets.df$Diameter.CV.track) & Droplets.df$Diameter.CV.track>0.10,]
cat('Tracks with within-track diameter CV > 10% (check for merged drops):',nrow(suspect),'\n')
if (nrow(suspect)>0) print(head(suspect[order(-suspect$Diameter.CV.track),
                                        c('Droplet.Designation','nFrames','Diameter.mm','Diameter.CV.track')],20))

################################################################################
## 2b. What the manual check removed (original dropData.csv vs Data Checks copy)
################################################################################
# Each stack folder still holds the original MATLAB export (dropData.csv).
# Rows are matched on frame number, centroid and diameters, so droplets that
# were renumbered during the check still match. Diameters of deleted droplets
# are as MATLAB measured them; some were deleted because that measurement was
# bad, so treat this as a check for size-selective deletion, not true sizes.

keyOf = function(df){
  if (all(!is.na(df$xCentroid.px))){
    paste(df$FrameNo,round(df$xCentroid.px,2),round(df$yCentroid.px,2),
          round(df$xDiameter.mm,3),round(df$yDiameter.mm,3))
  } else {
    paste(df$FrameNo,round(df$xDiameter.mm,3),round(df$yDiameter.mm,3))
  }
}

orig.list = list()
origEmpty = character()
for (id in stackFrames$StackID){
  origFile = list.files(path=file.path(stackRoot,id),pattern='^dropData\\.csv$',ignore.case=T,full.names=T)
  if (length(origFile)==0) next
  od = readDrops(origFile[1])
  if (nrow(od)==0){ origEmpty = c(origEmpty,id); next } #MATLAB found no in-focus droplets
  names(od) = c('FrameNo','DropletNo','xDiameter.mm','yDiameter.mm','xCentroid.px','yCentroid.px',
                'xCentroid.mm','yCentroid.mm','xDiameter.px','yDiameter.px','Gradient','Darkness')
  od$StackID = id
  checked = Frames.df[Frames.df$StackID==id,]
  useCentroid = all(!is.na(od$xCentroid.px)) & all(!is.na(checked$xCentroid.px))
  if (!useCentroid){ od$xCentroid.px = NA; checked$xCentroid.px = NA }
  od$Kept = keyOf(od) %in% keyOf(checked)
  extra = sum(!(keyOf(checked) %in% keyOf(od)))
  if (extra>0) warning(paste(id,': ',extra,' checked row(s) not found in the original dropData.csv (edited values?)',sep=''))
  orig.list[[id]] = od
}
noOrig = setdiff(stackFrames$StackID,c(names(orig.list),origEmpty))
if (length(noOrig)>0) warning(paste('Original dropData.csv not found for',length(noOrig),'stack(s):',paste(noOrig,collapse=', ')))

Deletion.Summary = NULL
OrigDroplets.df = NULL
if (length(orig.list)>0){
  Orig.df = do.call(rbind,orig.list)
  Orig.df$Diameter.mm = (Orig.df$xDiameter.mm^2*Orig.df$yDiameter.mm)^(1/3)

  #one row per ORIGINAL droplet: kept whole, trimmed (some frames removed) or deleted
  OrigDroplets.df = do.call(rbind,lapply(split(Orig.df,paste(Orig.df$StackID,Orig.df$DropletNo)),function(d){
    data.frame(StackID = d$StackID[1],
               DropletNo.Original = d$DropletNo[1],
               nFrames = nrow(d),
               nFramesKept = sum(d$Kept),
               Diameter.mm = mean(d$Diameter.mm),
               Status = ifelse(all(d$Kept),'kept',ifelse(any(d$Kept),'trimmed','deleted')))
  }))
  rownames(OrigDroplets.df) = NULL
  OrigDroplets.df$Dripping = OrigDroplets.df$Diameter.mm > drip.cutoff.mm

  Deletion.Summary = do.call(rbind,lapply(split(OrigDroplets.df,OrigDroplets.df$StackID),function(s){
    id = s$StackID[1]
    o = orig.list[[id]]
    del = s[s$Status=='deleted',]
    kep = s[s$Status!='deleted',]
    data.frame(StackID = id,
               OrigRows = nrow(o),
               RowsDeleted = sum(!o$Kept),
               OrigDroplets = nrow(s),
               DropletsDeleted = nrow(del),
               DropletsTrimmed = sum(s$Status=='trimmed'),
               CheckedDroplets = length(unique(Frames.df$DropletNo[Frames.df$StackID==id])), #differs if tracks were merged/split
               pct.DropletsDeleted = 100*nrow(del)/nrow(s),
               mean.d.kept = mean(kep$Diameter.mm),
               mean.d.deleted = ifelse(nrow(del)>0,mean(del$Diameter.mm),NA),
               drip.num.kept = mean(kep$Dripping),
               drip.num.deleted = ifelse(nrow(del)>0,mean(del$Dripping),NA))
  }))
  rownames(Deletion.Summary) = NULL

  del = OrigDroplets.df[OrigDroplets.df$Status=='deleted',]
  kep = OrigDroplets.df[OrigDroplets.df$Status!='deleted',]
  volFrac = function(s) sum(s$Diameter.mm[s$Dripping]^3)/sum(s$Diameter.mm^3)
  cat('\n---- Manual check: what was removed ----\n')
  cat('Original droplets:',nrow(OrigDroplets.df),' deleted:',nrow(del),
      paste0('(',round(100*nrow(del)/nrow(OrigDroplets.df),1),'%)'),
      ' trimmed:',sum(OrigDroplets.df$Status=='trimmed'),'\n')
  cat('Original rows:',nrow(Orig.df),' deleted:',sum(!Orig.df$Kept),'\n')
  cat('Mean diameter: kept',round(mean(kep$Diameter.mm),3),'mm, deleted',round(mean(del$Diameter.mm),3),'mm\n')
  cat('Above dripping cutoff (by number): kept',round(100*mean(kep$Dripping),1),'%, deleted',round(100*mean(del$Dripping),1),'%\n')
  cat('Share of all large (>',drip.cutoff.mm,'mm) droplets that were deleted:',
      round(100*sum(del$Dripping)/max(1,sum(OrigDroplets.df$Dripping)),1),'%\n')
  cat('Dripping volume fraction: kept only',round(volFrac(kep),3),', if deleted droplets were added back',
      round(volFrac(OrigDroplets.df),3),'\n')
  if (nrow(del)>=5 & nrow(kep)>=5){
    w = wilcox.test(del$Diameter.mm,kep$Diameter.mm,exact=F)
    cat('Wilcoxon test, deleted vs kept diameters: p =',signif(w$p.value,3),'\n')
  }
  write.csv(Deletion.Summary,outFile('Deletion_Summary.csv'),row.names=F)
  write.csv(OrigDroplets.df,outFile('OriginalDroplets_Status.csv'),row.names=F)
}

################################################################################
## 2c. Motion of deleted vs kept objects (splash fragments: smaller, lateral, faster)
################################################################################
# Lateral movers cannot be linked across frames (tracking gate < 1 px sideways),
# so they appear as single-frame objects. Speed = vertical travel per frame (mm),
# from objects tracked over 2+ frames.

if (length(orig.list)>0){
  motion = do.call(rbind,lapply(split(Orig.df,paste(Orig.df$StackID,Orig.df$DropletNo)),function(d){
    d = d[order(d$FrameNo),]
    data.frame(StackID = d$StackID[1],
               DropletNo.Original = d$DropletNo[1],
               dy.mm.per.frame = ifelse(nrow(d)>1,median(diff(d$yCentroid.mm)/diff(d$FrameNo)),NA))
  }))
  OrigDroplets.df = merge(OrigDroplets.df,motion,by=c('StackID','DropletNo.Original'))
  OrigDroplets.df$Check = ifelse(OrigDroplets.df$Status=='deleted','deleted','kept')
  
  Motion.Summary = do.call(rbind,lapply(split(OrigDroplets.df,OrigDroplets.df$Check),function(g){
    data.frame(Check = g$Check[1],
               nObjects = nrow(g),
               pct.SingleFrame = 100*mean(g$nFrames==1),
               median.TrackFrames = median(g$nFrames),
               median.d.mm = median(g$Diameter.mm),
               median.dy.mm.per.frame = median(g$dy.mm.per.frame,na.rm=T))
  }))
  cat('\n---- Manual check: deleted vs kept objects ----\n')
  print(Motion.Summary,digits=3)
  del = OrigDroplets.df[OrigDroplets.df$Check=='deleted',]
  kep = OrigDroplets.df[OrigDroplets.df$Check=='kept',]
  if (sum(!is.na(del$dy.mm.per.frame))>=5 & sum(!is.na(kep$dy.mm.per.frame))>=5){
    cat('Wilcoxon, vertical travel per frame (deleted vs kept): p =',
        signif(wilcox.test(del$dy.mm.per.frame,kep$dy.mm.per.frame,exact=F)$p.value,3),'\n')
  }
  cat('Fisher test, share of single-frame objects (deleted vs kept): p =',
      signif(fisher.test(table(OrigDroplets.df$Check,OrigDroplets.df$nFrames==1))$p.value,3),'\n')
  write.csv(Motion.Summary,outFile('ManualCheck_Motion.csv'),row.names=F)
  write.csv(OrigDroplets.df,outFile('OriginalDroplets_Status.csv'),row.names=F)
}

################################################################################
## 3. Stack locations (filled in by hand once)
################################################################################
# First run: writes a template with the stack IDs split at '_' and stops.
# Fill in Depth (front/middle/back), Zone (center/middle/edge), Side, Height and
# Uncertain (TRUE if the zone is a guess), save, and run again.

if (!file.exists(locationFile)){
  parts = strsplit(stackFrames$StackID,'_')
  template = data.frame(StackID = stackFrames$StackID,
                        Location.Label = sapply(parts,function(p) ifelse(length(p)>=2,p[2],'')),
                        Depth = '',
                        Zone = '',
                        Side = '',
                        Height = sapply(parts,function(p) ifelse(length(p)>=4,p[4],'top')), #no height in the name = top window
                        Uncertain = FALSE)
  #guess depth from the label where it is obvious; check every row
  template$Depth[grepl('front',template$Location.Label,ignore.case=T)] = 'front'
  template$Depth[grepl('back',template$Location.Label,ignore.case=T)] = 'back'
  write.csv(template,locationFile,row.names=F)
  stop(paste('Location template written to',locationFile,
             '\nFill in Depth, Zone, Side, Height and Uncertain, then run the script again.'))
}

Locations.df = read.csv(locationFile,stringsAsFactors=F)
Locations.df[is.na(Locations.df)] = ''
notListed = setdiff(stackFrames$StackID,Locations.df$StackID)
if (length(notListed)>0) warning(paste('Stacks missing from',basename(locationFile),':',paste(notListed,collapse=', ')))
noZone = Locations.df$StackID[Locations.df$Zone=='' & Locations.df$StackID %in% stackFrames$StackID]
if (length(noZone)>0) warning(paste(length(noZone),'stack(s) have no Zone assigned; they are left out of zone comparisons'))

locCols = c('StackID','Location.Label','Depth','Zone','Side','Height','Uncertain')
Droplets.df = merge(Droplets.df,Locations.df[,locCols],by='StackID',all.x=T)
Frames.df = merge(Frames.df,Locations.df[,c('StackID','Zone','Depth','Height')],by='StackID',all.x=T)
for (col in locCols[-1]) Droplets.df[[col]][is.na(Droplets.df[[col]])] = ''

################################################################################
## 4. Per-stack summary (Section 3C Step 1b; R1-6; R1-7)
################################################################################

byStack = split(Droplets.df,Droplets.df$StackID)
Stack.Summary = do.call(rbind,lapply(byStack,function(s){
  id = s$StackID[1]
  nTiff = stackFrames$nFramesStack[stackFrames$StackID==id]
  data.frame(StackID = id,
             Location.Label = s$Location.Label[1],
             Depth = s$Depth[1],
             Zone = s$Zone[1],
             Side = s$Side[1],
             Height = s$Height[1],
             Uncertain = s$Uncertain[1],
             n = nrow(s),
             nFramesStack = nTiff, #tiff files in the stack folder
             px.per.mm = median(Frames.df$xCentroid.px[Frames.df$StackID==id]/Frames.df$xCentroid.mm[Frames.df$StackID==id]), #calibration scale used for this stack
             FramesWithDrops = length(unique(Frames.df$FrameNo[Frames.df$StackID==id])),
             DropsPer100Frames = 100*nrow(s)/nTiff,
             MedianTrackFrames = median(s$nFrames),
             mean.d = mean(s$Diameter.mm),
             sd.d = sd(s$Diameter.mm),
             cv = sd(s$Diameter.mm)/mean(s$Diameter.mm),
             d32 = sum(s$Diameter.mm^3)/sum(s$Diameter.mm^2), #Sauter mean diameter
             min.d = min(s$Diameter.mm),
             max.d = max(s$Diameter.mm),
             mean.aspect = mean(s$AspectRatio),
             drip.num = mean(s$Dripping),
             drip.vol = sum(s$Volume.mm3[s$Dripping])/sum(s$Volume.mm3))
}))
rownames(Stack.Summary) = NULL
if (length(emptyStacks)>0){
  emptyRows = Stack.Summary[rep(1,length(emptyStacks)),]
  emptyRows[,] = NA
  emptyRows$StackID = emptyStacks
  loc = Locations.df[match(emptyStacks,Locations.df$StackID),]
  for (col in c('Location.Label','Depth','Zone','Side','Height','Uncertain')) emptyRows[[col]] = ifelse(is.na(loc[[col]]),'',loc[[col]])
  emptyRows$n = 0
  emptyRows$FramesWithDrops = 0
  emptyRows$nFramesStack = stackFrames$nFramesStack[match(emptyStacks,stackFrames$StackID)]
  emptyRows$DropsPer100Frames = 0
  Stack.Summary = rbind(Stack.Summary,emptyRows)
  rownames(Stack.Summary) = NULL
}
if (!is.null(Deletion.Summary)){
  Stack.Summary = merge(Stack.Summary,Deletion.Summary[,c('StackID','OrigDroplets','DropletsDeleted','pct.DropletsDeleted',
                                                          'mean.d.deleted','drip.num.deleted')],by='StackID',all.x=T)
}
Stack.Summary = Stack.Summary[order(Stack.Summary$Zone,Stack.Summary$StackID),]

cat('\n---- Per-stack summary ----\n')
print(Stack.Summary[,c('StackID','Zone','Depth','Height','n','nFramesStack','mean.d','cv','drip.num','drip.vol')],digits=3)

################################################################################
## 5. Pooled summary for this data set
################################################################################

Pooled.Summary = data.frame(Data = dataName,
                            FlowRate.LPM = FlowRate.LPM,
                            nStacks = nrow(Stack.Summary),
                            n = nrow(Droplets.df),
                            mean.d = mean(Droplets.df$Diameter.mm),
                            sd.d = sd(Droplets.df$Diameter.mm),
                            cv = sd(Droplets.df$Diameter.mm)/mean(Droplets.df$Diameter.mm),
                            d32 = sum(Droplets.df$Diameter.mm^3)/sum(Droplets.df$Diameter.mm^2),
                            min.d = min(Droplets.df$Diameter.mm),
                            max.d = max(Droplets.df$Diameter.mm),
                            drip.num = mean(Droplets.df$Dripping),
                            drip.vol = sum(Droplets.df$Volume.mm3[Droplets.df$Dripping])/sum(Droplets.df$Volume.mm3),
                            #spread of stack means: orifice-to-orifice / location variation (R1-7)
                            sd.of.stack.means = sd(Stack.Summary$mean.d,na.rm=T))
cat('\n---- Pooled ----\n')
print(Pooled.Summary,digits=3)

################################################################################
## 6. Dripping volume fraction four ways (Section 3C Step 1b)
################################################################################
# (a) pooled as the model uses it, (b) equal weight per stack, (c) center stacks
# only, (d) edge stacks only. Compare each with the outlet-concentration bound
# (~10% of liquid by volume, Step 4.4).

pooledVol = function(s) ifelse(nrow(s)==0,NA,sum(s$Volume.mm3[s$Dripping])/sum(s$Volume.mm3))
pooledNum = function(s) ifelse(nrow(s)==0,NA,mean(s$Dripping))
ctr = Droplets.df[Droplets.df$Zone=='center',]
edg = Droplets.df[Droplets.df$Zone=='edge',]

Drip.FourWays = data.frame(Method = c('a. pooled','b. equal weight per stack','c. center stacks only','d. edge stacks only'),
                           nStacks = c(sum(Stack.Summary$n>0),sum(Stack.Summary$n>0),length(unique(ctr$StackID)),length(unique(edg$StackID))),
                           drip.vol = c(pooledVol(Droplets.df),mean(Stack.Summary$drip.vol,na.rm=T),pooledVol(ctr),pooledVol(edg)),
                           drip.num = c(pooledNum(Droplets.df),mean(Stack.Summary$drip.num,na.rm=T),pooledNum(ctr),pooledNum(edg)))
cat('\n---- Dripping fraction (cutoff',drip.cutoff.mm,'mm) ----\n')
print(Drip.FourWays,digits=3)

################################################################################
## 7. Spatial patterns: zone, depth and height (stack is the unit of analysis)
################################################################################
# Droplets in one stack share an orifice/location, so comparisons use stack
# values, not individual droplets.

groupSummary = function(st,group){
  st = st[st[[group]]!='',]
  if (nrow(st)==0) return(NULL)
  out = do.call(rbind,lapply(split(st,st[[group]]),function(g){
    data.frame(Group = g[[group]][1],
               nStacks = nrow(g),
               nDroplets = sum(g$n),
               mean.d = mean(g$mean.d,na.rm=T),
               se.d = ifelse(sum(!is.na(g$mean.d))>1,sd(g$mean.d,na.rm=T)/sqrt(sum(!is.na(g$mean.d))),NA),
               mean.cv = mean(g$cv,na.rm=T),
               drip.vol = mean(g$drip.vol,na.rm=T),
               DropsPer100Frames = mean(g$DropsPer100Frames,na.rm=T),
               pct.DropletsDeleted = ifelse(is.null(g$pct.DropletsDeleted),NA,mean(g$pct.DropletsDeleted,na.rm=T)))
  }))
  names(out)[1] = group
  rownames(out) = NULL
  out
}

kwTest = function(st,group){
  #Kruskal-Wallis on stack values
  st = st[st[[group]]!='',]
  if (length(unique(st[[group]]))<2) return(NULL)
  data.frame(Grouping = group,
             p.mean.d = kruskal.test(st$mean.d,factor(st[[group]]))$p.value,
             p.cv = kruskal.test(st$cv,factor(st[[group]]))$p.value,
             p.drip.vol = kruskal.test(st$drip.vol,factor(st[[group]]))$p.value)
}

Zone.Summary = groupSummary(Stack.Summary,'Zone')
Depth.Summary = groupSummary(Stack.Summary,'Depth')
Height.Summary = groupSummary(Stack.Summary,'Height')
Spatial.Tests = rbind(kwTest(Stack.Summary,'Zone'),kwTest(Stack.Summary,'Depth'),kwTest(Stack.Summary,'Height'))

cat('\n---- By zone (center/middle/edge) ----\n'); print(Zone.Summary,digits=3)
cat('\n---- By depth (front/middle/back; R1-6) ----\n'); print(Depth.Summary,digits=3)
cat('\n---- By height ----\n'); print(Height.Summary,digits=3)
cat('\n---- Kruskal-Wallis on stack values ----\n'); print(Spatial.Tests,digits=3)

#zone sensitivity: stacks marked Uncertain left out
if (any(Stack.Summary$Uncertain %in% c(TRUE,'TRUE'))){
  cat('\n---- By zone, uncertain stacks left out ----\n')
  print(groupSummary(Stack.Summary[!(Stack.Summary$Uncertain %in% c(TRUE,'TRUE')),],'Zone'),digits=3)
}

#paired heights: same location, two heights (same orifices feed both)
Height.Pairs = NULL
hs = Stack.Summary[Stack.Summary$Height!='' & Stack.Summary$Location.Label!='',]
heights = sort(unique(hs$Height))
if (length(heights)==2){
  keep = c('Location.Label','n','mean.d','cv','drip.vol')
  Height.Pairs = merge(hs[hs$Height==heights[1],keep],hs[hs$Height==heights[2],keep],
                       by='Location.Label',suffixes=paste0('.',heights))
  Height.Pairs$diff.mean.d = Height.Pairs[[paste0('mean.d.',heights[2])]]-Height.Pairs[[paste0('mean.d.',heights[1])]]
  cat('\n---- Paired heights ----\n'); print(Height.Pairs,digits=3)
  if (nrow(Height.Pairs)>=5){
    w = wilcox.test(Height.Pairs[[paste0('mean.d.',heights[1])]],Height.Pairs[[paste0('mean.d.',heights[2])]],paired=T,exact=F)
    cat('Paired Wilcoxon on mean diameter, p =',signif(w$p.value,3),'\n')
  }
} else if (length(heights)>2){
  cat('\nMore than two heights found (',paste(heights,collapse=', '),'); paired-height table skipped.\n')
}

################################################################################
## 8. Figures
################################################################################

breaks = seq(floor(min(Droplets.df$Diameter.mm)/binwidth.mm)*binwidth.mm,
             ceiling(max(Droplets.df$Diameter.mm)/binwidth.mm)*binwidth.mm+binwidth.mm,by=binwidth.mm)
plotTitle = paste0(dataName,', ',FlowRate.LPM,' L/min (n = ',nrow(Droplets.df),')')

#pooled distribution, number- and volume-weighted
SizeDist.plot = ggplot(Droplets.df,aes(x=Diameter.mm))+
  geom_histogram(aes(y=after_stat(density)),breaks=breaks,fill='grey60',color='black')+
  geom_vline(xintercept=drip.cutoff.mm,linetype='dashed')+
  labs(x='Droplet diameter (mm)',y='Probability density',title=plotTitle)+
  theme_bw()
ggsave(outFile('SizeDistribution.png'),SizeDist.plot,width=7,height=5,dpi=300)

VolDist.plot = ggplot(Droplets.df,aes(x=Diameter.mm,weight=Volume.mm3))+
  geom_histogram(aes(y=after_stat(density)),breaks=breaks,fill='grey60',color='black')+
  geom_vline(xintercept=drip.cutoff.mm,linetype='dashed')+
  labs(x='Droplet diameter (mm)',y='Volume-weighted probability density',title=plotTitle)+
  theme_bw()
ggsave(outFile('SizeDistribution_VolumeWeighted.png'),VolDist.plot,width=7,height=5,dpi=300)

#per-stack distributions (Supplement; R1-7, R2-9)
Droplets.df$Panel = paste0(Droplets.df$StackID,'\n',Droplets.df$Zone,' ',Droplets.df$Depth,
                           ' (n = ',ave(Droplets.df$Diameter.mm,Droplets.df$StackID,FUN=length),')')
PerStack.plot = ggplot(Droplets.df,aes(x=Diameter.mm))+
  geom_histogram(aes(y=after_stat(density)),breaks=breaks,fill='grey60',color='black')+
  geom_vline(xintercept=drip.cutoff.mm,linetype='dashed')+
  facet_wrap(~Panel,ncol=4)+
  labs(x='Droplet diameter (mm)',y='Probability density',title=paste(dataName,',',FlowRate.LPM,'L/min'))+
  theme_bw()+theme(strip.text=element_text(size=6))
ggsave(outFile('SizeDistribution_PerStack.png'),PerStack.plot,width=10,
       height=2+2*ceiling(nrow(Stack.Summary)/4),dpi=300,limitsize=F)

#diameters of droplets kept vs deleted in the manual check (counts, original MATLAB sizes)
if (!is.null(OrigDroplets.df)){
  od = OrigDroplets.df
  od$Check = ifelse(od$Status=='deleted','Deleted','Kept (incl. trimmed)')
  obreaks = seq(floor(min(od$Diameter.mm)/binwidth.mm)*binwidth.mm,
                ceiling(max(od$Diameter.mm)/binwidth.mm)*binwidth.mm+binwidth.mm,by=binwidth.mm)
  Deleted.plot = ggplot(od,aes(x=Diameter.mm,fill=Check))+
    geom_histogram(breaks=obreaks,color='black')+
    geom_vline(xintercept=drip.cutoff.mm,linetype='dashed')+
    scale_fill_manual(values=c('Deleted'='#CC6677','Kept (incl. trimmed)'='grey70'))+
    labs(x='Droplet diameter (mm, original MATLAB measurement)',y='Number of droplets',fill='',
         title=paste(dataName,',',FlowRate.LPM,'L/min: droplets removed in manual check'))+
    theme_bw()+theme(legend.position=c(0.98,0.98),legend.justification=c(1,1))
  ggsave(outFile('ManualCheck_KeptVsDeleted.png'),Deleted.plot,width=7,height=5,dpi=300)
}

#stack mean diameter and dripping fraction by zone (one point per stack)
st = Stack.Summary[Stack.Summary$Zone!='' & Stack.Summary$n>0,]
if (nrow(st)>0){
  st$Zone = factor(st$Zone,levels=c(intersect(c('center','middle','edge'),unique(st$Zone)),
                                    setdiff(unique(st$Zone),c('center','middle','edge'))))
  Zone.plot = ggplot(st,aes(x=Zone,y=mean.d,shape=Depth))+
    geom_point(size=3,position=position_jitter(width=0.1,height=0,seed=1))+
    labs(x='Zone',y='Stack mean diameter (mm)',title=paste(dataName,',',FlowRate.LPM,'L/min'))+
    theme_bw()
  ggsave(outFile('StackMeans_ByZone.png'),Zone.plot,width=6,height=4.5,dpi=300)

  Drip.plot = ggplot(st,aes(x=Zone,y=100*drip.vol,shape=Depth))+
    geom_point(size=3,position=position_jitter(width=0.1,height=0,seed=1))+
    geom_hline(yintercept=10,linetype='dashed')+ #outlet-concentration bound (Step 4.4)
    labs(x='Zone',y='Dripping fraction (% of liquid volume)',title=paste(dataName,',',FlowRate.LPM,'L/min'))+
    theme_bw()
  ggsave(outFile('DrippingFraction_ByZone.png'),Drip.plot,width=6,height=4.5,dpi=300)
}

################################################################################
## 9. Exports
################################################################################

write.csv(Droplets.df[,setdiff(names(Droplets.df),'Panel')],outFile('Droplets.csv'),row.names=F)
write.csv(Stack.Summary,outFile('Stack_Summary.csv'),row.names=F)
write.csv(Pooled.Summary,outFile('Pooled_Summary.csv'),row.names=F)
write.csv(Drip.FourWays,outFile('DrippingFraction_FourWays.csv'),row.names=F)
if (!is.null(Zone.Summary)) write.csv(Zone.Summary,outFile('Zone_Summary.csv'),row.names=F)
if (!is.null(Depth.Summary)) write.csv(Depth.Summary,outFile('Depth_Summary.csv'),row.names=F)
if (!is.null(Height.Summary)) write.csv(Height.Summary,outFile('Height_Summary.csv'),row.names=F)
if (!is.null(Spatial.Tests)) write.csv(Spatial.Tests,outFile('Spatial_Tests.csv'),row.names=F)
if (!is.null(Height.Pairs)) write.csv(Height.Pairs,outFile('Height_Pairs.csv'),row.names=F)

#frame-level data with centroids, for the velocity/acceleration analysis later
write.csv(Frames.df,outFile('Frames.csv'),row.names=F)

#AllDropletData, same columns as dataCleaning.R so the MATLAB scripts read it unchanged
compositeData = data.frame(Droplet.Designation = Frames.df$Droplet.Designation,
                           FrameNo = Frames.df$FrameNo,
                           xDiameter.mm = Frames.df$xDiameter.mm,
                           yDiameter.mm = Frames.df$yDiameter.mm)
compositeData = compositeData[order(compositeData$Droplet.Designation,compositeData$FrameNo),]
compositeData$Droplet.No = match(compositeData$Droplet.Designation,unique(compositeData$Droplet.Designation))
rownames(compositeData) = NULL
write.csv(compositeData,outFile('AllDropletData.csv'))

cat('\nOutputs written to',outFolder,'\n')
