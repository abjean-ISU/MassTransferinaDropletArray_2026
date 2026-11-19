#SupportingCalculations.R
#Short calculations behind the reviewer answers that do not need the model or the images.
#Each section names the comment it supports. Run top to bottom; everything prints to the console
#and the main numbers are collected in SupportingCalculations.csv.

rm(list=ls())

results = data.frame(Section=character(),Quantity=character(),Value=character(),stringsAsFactors=F)
report = function(section,quantity,value){
  results[nrow(results)+1,] <<- c(section,quantity,value)
  cat(sprintf('  %-60s %s\n',quantity,value))
}

#common inputs
g = 9.8                 #m/s2
rho.l = 997             #water density (kg/m3)
mu.l = 0.89e-3          #water viscosity at 25 C (Pa s)
d.o = 0.8e-3            #orifice (nub bore) diameter (m), re-measured 27 Oct 2026
N0 = 200                #number of orifices
Q.LPM = c(3.0,4.5)      #flow rates (L/min)
sigma = function(T.C){T = T.C+273.15; 235.8e-3*(1-T/647.15)^1.256*(1-0.625*(1-T/647.15))}  #surface tension (N/m), as in the model
pixel.m = 10e-6         #Chronos 2.1-HD pixel size
f.lens = 0.100          #lens focal length (m)
scale.range = c(24,30,39)   #px/mm, smallest/typical/largest calibration across all stacks

################################################################################
cat('\n==== 2. Shape filters (R1-5) ====\n')
################################################################################
#MATLAB regionprops Circularity = 4*pi*Area/Perimeter^2; ellipse perimeter from Ramanujan's approximation.
#Find the aspect ratio at which the 0.8 lower limit is reached.
circ.ellipse = function(AR){a = AR; b = 1; P = pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b))); 4*pi*(pi*a*b)/P^2}
AR.limit = uniroot(function(AR) circ.ellipse(AR)-0.8,c(1,5))$root
report('2 Filters','Aspect ratio at the circularity limit 0.8',sprintf('%.2f',AR.limit))

#n = 2 oscillation about the vertical: r = R[1 + eps P2(cos theta)], so aspect ratio = (1+eps)/(1-eps/2)
#for a tall drop and (1+eps/2)/(1-eps) for a flat drop. For small eps, |ln AR| = 1.5 eps.
eps.flat = uniroot(function(e) (1+e/2)/(1-e)-AR.limit,c(0,0.99))$root
eps.tall = uniroot(function(e) (1+e)/(1-e/2)-AR.limit,c(0,1.9))$root
report('2 Filters','n = 2 amplitude eps that reaches the circularity limit: flat / tall',
       sprintf('%.2f / %.2f',eps.flat,eps.tall))

#measured values from AspectRatioAnalysis.R, 4.5 L/min, top window:
#within-track SD of log(AR), noise-corrected; sinusoid amplitude = sqrt(2) x SD; eps = amplitude/1.5
sd.top = 0.0320; noise.top = 0.0088
eps.top = sqrt(2)*sqrt(sd.top^2-noise.top^2)/1.5
decay.top = 0.422                                #exp(-t/tau_damp) at the top window, AspectRatioAnalysis.R
eps.max = log(1.15)/1.5                          #largest single-drop swing, about +/-15% in AR
report('2 Filters','Measured eps: median at top window / back-extrapolated to breakup / largest drop',
       sprintf('%.3f / %.3f / %.3f',eps.top,eps.top/decay.top,eps.max))
report('2 Filters','Margin, limit / measured: top window; breakup; largest drop',
       sprintf('%.0f-%.0fx; %.0f-%.0fx; %.0f-%.0fx',eps.flat/eps.top,eps.tall/eps.top,
               eps.flat/(eps.top/decay.top),eps.tall/(eps.top/decay.top),eps.flat/eps.max,eps.tall/eps.max))

################################################################################
cat('\n==== 3. Calibration precision (R1-6) ====\n')
################################################################################
#ASSUMPTION: each end of the 1 mm calibration line placed to +/-1 px (rectangular): u(length) = sqrt(2)/sqrt(3) px
u.len.px = sqrt(2)/sqrt(3)
report('3 Calibration','Scale standard uncertainty per stack at 24 / 30 / 39 px/mm (%)',
       paste(sprintf('%.1f',100*u.len.px/scale.range),collapse=' / '))
report('3 Calibration','Effect on the pooled mean over 24 stacks (typical scale) (%)',sprintf('%.1f',100*u.len.px/30/sqrt(24)))

################################################################################
cat('\n==== 4. Depth of field and sizing error (R1-6) ====\n')
################################################################################
#m = P x pixel size; focus distance u = f(1+1/m); DOF = 2 N c (1+m)/m^2; max size error = N c/(m f)
#Acceptable blur c taken as 1-2 pixels (assumption).
N = 2.8
for (P in scale.range){
  m = P*1000*pixel.m
  u = f.lens*(1+1/m)
  c.blur = c(1,2)*pixel.m
  DOF = 2*N*c.blur*(1+m)/m^2
  err = N*c.blur/(m*f.lens)
  report('4 Depth of field',sprintf('f/%g, %d px/mm: focus distance (cm); DOF (mm); max size error (%%)',N,P),
         sprintf('%.0f; %.1f-%.1f; %.2f-%.2f',100*u,1000*DOF[1],1000*DOF[2],100*err[1],100*err[2]))
  #conservative case: distance from the focal plane and edge blur that give a 1% size error
  report('4 Depth of field',sprintf('%d px/mm: for a 1%% size error, distance from focal plane (mm); blur (px)',P),
         sprintf('%.1f; %.1f',1000*0.01*u,0.01*m*f.lens/N/pixel.m))
}

################################################################################
cat('\n==== 6. Residence time bounds and oscillation (R2-4, R2-7) ====\n')
################################################################################
for (H in c(0.53,0.68,0.98)) report('6 Free fall',sprintf('Free fall from rest over %.2f m (s)',H),sprintf('%.3f',sqrt(2*H/g)))
r.typ = 0.8e-3
report('6 Oscillation','Lamb damping time rho r^2/(5 mu), r = 0.8 mm (s)',sprintf('%.3f',rho.l*r.typ^2/(5*mu.l)))
report('6 Oscillation','Rayleigh n = 2 frequency, r = 0.8 mm (Hz)',sprintf('%.0f',sqrt(8*sigma(22)/(rho.l*r.typ^3))/(2*pi)))

#vessel conditions and model properties (same equations as the MATLAB model, 27 Oct rerun)
T.vessel = 22                                    #typical trial temperature (C); trials span 18-25 C
P.vessel = (14.2+5)/14.7                         #atm
Tk = T.vessel+273.15
rho.CO2 = P.vessel*101325*0.04401/(8.314462*Tk)  #CO2 density, ideal gas (kg/m3)
mu.CO2 = 1.370e-5*(Tk/273.15)^1.5*(273.15+240)/(Tk+240)  #CO2 viscosity, Sutherland (Pa s)
rho.air = 1.184; mu.air = 1.86e-5                #imaging rig (room air)
D.AB = 13.942e-9*(Tk/227-1)^1.7094               #CO2 diffusivity in water (m2/s)
Csat = P.vessel/(1/(exp(-159.854+8741.68/Tk+21.6694*log(Tk)-1.10261e-3*Tk))/55342)*44.01/1000 #g/L
#model quadrature nodes (RunTrial, 28 Oct): radius (m) and number weights
r.nodes = c(0.5568,0.7515,0.9521,1.2635)*1e-3; w.nodes = c(0.0794,0.5745,0.3381,0.0080)
V.mean = sum(w.nodes*4/3*pi*r.nodes^3)/sum(w.nodes)   #Eq. 24 mean drop volume (m3)
d.drop = (6*V.mean/pi)^(1/3)                          #diameter of the mean volume (m), 1.65 mm
report('7 Geometry','Eq. 24 mean drop volume (mm3); diameter of mean volume (mm); / orifice diameter',
       sprintf('%.3f; %.3f; %.2f',V.mean*1e9,d.drop*1000,d.drop/d.o))

#drag law of the model and a simple trajectory (speed and time at distance x)
accel = function(v,d,rho.g,mu.g){r = d/2; g-(9/2)*(mu.g/(r^2*rho.l))*v*(1+0.158*(2*rho.g*v*r/mu.g)^(2/3))}
traj = function(d,v0,x.end,rho.g,mu.g,dx=1e-3){x = 0; v = v0; t = 0
while (x < x.end-1e-12){t = t+dx/v; v = v+accel(v,d,rho.g,mu.g)/v*dx; x = x+dx}; c(v=v,t=t)}
vTerm = function(d,rho.g,mu.g) uniroot(function(v) accel(v,d,rho.g,mu.g),c(1e-6,50))$root
v0of = function(Q) Q/1000/60/(N0*pi*d.o^2/4)

################################################################################
cat('\n==== 7. Manifold geometry and droplet spacing (R2-3, R1-8) ====\n')
################################################################################
#manifold: 9 in showerhead, 200 orifices; field radii and pitch ESTIMATED from a photo (+/-15%)
D.manifold = 9*0.0254                            #m
r.field = c(0.016,0.100)                         #inner and outer radius of the orifice field (m), photo estimate
pitch = c(min=0.006,mean=0.013,max=0.030)        #nearest-neighbour spacing (m), photo estimate
A.field = pi*(r.field[2]^2-r.field[1]^2)
report('7 Geometry','Orifice field area (m2) [photo estimate]',sprintf('%.3f',A.field))
report('7 Geometry','Orifice pitch min / mean / max (mm) [photo estimate]',paste(1000*pitch,collapse=' / '))
Oh = mu.l/sqrt(rho.l*sigma(T.vessel)*d.o)        #Ohnesorge number, used in Section 8
report('7 Regime','Ohnesorge number at the orifice',sprintf('%.4f',Oh))
for (Q in Q.LPM){
  v0 = v0of(Q)
  q = Q/1000/60/N0                                #flow per orifice (m3/s)
  f.drop = q/V.mean                               #drops per second per orifice, from the Eq. 24 mean volume
  report('7 Regime',sprintf('%.1f L/min: v0 (m/s); drops per second per orifice',Q),sprintf('%.3f; %.0f',v0,f.drop))
  #axial spacing between successive drops from one orifice = v/f, at the orifice and at 0.15, 0.375, 0.98 m (CO2)
  xs = c(0,0.15,0.375,0.98)
  vs = sapply(xs,function(x) if (x==0) v0 else traj(d.drop,v0,x,rho.CO2,mu.CO2)['v'])
  report('7 Spacing',sprintf('%.1f L/min: speed at 0 / 0.15 / 0.375 / 0.98 m (m/s, CO2)',Q),paste(sprintf('%.2f',vs),collapse=' / '))
  report('7 Spacing',sprintf('%.1f L/min: axial spacing in drop diameters at the same points',Q),paste(sprintf('%.1f',vs/f.drop/d.drop),collapse=' / '))
  #gas Weber number at the fastest point in the column
  report('7 Regime',sprintf('%.1f L/min: gas Weber number at 0.98 m (CO2)',Q),
         sprintf('%.2f',rho.CO2*vs[4]^2*d.drop/sigma(T.vessel)))
}

################################################################################
cat('\n==== 8. Jet breakup length (R2-2) ====\n')
################################################################################
#Kooij et al. (2025), Phys. Rev. Lett. 135, 214001:
#Eq. 4: L/D = ln(D/(2*delta0)) * We^0.5 * (1+3*Z), with We = rho*D*V^2/sigma and Z = Ohnesorge number
#Eq. 6: delta0 = sqrt(kB*T/(sqrt(8)*pi^2*sigma)), the thermal capillary-wave amplitude
#Constant jet speed; gravity is not included. Kooij et al. do not address the regime near the
#dripping-jetting transition, so the 3.0 L/min value is less certain. No images of the jets exist.
kB = 1.380649e-23                                                   #Boltzmann constant (J/K)
delta0 = sqrt(kB*Tk/(sqrt(8)*pi^2*sigma(T.vessel)))                 #thermal amplitude (m)
report('8 Breakup','Thermal amplitude delta0 (m)',sprintf('%.2e',delta0))
for (Q in Q.LPM){
  v0 = v0of(Q)
  We = rho.l*v0^2*d.o/sigma(T.vessel)
  L = d.o*log(d.o/(2*delta0))*sqrt(We)*(1+3*Oh)  #breakup length (m)
  tb = L/v0                                      #time on the intact jet (s)
  report('8 Breakup',sprintf('%.1f L/min: breakup length (mm); time on the jet (ms)',Q),sprintf('%.0f; %.0f',1000*L,1000*tb))
}

################################################################################
cat('\n==== 9. Exit velocity vs terminal velocity (R1-3) ====\n')
################################################################################
#Terminal velocity from the model drag law (this is how U_inf is obtained for Amokrane: the speed at which
#the model acceleration is zero, at the trial's gas properties). Compare with Gunn & Kinzer (1949) for air:
#fill gk.U from their Table 2 to check the drag law (values not entered here).
gk = data.frame(d.mm = c(1.0,1.6,2.0,3.0), U = NA)
for (dd in c(1.1,1.5,1.6,1.9,2.5)){
  Ua = vTerm(dd/1000,rho.air,mu.air); Uc = vTerm(dd/1000,rho.CO2,mu.CO2)
  report('9 Terminal',sprintf('d = %.1f mm: terminal velocity air / CO2 (m/s)',dd),sprintf('%.2f / %.2f',Ua,Uc))
}
for (Q in Q.LPM){
  v0 = v0of(Q); Uc = vTerm(d.drop,rho.CO2,mu.CO2)
  v98 = traj(d.drop,v0,0.98,rho.CO2,mu.CO2)['v']
  report('9 Terminal',sprintf('%.1f L/min, 1.6 mm, CO2: v0 / U_inf (%%); speed at 0.98 m / U_inf (%%)',Q),
         sprintf('%.0f; %.0f',100*v0/Uc,100*v98/Uc))
}
if (all(!is.na(gk$U))) for (i in 1:nrow(gk)) report('9 Terminal',sprintf('Gunn & Kinzer d = %.1f mm: measured vs model (m/s)',gk$d.mm[i]),
                                                    sprintf('%.2f vs %.2f',gk$U[i],vTerm(gk$d.mm[i]/1000,rho.air,mu.air)))

################################################################################
cat('\n==== 10. Holdup gap and drainage (R2-4) ====\n')
################################################################################
#Table 3 values (ManuscriptTables.xlsx, 27 Oct rerun): measured and model holdup (L)
LH = data.frame(Q = rep(c(3.0,4.5),each=3), H = rep(c(0.53,0.68,0.98),2),
                VL.meas = c(0.034,0.047,0.050,0.044,0.043,0.074),
                VL.model = c(0.015,0.017,0.022,0.021,0.025,0.031))
LH$gap = LH$VL.meas-LH$VL.model
LH$gap.flow.s = LH$gap/(LH$Q/60)                 #seconds of full flow that the gap represents
LH$gap.pct.trapped.low = 100*LH$gap/0.6           #% of the 0.2-0.6 L trapped between valve and manifold
LH$gap.pct.trapped.high = 100*LH$gap/0.2
print(LH,row.names=F)
report('10 Drainage','Holdup gap, measured - model (L)',sprintf('%.3f to %.3f',min(LH$gap),max(LH$gap)))
report('10 Drainage','Gap as seconds of full flow after closure',sprintf('%.2f to %.2f',min(LH$gap.flow.s),max(LH$gap.flow.s)))
report('10 Drainage','Gap as % of the liquid trapped upstream (0.2-0.6 L)',
       sprintf('%.0f to %.0f',min(LH$gap.pct.trapped.low),max(LH$gap.pct.trapped.high)))
report('10 Drainage','Floor-scale resolution 0.02 kg: SD per trial (kg) and as % of mean holdup',
       sprintf('%.4f; %.0f to %.0f',0.02/sqrt(12),100*0.02/sqrt(12)/max(LH$VL.meas),100*0.02/sqrt(12)/min(LH$VL.meas)))

################################################################################
cat('\n==== 11. kLa per unit column volume (R2-1) ====\n')
################################################################################
#AUTHOR TO CHECK: proxy column = manifold diameter + 0.1 m (no drop-wall contact), height 0.98 m
D.col = D.manifold+0.1; H.col = 0.98
V.col = pi/4*D.col^2*H.col*1000                  #L
report('11 Column','Proxy column diameter (m); height (m); volume (L) [CHECK]',sprintf('%.4f; %.2f; %.1f',D.col,H.col,V.col))
for (Q in Q.LPM) report('11 Column',sprintf('%.1f L/min: liquid loading (m3/m2/h)',Q),sprintf('%.2f',Q/1000*60/(pi/4*D.col^2)))
#Table 4 kLa (s-1, liquid-volume basis) and Table 3 holdup; kLa_col = kLa * V_L / V_col = (Q/V_col) ln(Csat/(Csat-C))
kLa.tab = data.frame(Q = rep(c(3.0,4.5),each=3), H = rep(c(0.53,0.68,0.98),2),
                     kLa.model = c(2.12,2.00,2.00,2.15,2.01,2.08), kLa.exp = c(1.06,0.82,0.84,1.38,1.54,1.02))
kLa.tab = merge(kLa.tab,LH[,c('Q','H','VL.meas','VL.model')])
kLa.tab$holdup.frac.model = kLa.tab$VL.model/V.col
kLa.tab$kLa.col.model = kLa.tab$kLa.model*kLa.tab$VL.model/V.col
kLa.tab$kLa.col.exp = kLa.tab$kLa.exp*kLa.tab$VL.meas/V.col
print(kLa.tab,row.names=F,digits=3)
report('11 Column','Holdup fraction (model V_L / V_col)',sprintf('%.1e to %.1e',min(kLa.tab$holdup.frac.model),max(kLa.tab$holdup.frac.model)))
report('11 Column','kLa per column volume, model (s-1)',sprintf('%.1e to %.1e',min(kLa.tab$kLa.col.model),max(kLa.tab$kLa.col.model)))
report('11 Column','kLa per column volume, experimental (s-1)',sprintf('%.1e to %.1e',min(kLa.tab$kLa.col.exp),max(kLa.tab$kLa.col.exp)))
report('11 Column','Note','the holdup cancels on this basis: kLa_col depends only on Q, V_col and the outlet C')
write.csv(results,'SupportingCalculations.csv',row.names=F)
cat('\nWritten: SupportingCalculations.csv\n')
