%runSensitivityBatch.m
%One MATLAB batch for the remaining model runs (28 Oct 2026). Self-contained: needs only this file,
%PD_gqmom_0_1.m, the experimental Excel file and the two droplet data files. Run time: a few minutes.
%
%  A. Reproduction check   Angelo, 4.5 L/min, 0.98 m, trial 1 must give 1.1707 g/L, tau 0.417 s (Table 1)
%  B. Hybrid runs (R2-7)   below Re = 200 use a circulating-drop kL, above it Angelo:
%                            Amokrane (realistic circulation, low end),
%                            Ruckenstein Eq. 49 (fully mixed potential-flow circulation, upper bound)
%                          + where each model drop size crosses Re = 200 (height, time)
%  C. Size sensitivity     all radii x0.8/0.9/1.1/1.2; distribution width x0.5/1.5 (R2-9, Section 3C)
%                          + the radius scale factor that would match the measured mean at each flow rate
%  D. Dripping (3C/R2-3)   3.0 L/min trials with 10/20/30% of the volume as 2.5/2.75/3.1 mm dripping drops
%                          leaving the orifice at 0.05 or 0.2 m/s; jets keep v0 (independent trajectories)
%  E. Rig distribution     3.0 L/min trials with the imaging-rig 3.0 L/min distribution (09232024),
%                          the bracket for using the 4.5 L/min distribution at 3.0 L/min
%
%Physics and settings are those of the 27 Oct rerun (RerunResults_27Oct.md): 0.8 mm orifice, 200 orifices,
%CO2 gas at 1.306 atm (ideal-gas density, Sutherland viscosity), trial-mean temperature from the Excel file,
%4-node QMOM from the 20-bin histogram, volume-weighted (mixing-cup) outlet concentration and tau.
%Outputs (CSV, all tagged with runTag): Check, Hybrid, Re200Crossing, SizeSensitivity, SizeScaleToMatch,
%Dripping, RigDistribution3LPM.

clear all; close all; clc;

%% ---------------- settings ----------------
dropletFile45 = '../Experimental Data/Stack Analysis/09242024/09242024_AllDropletData.csv';   %corrected 4.5 L/min distribution (used for both flow rates in the paper)
dropletFile30 = '../Experimental Data/Stack Analysis/09232024/09232024_AllDropletData.csv';   %imaging-rig 3.0 L/min distribution (part E only)
experimentalFile = '../Experimental Data/MassTransfer Experimental Data 10222025.xlsx';
runTag = '10282026';
outDir = '../Model Results/';   %all CSV outputs go to the Model Results folder
d_orifice = 0.8/1000;          %orifice diameter (m), re-measured 27 Oct
N0 = 200;                      %number of orifices
P = (14.2+5)/14.7;             %operating pressure (atm), as in the 27 Oct rerun
nNodes = 4;                    %quadrature nodes
xspan = 0:0.01:1;              %location (m)
Re_switch = 200;               %circulating -> oscillating transition
hybridOpts = odeset('MaxStep',0.002); %small steps so the kL switch at Re = 200 is resolved
sizeScales = [0.8 0.9 1.1 1.2];
widthScales = [0.5 1.5];
dripFractions = [0.1 0.2 0.3];
dripDiameters_mm = [2.5 2.75 3.1];  %2.5-3 mm from the rig images; 3.1 mm Tate/Harkins-Brown estimate
dripV0 = [0.05 0.2];                %dripping drops leave the orifice slowly (m/s)

%% ---------------- quadrature nodes ----------------
[r45,w45,dmax45] = quadratureNodes(dropletFile45,nNodes);
fprintf('4.5 L/min distribution: radius (mm) %s, number weights %s, d_max %.3f mm\n', ...
    mat2str(r45'*1000,4),mat2str(w45',3),dmax45);
%expected 0.553 / 0.748 / 0.948 / 1.249 mm and 0.075 / 0.570 / 0.346 / 0.010

%% ---------------- trial list from the experimental file ----------------
sheets = sheetnames(experimentalFile);
trials = zeros(0,5); %Q (L/min), H (m), trial, T (C), measured C (g/L)
for s = 1:length(sheets)
    e = readtable(experimentalFile,'Sheet',sheets{s},'VariableNamingRule','preserve');
    for k = unique(e.Trial)'
        in = e.Trial==k;
        trials(end+1,:) = [e.('Flowrate (L/min)')(find(in,1)), e.('Height.m')(find(in,1)), k, ...
            mean(e.('Liquid Temperature (C)')(in)), mean(e.('CO2 (g/L)')(in))]; %#ok<SAGROW>
    end
end
trials = sortrows(trials,[1 2 3]);
nT = size(trials,1);

%% ---------------- A. reproduction check ----------------
i0 = find(trials(:,1)==4.5 & trials(:,2)==0.98 & trials(:,3)==1);
[C0,t0] = runModel(trials(i0,:),P,d_orifice,N0,r45,w45,xspan,'none',Re_switch,[],[]);
fprintf('\nA. Check: Angelo 4.5 L/min 0.98 m T1 = %.4f g/L (expect 1.1707), tau = %.3f s (expect 0.417)\n',C0,t0);
if abs(C0-1.1707)>0.002
    warning('Check value not reproduced: compare settings (orifice, pressure, gas properties, droplet file) with your RunTrial.m before using the results.')
end
writetable(table(C0,t0,'VariableNames',{'Angelo_gL','tau_s'}),[outDir 'Check_' runTag '.csv']);

%% ---------------- B. hybrid runs ----------------
lowModels = {'none','Amokrane','Ruckenstein'};
lowLabels = {'Angelo','Hybrid_Amokrane','Hybrid_Ruckenstein'};
Hyb = zeros(nT,length(lowModels));
for n = 1:nT
    for m = 1:length(lowModels)
        Hyb(n,m) = runModel(trials(n,:),P,d_orifice,N0,r45,w45,xspan,lowModels{m},Re_switch,[],hybridOpts);
    end
end
hybT = array2table([trials Hyb 100*(Hyb(:,2:end)./Hyb(:,1)-1)], 'VariableNames', ...
    [{'FlowRate_LPM','Height_m','Trial','Temp_C','Measured_gL'},lowLabels, ...
     {'Amokrane_pct','Ruckenstein_pct'}]);
disp(hybT); writetable(hybT,[outDir 'Hybrid_' runTag '.csv']);
for Q = [3 4.5]
    q = hybT.FlowRate_LPM==Q;
    fprintf('B. %.1f L/min: hybrid change vs Angelo, mean (range) %%: Amokrane %.2f (%.2f to %.2f), Ruckenstein %.2f (%.2f to %.2f)\n', Q, ...
        mean(hybT.Amokrane_pct(q)),min(hybT.Amokrane_pct(q)),max(hybT.Amokrane_pct(q)), ...
        mean(hybT.Ruckenstein_pct(q)),min(hybT.Ruckenstein_pct(q)),max(hybT.Ruckenstein_pct(q)));
    fprintf('   measured/model: Angelo %.3f, Amokrane hybrid %.3f, Ruckenstein hybrid %.3f\n', ...
        mean(hybT.Measured_gL(q)./hybT.Angelo(q)),mean(hybT.Measured_gL(q)./hybT.Hybrid_Amokrane(q)), ...
        mean(hybT.Measured_gL(q)./hybT.Hybrid_Ruckenstein(q)));
end

%where each node crosses Re = 200 (trial 1 at 0.98 m for each flow rate)
cross = zeros(0,7);
for Q = [3 4.5]
    tr = trials(trials(:,1)==Q & trials(:,2)==0.98 & trials(:,3)==1,:);
    p = props(tr(4),P);
    v0 = (Q/60000)/(pi*(d_orifice/2)^2*N0);
    for i = 1:nNodes
        r = r45(i);
        [x,y] = ode45(@(x,y) dropletODE(x,y,r,p,'none',Re_switch),xspan,[v0 0 0],hybridOpts);
        Re = 2*p.rho_g*y(:,1)*r/p.mu_g;
        k = find(Re>=Re_switch,1);
        if isempty(k); xc = NaN; tc = NaN; else; xc = x(k); tc = y(k,3); end
        cross(end+1,:) = [Q, 2*r*1000, w45(i), Re(1), xc, tc, tc/y(end,3)]; %#ok<SAGROW>
    end
end
crossT = array2table(cross,'VariableNames',{'FlowRate_LPM','Diameter_mm','NumberWeight','Re_at_orifice', ...
    'Height_Re200_m','Time_Re200_s','FractionOfFallTime'});
disp(crossT); writetable(crossT,[outDir 'Re200Crossing_' runTag '.csv']);
%note: heights are on the 0.01 m output grid

%% ---------------- C. size and width sensitivity (Angelo) ----------------
rMean = sum(w45.*r45)/sum(w45);
variants = {}; rList = {};
for s = sizeScales;  variants{end+1} = sprintf('size_x%.1f',s);  rList{end+1} = r45*s; end %#ok<SAGROW>
for s = widthScales; variants{end+1} = sprintf('width_x%.1f',s); rList{end+1} = rMean+s*(r45-rMean); end %#ok<SAGROW>
Sz = zeros(nT,length(variants));
for n = 1:nT
    base = Hyb(n,1);
    for v = 1:length(variants)
        Sz(n,v) = 100*(runModel(trials(n,:),P,d_orifice,N0,rList{v},w45,xspan,'none',Re_switch,[],hybridOpts)/base-1);
    end
end
szT = array2table([trials(:,1:4) Sz],'VariableNames',[{'FlowRate_LPM','Height_m','Trial','Temp_C'},strcat(variants,'_pct')]);
disp(szT); writetable(szT,[outDir 'SizeSensitivity_' runTag '.csv']);

%radius scale factor that makes the model match the measured mean ratio at each flow rate
match = zeros(0,3);
for Q = [3 4.5]
    q = trials(:,1)==Q;
    target = mean(trials(q,5)./Hyb(q,1));
    f = @(s) mean(arrayfun(@(n) runModel(trials(n,:),P,d_orifice,N0,r45*s,w45,xspan,'none',Re_switch,[],hybridOpts), ...
        find(q))./Hyb(q,1)) - target;
    sMatch = fzero(f,[0.6 1.2]);
    match(end+1,:) = [Q target sMatch]; %#ok<SAGROW>
    fprintf('C. %.1f L/min: measured/model %.3f would be matched by drops %.1f%% %s (all radii x%.3f)\n', ...
        Q,target,100*abs(1-sMatch),ternary(sMatch<1,'smaller','larger'),sMatch);
end
writetable(array2table(match,'VariableNames',{'FlowRate_LPM','MeasuredOverModel','RadiusScale'}),[outDir 'SizeScaleToMatch_' runTag '.csv']);

%% ---------------- D. synthetic dripping at 3.0 L/min ----------------
drip = zeros(0,10);
for n = find(trials(:,1)==3)'
    [Cb,tb] = runModel(trials(n,:),P,d_orifice,N0,r45,w45,xspan,'none',Re_switch,[],hybridOpts);
    for dd = dripDiameters_mm
        for vd = dripV0
            [Cd,td] = runModel(trials(n,:),P,d_orifice,N0,dd/2000,1,xspan,'none',Re_switch,vd,hybridOpts);
            for f = dripFractions
                Cmix = (1-f)*Cb + f*Cd; tmix = (1-f)*tb + f*td; %volume-weighted: f is a volume fraction
                drip(end+1,:) = [trials(n,1:3) dd vd f Cb Cmix 100*(Cmix/Cb-1) tmix]; %#ok<SAGROW>
            end
        end
    end
end
dripT = array2table(drip,'VariableNames',{'FlowRate_LPM','Height_m','Trial','DripDiameter_mm','DripV0_ms', ...
    'DripVolumeFraction','Angelo_gL','WithDripping_gL','Change_pct','tau_s'});
writetable(dripT,[outDir 'Dripping_' runTag '.csv']);
for f = dripFractions
    q = dripT.DripVolumeFraction==f;
    fprintf('D. %.0f%% of volume dripping: outlet C changes %.1f to %.1f%%\n',100*f,min(dripT.Change_pct(q)),max(dripT.Change_pct(q)));
end

%% ---------------- E. imaging-rig 3.0 L/min distribution ----------------
if exist(dropletFile30,'file')
    [r30,w30,dmax30] = quadratureNodes(dropletFile30,nNodes);
    fprintf('\nE. 3.0 L/min rig distribution: radius (mm) %s, number weights %s, d_max %.3f mm\n', ...
        mat2str(r30'*1000,4),mat2str(w30',3),dmax30);
    rig = zeros(0,8);
    for n = find(trials(:,1)==3)'
        [Cr,tr_] = runModel(trials(n,:),P,d_orifice,N0,r30,w30,xspan,'none',Re_switch,[],hybridOpts);
        rig(end+1,:) = [trials(n,1:5) Hyb(n,1) Cr 100*(Cr/Hyb(n,1)-1)]; %#ok<SAGROW>
    end
    rigT = array2table(rig,'VariableNames',{'FlowRate_LPM','Height_m','Trial','Temp_C','Measured_gL', ...
        'Angelo_45dist_gL','Angelo_rig30dist_gL','Change_pct'});
    disp(rigT); writetable(rigT,[outDir 'RigDistribution3LPM_' runTag '.csv']);
    fprintf('E. measured/model at 3.0 L/min: 4.5 L/min distribution %.3f, rig 3.0 L/min distribution %.3f\n', ...
        mean(rigT.Measured_gL./rigT.Angelo_45dist_gL),mean(rigT.Measured_gL./rigT.Angelo_rig30dist_gL));
else
    fprintf('\nE. skipped: %s not found\n',dropletFile30);
end
fprintf('\nDone. MATLAB %s on %s\n',version,computer);

%% ================= local functions =================
function [r,w,dmax] = quadratureNodes(file,nNodes)
%4-node QMOM from the 20-bin histogram of normalized droplet diameters (same steps as RunTrial)
d = readtable(file);
d.trueDiameter_mm = ((d.xDiameter_mm.^2).*d.yDiameter_mm).^(1/3);
[u,~,id] = unique(d.Droplet_Designation);
dm = accumarray(id,d.trueDiameter_mm,[length(u),1],@mean);
dmax = max(dm);
[hw,edges] = histcounts(dm/dmax,20,'Normalization','probability');
hx = edges(1:end-1)+diff(edges)/2;
mom = zeros(1,11);
for i = 0:10; mom(i+1) = sum(hw.*hx.^i); end
[w,x,werr] = PD_gqmom_0_1(mom,floor(length(mom)/2-1),nNodes);
if any(werr(:)~=0); warning('Nonzero realizability flag for %s',file); end
w = w(:)/sum(w); r = x(:)*dmax/2/1000; %radius (m)
end

function [C,tau] = runModel(trial,P,d_orifice,N0,r,w,xspan,lowModel,Re_switch,v0,opts)
%volume-weighted outlet concentration (g/L) and residence time (s) at the trial's height
Q = trial(1)/60000; H = trial(2);
p = props(trial(4),P);
if isempty(v0); v0 = Q/(pi*(d_orifice/2)^2*N0); end
iH = find(abs(xspan-H)<1e-9,1);
Cr = zeros(numel(r),1); tr = zeros(numel(r),1);
for i = 1:numel(r)
    if isempty(opts)
        [~,y] = ode45(@(x,y) dropletODE(x,y,r(i),p,lowModel,Re_switch),xspan,[v0 0 0]);
    else
        [~,y] = ode45(@(x,y) dropletODE(x,y,r(i),p,lowModel,Re_switch),xspan,[v0 0 0],opts);
    end
    Cr(i) = y(iH,2); tr(i) = y(iH,3);
end
vw = w(:).*r(:).^3; vw = vw/sum(vw);
C = sum(vw.*Cr); tau = sum(vw.*tr);
end

function dydx = dropletODE(~,y,r,p,lowModel,Re_switch)
%y = [v; c; t]; Angelo above Re_switch, lowModel below it ('none' = Angelo everywhere)
v = y(1);
Re = 2*p.rho_g*v*r/p.mu_g;
tau_osc = (pi/4)*sqrt(((p.rho_l*1000)*((2*r)^3))/p.s_tension); %oscillation timescale (s)
kL = 2*sqrt(p.DAB/(pi*tau_osc));                                %Angelo (m/s)
if Re < Re_switch
    switch lowModel
        case 'none'
        case 'Amokrane'
            kL = amokraneKL(r,p);                    %friction velocity at terminal velocity
        case 'Ruckenstein'
            kL = 2*sqrt(p.DAB*v/(pi*2*r));           %Ruckenstein (1967) Eq. 49
        otherwise
            error('Unknown low-Re model %s',lowModel)
    end
end
Gv = p.g - (9/2)*(p.mu_g/(r^2*p.rho_l))*v*(1+0.158*Re^(2/3)); %acceleration (m/s2)
Gc = kL*(3/r)*(p.Csat-y(2));                                  %g/L/s
dydx = [Gv/v; Gc/v; 1/v];
end

function p = props(T_C,P)
%CO2 gas and water properties at T_C (C) and P (atm), as in the 27 Oct rerun
T = T_C+273.15;
p.g = 9.8; p.rho_l = 997;
p.rho_g = P*101325*0.04401/(8.314462*T);                         %CO2 density, ideal gas (kg/m3)
p.mu_g = 1.370E-5*(T/273.15)^1.5*(273.15+240)/(T+240);          %CO2 viscosity, Sutherland (Pa s)
H = 1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T))/55342; %Henry (atm/(mol/m3))
p.Csat = P/H*44.01/1000;                                         %g/L
p.DAB = 13.942E-9*(T/227-1)^1.7094;                              %m2/s
p.s_tension = 235.8*(1-T/647.15)^1.256*(1-0.625*(1-T/647.15));   %g/s2
end

function kL = amokraneKL(r,p)
accel = @(v) p.g - (9/2)*(p.mu_g/(r^2*p.rho_l))*v*(1+0.158*(2*p.rho_g*v*r/p.mu_g)^(2/3));
vt = fzero(accel,[1E-6 50]);
Re = 2*p.rho_g*vt*r/p.mu_g;
Cd = 24/Re*(1+0.158*Re^(2/3));
v_star = sqrt(0.5*p.rho_g*vt^2*Cd/p.rho_l);
kL = 0.8*sqrt(p.DAB*v_star/(2*r));
end

function out = ternary(c,a,b)
if c; out = a; else; out = b; end
end
