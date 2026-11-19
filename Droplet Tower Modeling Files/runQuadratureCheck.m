%runQuadratureCheck.m
%R2-8: how much does the 4-node quadrature matter? Compares outlet concentration, residence time,
%holdup and absorption rate from
%  - QMOM with 2, 3, 4 (the model) and 5 nodes
%  - beta-GQMOM with 8 and 20 nodes (4 moment nodes + closure)
%  - direct integration over the 20-bin histogram (bin centres and probabilities)
%  - direct integration over every measured droplet (the reference)
%  - single mean diameter and single Sauter mean diameter (d32)
%  - the 4-node model with a 1000x tighter ODE tolerance (numerical check)
%for one trial per flow rate and all four kL correlations, with the run time of each method.
%Outlet concentration and residence time are volume-weighted, as in RunTrial.
%Uses the same ODE files as RunTrial (inputs x,y,r,T_C,P,fH,fD) and PD_gqmom_0_1.m.

clear all; close all; clc;

%% ---------------- settings ----------------
dropletFile = '../Experimental Data/Stack Analysis/09242024/09242024_AllDropletData.csv';   %corrected 4.5 L/min distribution
runTag = '10272026';
%one trial per flow rate: flow rate (L/min), measurement height (m), trial-mean temperature (C) from the Excel file
cases = [4.5 0.98 24.1438;    %4.5 L/min High, trial 1
         3.0 0.98 19.3530];   %3.0 L/min High, trial 1
P = (14.2+5)/14.7;            %operating pressure (atm)
d_orifice = 0.8/1000;         %orifice diameter (m)
N0 = 200;                     %number of orifices
xspan = 0:0.01:1;             %location (m)
corrNames = {'Angelo','Hsu','Ruckenstein','Amokrane'};
odeFuns = {@ODEsys_Angelo,@ODEsys_Hsu,@ODEsys_Ruckenstein,@ODEsys_Amokrane};
tightOpts = odeset('RelTol',1e-6,'AbsTol',1e-9); %ode45 default is RelTol 1e-3, AbsTol 1e-6
%Amokrane's ODE file solves for terminal velocity at every step, so 3,670 droplets take ~1 h per case;
%it is skipped for the all-droplet reference and compared with the 20-bin histogram instead
skipAllDrops = {'Amokrane'};

%% ---------------- droplet sizes (same steps as RunTrial) ----------------
dropletData = readtable(dropletFile);
dropletData.trueDiameter_mm = ((dropletData.xDiameter_mm.^2).*dropletData.yDiameter_mm).^(1/3);
[unique_droplet,~,identifier_index] = unique(dropletData.Droplet_Designation);
average_diameter_mm = accumarray(identifier_index, dropletData.trueDiameter_mm, [length(unique_droplet),1], @mean);
diameter_max = max(average_diameter_mm);
average_diameter_norm = average_diameter_mm/diameter_max;
[hist_weights,edges] = histcounts(average_diameter_norm,20,'Normalization','probability');
hist_abscissae = edges(1:end-1)+ diff(edges)/2;
hist_moments = zeros(1,11);
for i = 0:10
    hist_moments(i+1) = sum(hist_weights.*(hist_abscissae.^i));
end
fprintf('%d droplets, d_max = %.3f mm, mean d = %.3f mm, d32 = %.3f mm\n',length(average_diameter_mm),diameter_max, ...
    mean(average_diameter_mm),sum(average_diameter_mm.^3)/sum(average_diameter_mm.^2));

%% ---------------- size representations: radius (m) and number weight ----------------
methods = {};
for n = 2:5
    [w,x,werr] = PD_gqmom_0_1(hist_moments,n,n);
    methods(end+1,:) = {sprintf('QMOM %d nodes',n),x(:)*diameter_max/2/1000,w(:)/sum(w),werr,[]}; %#ok<SAGROW>
end
for Nn = [8 20]
    [w,x,werr] = PD_gqmom_0_1(hist_moments,4,Nn);
    methods(end+1,:) = {sprintf('beta-GQMOM %d nodes',Nn),x(:)*diameter_max/2/1000,w(:)/sum(w),werr,[]}; %#ok<SAGROW>
end
methods(end+1,:) = {'Histogram 20 bins',hist_abscissae(:)*diameter_max/2/1000,hist_weights(:)/sum(hist_weights),0,[]};
methods(end+1,:) = {'All droplets (reference)',average_diameter_mm(:)/2/1000,ones(length(average_diameter_mm),1)/length(average_diameter_mm),0,[]};
methods(end+1,:) = {'Single mean diameter',mean(average_diameter_mm)/2/1000,1,0,[]};
methods(end+1,:) = {'Single Sauter d32',sum(average_diameter_mm.^3)/sum(average_diameter_mm.^2)/2/1000,1,0,[]};
[w,x,werr] = PD_gqmom_0_1(hist_moments,4,4);
methods(end+1,:) = {'QMOM 4 nodes (tight ODE tolerance)',x(:)*diameter_max/2/1000,w(:)/sum(w),werr,tightOpts};
iRef = find(strcmp(methods(:,1),'All droplets (reference)'));
iHist = find(strcmp(methods(:,1),'Histogram 20 bins'));

%nodes of the model quadrature for the Supplement
r4 = methods{3,2}; w4 = methods{3,3};
fprintf('\nModel quadrature (QMOM 4 nodes): radius (mm), number weight, volume weight\n');
vw4 = w4.*r4.^3/sum(w4.*r4.^3);
for i = 1:4; fprintf('  %.4f  %.4f  %.4f\n',r4(i)*1000,w4(i),vw4(i)); end
fprintf('Realizability flags (should all be 0): %s\n',mat2str(cell2mat(methods(:,4))'));

%% ---------------- run every method ----------------
outFile = sprintf('../Model Results/QuadratureCheck_%s.csv',runTag);
fid = fopen(outFile,'w');
fprintf(fid,'FlowRate.LPM,Height.m,Temp.C,Method,nNodes,Correlation,C.gL,C.diff.pct,tau.s,Holdup.L,Absorption.mg.s,Runtime.s\n');
for c = 1:size(cases,1)
    Q_LPM = cases(c,1); H = cases(c,2); T_C = cases(c,3);
    Q = Q_LPM/1000/60;
    v0 = Q/(pi*(d_orifice/2)^2*N0);
    iH = round(H/0.01)+1;
    res = zeros(size(methods,1),length(corrNames),4); %C, tau, runtime, nNodes
    for m = 1:size(methods,1)
        r = methods{m,2}; w = methods{m,3};
        vw = w.*r.^3/sum(w.*r.^3); %volume fraction carried by each size
        for j = 1:length(corrNames)
            if m==iRef && any(strcmp(corrNames{j},skipAllDrops))
                res(m,j,:) = [NaN, NaN, NaN, length(r)];
                continue
            end
            tic
            C_r = zeros(length(r),1); t_r = zeros(length(r),1);
            for i = 1:length(r)
                f = @(x,y) odeFuns{j}(x,y,r(i),T_C,P,1,1);
                if isempty(methods{m,5})
                    [~,y] = ode45(f,xspan,[v0,0,0]);
                else
                    [~,y] = ode45(f,xspan,[v0,0,0],methods{m,5});
                end
                C_r(i) = y(iH,2); t_r(i) = y(iH,3);
            end
            res(m,j,:) = [sum(vw.*C_r), sum(vw.*t_r), toc, length(r)];
        end
    end
    %write and print
    fprintf('\n==== %.1f L/min, %.2f m, %.2f C (v0 = %.3f m/s) ====\n',Q_LPM,H,T_C,v0);
    fprintf('diff (%%) is relative to the all-droplet integration (20-bin histogram for %s)\n',strjoin(skipAllDrops,', '));
    fprintf('%-36s %6s %-12s %8s %9s %7s %9s %9s %8s\n','Method','nodes','Correlation','C (g/L)','diff (%)','tau (s)','Holdup(L)','mg/s','time (s)');
    for m = 1:size(methods,1)
        for j = 1:length(corrNames)
            Cm = res(m,j,1); tau = res(m,j,2);
            Cref = res(iRef,j,1);
            if isnan(Cref); Cref = res(iHist,j,1); end %skipped correlation: compare with the 20-bin histogram
            dpct = 100*(Cm/Cref-1);
            fprintf(fid,'%.1f,%.2f,%.4f,%s,%d,%s,%.5f,%.3f,%.5f,%.6f,%.4f,%.3f\n',Q_LPM,H,T_C,methods{m,1},res(m,j,4), ...
                corrNames{j},Cm,dpct,tau,1000*Q*tau,1000*Q*1000*Cm,res(m,j,3));
            fprintf('%-36s %6d %-12s %8.4f %9.2f %7.4f %9.4f %9.3f %8.2f\n',methods{m,1},res(m,j,4),corrNames{j}, ...
                Cm,dpct,tau,1000*Q*tau,1000*Q*1000*Cm,res(m,j,3));
        end
    end
end
fclose(fid);
fprintf('\nWritten: %s\n',outFile);
fprintf('Hardware for the runtime column (Methods, R2-8): %s\n',computer);
