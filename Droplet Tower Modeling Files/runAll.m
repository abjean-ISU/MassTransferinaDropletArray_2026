clear all; close all; clc;

experimentalFile = '../Experimental Data/MassTransfer Experimental Data 10222025.xlsx';
dropletFile = '../Experimental Data/Stack Analysis/09242024/09242024_AllDropletData.csv';   
runTag = '10272026';
outDir = '../Model Results';   %all outputs go to the Model Results folder

sheets = sheetnames(experimentalFile);
holdup = []; %[Q_LPM, height, trial, tau (s), holdup volume (L)]
inputs = []; %one row of model inputs per trial (Supplement, R2-5)
for s = 1:length(sheets)
    expData = readtable(experimentalFile,'Sheet',sheets{s},'VariableNamingRule','preserve');
    Q_LPM = expData.('Flowrate (L/min)')(1);
    for k = unique(expData.Trial)'
        inTrial = expData.Trial==k;
        T_C = mean(expData.('Liquid Temperature (C)')(inTrial));          %trial mean, as in DataAnalysis.R
        P = (mean(expData.('CO2 Pressure (psig)')(inTrial))+14.2)/14.7;   %atm
        folder = fullfile(outDir,[strrep(sheets{s},'LPM ','LPM CO2 ') ' ' runTag],sprintf('T%d',k));
        if ~exist(folder,'dir'); mkdir(folder); end
        t_x = RunTrial(Q_LPM,T_C,P,dropletFile,fullfile(folder,'MT_Models.csv'),0.8,1,1);
        H = expData.('Height.m')(1); %measurement height (m)
        iH = round(H/0.01)+1; %row of xspan = 0:0.01:1 at that height
        holdup(end+1,:) = [Q_LPM, H, k, t_x(iH), Q_LPM/60*t_x(iH)];
        p = trialInputs(T_C,P,Q_LPM,0.8);
        inputs(end+1,:) = [Q_LPM, H, k, T_C, P, p.Csat_gL, p.H_atm, p.DAB, p.sigma_Nm, p.rho_g, p.mu_g, p.v0]; %#ok<SAGROW>
        close all   %RunTrial draws 4 figures per trial
        fprintf('%s  T%d  %.2f C  done\n',sheets{s},k,T_C);
    end
end

Hs = [0.53;0.68;0.98];
for Qval = [3 4.5]
    rows = holdup(holdup(:,1)==Qval,:);
    tau = arrayfun(@(h) mean(rows(abs(rows(:,2)-h)<1e-6,4)),Hs);
    VL = arrayfun(@(h) mean(rows(abs(rows(:,2)-h)<1e-6,5)),Hs);
    LH = table(Hs,tau,VL,'VariableNames',{'Height.m','tau.s','volume.holdup.L'});
    writetable(LH,fullfile(outDir,sprintf('%sLPM_LH_t_model_%s.csv',num2str(Qval),runTag)));
end

inT = array2table(inputs,'VariableNames',{'FlowRate_LPM','Height_m','Trial','Temp_C','Pressure_atm','Csat_gL', ...
    'Henry_atm','DAB_m2s','SurfaceTension_Nm','rho_gas_kgm3','mu_gas_Pas','v0_ms'});
inT = sortrows(inT,{'FlowRate_LPM','Height_m','Trial'});
writetable(inT,fullfile(outDir,sprintf('PerTrialInputs_%s.csv',runTag)));