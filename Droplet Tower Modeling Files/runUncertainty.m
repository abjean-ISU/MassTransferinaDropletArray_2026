clear all; close all; clc;

experimentalFile = '../Experimental Data/MassTransfer Experimental Data 10222025.xlsx';
dropletFile = '../Experimental Data/Stack Analysis/09242024/09242024_AllDropletData.csv';   %corrected 4.5 L/min distribution
runTag = '10272026';
d0 = 0.8; %nominal orifice diameter (mm)

%standard uncertainties of the model inputs
uT = 0.2/sqrt(3);  %temperature (C)
uP = 0.2/14.7;     %pressure, day-to-day barometric (atm)
uQ = 0.5/sqrt(3);  %flowmeter (L/min)
uH = 0.02;         %Henry's law, fraction (Perry's vs Carroll et al. 1991)
uD = 0.05;         %diffusivity, fraction (conservative)
ud = 0.05;         %orifice diameter (mm)
%saturation concentration (g/L), same as the ODE files
Csat_fun = @(T_C,P,fH) fH*P/((1/exp(-159.854+8741.68/(T_C+273.15)+21.6694*log(T_C+273.15)-1.10261E-3*(T_C+273.15)))/55342)*44.01/1000;
inputNames = {'T','P','Q','Henry','DAB','d_o'};
corrNames = {'Amokrane','Angelo','Hsu','Ruckenstein'};

out = [];
sheets = sheetnames(experimentalFile);
for s = 1:length(sheets)
    expData = readtable(experimentalFile,'Sheet',sheets{s},'VariableNamingRule','preserve');
    Q_LPM = expData.('Flowrate (L/min)')(1);
    H = expData.('Height.m')(1); %measurement height (m)
    iH = round(H/0.01)+1; %row of xspan = 0:0.01:1 at that height
    for k = unique(expData.Trial)'
        inTrial = expData.Trial==k;
        T_C = mean(expData.('Liquid Temperature (C)')(inTrial));
        P = (mean(expData.('CO2 Pressure (psig)')(inTrial))+14.2)/14.7;
        base = {Q_LPM,T_C,P,dropletFile,'',d0,1,1}; %'' = don't write a file or plot
        dC = zeros(length(inputNames),length(corrNames));
        dK = zeros(length(inputNames),1); %change in Angelo model kLa for each input (1/s)
        for j = 1:length(inputNames)
            up = base; dn = base;
            switch j
                case 1; up{2} = T_C+uT;   dn{2} = T_C-uT;
                case 2; up{3} = P+uP;     dn{3} = P-uP;
                case 3; up{1} = Q_LPM+uQ; dn{1} = Q_LPM-uQ;
                case 4; up{7} = 1+uH;     dn{7} = 1-uH;
                case 5; up{8} = 1+uD;     dn{8} = 1-uD;
                case 6; up{6} = d0+ud;    dn{6} = d0-ud;
            end
            [t_up,Rup] = RunTrial(up{:});
            [t_dn,Rdn] = RunTrial(dn{:});
            dC(j,:) = (Rup{iH,3:6}-Rdn{iH,3:6})/2; %change in C for each correlation (g/L)
            
            %Angelo model kLa = ln(C*/(C*-C))/tau at the measurement height
            Cs_up = Csat_fun(up{2},up{3},up{7});
            Cs_dn = Csat_fun(dn{2},dn{3},dn{7});
            kLa_up = log(Cs_up/(Cs_up-Rup.Angelo(iH)))/t_up(iH);
            kLa_dn = log(Cs_dn/(Cs_dn-Rdn.Angelo(iH)))/t_dn(iH);
            dK(j) = (kLa_up-kLa_dn)/2;
        end
        out(end+1,:) = [Q_LPM,H,k,sqrt(sum(dC.^2,1)),dC(:,2)',sqrt(sum(dK.^2))]; %#ok<SAGROW>
        fprintf('%s  T%d  Angelo u = %.4f g/L\n',sheets{s},k,sqrt(sum(dC(:,2).^2)));
    end
end

%one row per trial: combined uncertainty for each correlation, then each input's share for Angelo
uTable = array2table(out,'VariableNames',[{'FlowRate.LPM','Height.m','Trial'}, ...
     strcat(corrNames,'.u.gpL'),strcat('Angelo.d',inputNames),{'Angelo.kLa.u'}]);
writetable(uTable,['../Model Results/ModelUncertainty_' runTag '.csv']);