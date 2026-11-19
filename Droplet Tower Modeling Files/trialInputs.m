function p = trialInputs(T_C,P,Q_LPM,d_orifice_mm)
%Inputs of one trial for the Supplement table (R2-5).
%Property equations copied from ODEsys_*.m (lines 17-29); keep identical if the ODE files change.
T = T_C+273.15;                                                     %K
H_atm = 1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T));     %Henry's constant (atm, mole-fraction basis)
H = H_atm/55342;                                                    %atm/(mol/m3)
p.Csat_gL  = P/H*44.01/1000;                                        %saturation concentration (g/L)
p.H_atm    = H_atm;
p.rho_g    = P*101325*0.04401/(8.314462*T);                         %CO2 density (kg/m3)
p.mu_g     = 1.370E-5*(T/273.15)^1.5*(273.15+222)/(T+222);          %CO2 viscosity (Pa s)
p.DAB      = 13.942E-9*(T/227-1)^1.7094;                            %CO2 diffusivity in water (m2/s)
p.sigma_Nm = 235.8*(1-T/647.15)^1.256*(1-0.625*(1-T/647.15))/1000;  %surface tension (N/m)
N0 = 200;                                                           %number of orifices
p.v0       = Q_LPM/1000/60/(pi*(d_orifice_mm/1000/2)^2*N0);         %orifice exit velocity (m/s)
end