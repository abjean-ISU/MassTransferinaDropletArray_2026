function dydx = ODEsys_Ruckenstein(x,y,r,T_C,P,fH,fD)
%inputs
%x: position 
%y[v,c,t]
%v: velocity (m/s) v(0)=Q/A
%c: concentration (g/L) c(0)=0
%t: time (s) t(0)=0
%r: droplet radius (m)
%T_C: experimental temperature (C)
%P: Experimental Pressure (atm)

%constants 
g = 9.8; %acceleration due to gravity (m/s2)
rho_l = 997; %density of water (kg/m3) 

%saturation concentration 
T = T_C+273.15; %Temperature (K)
H = 1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T)); %Henry's Constant for CO2 (atm)
H = H/55342; %atm/(mol/m3)

Csat = P/H; %mol/m3
Csat = Csat*44.01/1000; %g/L
Csat = Csat*fH; %Henry's law uncertainty scale factor (1 = nominal)

rho_g = P*101325*0.04401/(8.314462*T); %CO2 density, ideal gas (kg/m3)
mu_g = 1.370E-5*(T/273.15)^1.5*(273.15+222)/(T+222); %CO2 viscosity, Sutherland (Pa s)
DAB = 13.942E-9*(T/227-1)^1.7094; %diffusivity of CO2 in water (m2/s)
DAB = DAB*fD; %diffusivity uncertainty scale factor (1 = nominal)
s_tension = 235.8*(1-T/647.15)^1.256*(1-0.625*(1-T/647.15)); %water surface tenstion (g/s2)

%dynamic values 
Re = 2*rho_g*y(1)*r/mu_g; %Reynold's Number 
a = 3/r; %surface area to volume ratio 1/m
kL = 2*sqrt(DAB*y(1)/(pi*2*r)); %mass transfer coefficient m/s

Gv = g - (9/2)*(mu_g/(r^2*rho_l))*y(1)*(1+0.158*Re^(2/3)); %acceleration (m/s2)
Gc = kL*a*(Csat-y(2)); %rate of change of concentration (g/L/s)

%system of odes 
dydx(1,:) = Gv./y(1); %dvdx velocity ode (1/s)
dydx(2,:) = Gc./y(1); %dcdx concentration ode (g/L/m)
dydx(3,:) = 1/y(1); %time (s/m)

end


