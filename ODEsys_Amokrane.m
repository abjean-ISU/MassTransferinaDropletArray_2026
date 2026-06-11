function dydx = ODEsys_Amokrane(x,y,r)
%inputs
%x: position 
%y[v,c,t]
%v: velocity (m/s) v(0)=Q/A
%c: concentration (g/L) c(0)=0
%t: time (s) t(0)=0
%r: droplet radius (m)

%constants 
g = 9.8; %acceleration due to gravity (m/s2)
rho_g = 1.184; %density of air (kg/m3)
rho_l = 997; %density of water (kg/m3) 
mu_g = 1.53E-5; %viscocity of air (m2/s)
DAB = 1.91E-9; %diffusion coefficient of CO2 in water at 25C (m2/s)
s_tension = 72.8; %water surface tenstion (g/s2)

%saturation concentration 
T = 298; %Temperature (K)
H = 1/(exp(-159.854+8741.68/T+21.6694*log(T)-1.10261E-3*T)); %Henry's Constant for CO2 (atm)
H = H/55342; %atm/(mol/m3)
P = 1; %operating pressure (atm)
Csat = P/H; %mol/m3
Csat = Csat*44.01/1000; %g/L

%finding the terminal velocity 
a = @(v) g - (9/2)*(mu_g/(r^2*rho_l))*v*(1+0.158*(2*rho_g*v*r/mu_g)^(2/3));
v_terminal = fsolve(a,5);

%dynamic values 
Re_term = 2*rho_g*v_terminal*r/mu_g; % terminal velocity Reynold's Number 
Re = 2*rho_g*y(1)*r/mu_g %dynamic Reynold's Number 
Cd = 24/Re_term*(1+0.158*Re_term^(2/3)); %drag coefficinet
tau_s = 1/2*rho_g*v_terminal^2*Cd; %interfacial shear stress (kg/m-s2)
a = 3/r; %surface area to volume ratio 1/m
v_star = sqrt(tau_s/rho_l); %interfacial friction velocity (m/s)
kL = 0.8*sqrt(DAB*v_star/(2*r)); %mass transfer coefficient m/s

Gv = g - (9/2)*(mu_g/(r^2*rho_l))*y(1)*(1+0.158*Re^(2/3)); %acceleration (m/s2)
Gc = kL*a*(Csat-y(2)); %rate of change of concentration (g/L/s)




%system of odes 
dydx(1,:) = Gv./y(1); %dvdx velocity ode (1/s)
dydx(2,:) = Gc./y(1); %dcdx concentration ode (g/L/m)
dydx(3,:) = 1./y(1); %time (s/m)

end


