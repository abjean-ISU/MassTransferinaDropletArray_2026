function [t_angelo_x,resultsTable] = RunTrial(Q_LPM,T_C,P,dropletFile,csvFile,d_orifice_mm,fH,fD)

%import cleaned droplet data 
dropletData = readtable(dropletFile);

%compute average diameter from x and y diameters 
dropletData.trueDiameter_mm = ((dropletData.xDiameter_mm.^2).*dropletData.yDiameter_mm).^(1/3);

%designate the average droplet diameter for each distinct droplet 
%identify each unique droplet identifier tag
[unique_droplet,~,identifier_index] = unique(dropletData.Droplet_Designation);
%take the mean of each uniqud droplet
average_diameter_mm = accumarray(identifier_index, dropletData.trueDiameter_mm, [length(unique_droplet),1], @mean); 
average_diameter_mm = sort(average_diameter_mm);

%normalizing diameters against maximum diameter 
diameter_max = max(average_diameter_mm);
average_diameter_norm = average_diameter_mm/diameter_max; 

%computing the histogram and moments
[hist_weights,edges] = histcounts(average_diameter_norm,20,'Normalization','probability');
hist_abscissae = edges(1:end-1)+ diff(edges)/2;
for i = 0:10
    calc = hist_weights.*(hist_abscissae.^i);
    hist_moments(i+1) = sum(calc);
end

%Using GQMOM with 4 abscissae to determine odes 
[GQMOM_model_w,GQMOM_model_x,GQMOM_model_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),4);

%converting x to radii in m
GQMOM_model_x = GQMOM_model_x*diameter_max/2/1000; %m
V_ri = (4/3)*pi*GQMOM_model_x(:).^3; %droplet volume at each node (m3)
volW = GQMOM_model_w(:).*V_ri/sum(GQMOM_model_w(:).*V_ri);  %volume fraction carried by each node

%system parameters
d_manifold = 9*0.0254; %9" diameter to m
d_orifice = d_orifice_mm/1000; %0.9mm orifice diameter to m
A_manifold = pi*(d_manifold/2)^2; %cross sectional manifold area (m2)
Q = Q_LPM/1000/60; %LPM to m3/s
N0 = 200; %initial number of droplets 
A_orifice = pi*(d_orifice/2)^2*N0; %cross sectional orifice area (m2)


%initial values for system solutions 
xspan = 0:0.01:1; %location (m)
v0 = Q/A_orifice; %initial velocity (m/s)
v0N0 = N0*v0; %constant value for computing <N(x|r)> and <N(x)>
c0 = 0; %initial concentration (g/L)
t0 = 0; %initial time (s)
y0 = [v0,c0,t0]; %initial conditions 


%Angelo modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_angelo,y_angelo] = ode45(@(x,y) ODEsys_Angelo(x,y,GQMOM_model_x(i),T_C,P,fH,fD), xspan, y0);
    v_angelo_ri(:,i) = y_angelo(:,1);
    c_angelo_ri(:,i) = y_angelo(:,2);
    t_angelo_ri(:,i) = y_angelo(:,3);
    v_angelo_ri_weighted(:,i) = v_angelo_ri(:,i)*GQMOM_model_w(i);
    c_angelo_ri_weighted(:,i) = c_angelo_ri(:,i)*GQMOM_model_w(i);
    t_angelo_ri_weighted(:,i) = t_angelo_ri(:,i)*GQMOM_model_w(i);
end

%Angelo modeling <v(x)> and <c(x)> using weights
v_angelo_x = sum(v_angelo_ri_weighted,2);
c_angelo_x = c_angelo_ri*volW;
t_angelo_x = t_angelo_ri*volW;

%Hsu modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_hsu,y_hsu] = ode45(@(x,y) ODEsys_Hsu(x,y,GQMOM_model_x(i),T_C,P,fH,fD), xspan, y0);
    v_hsu_ri(:,i) = y_hsu(:,1);
    c_hsu_ri(:,i) = y_hsu(:,2);
    t_hsu_ri(:,i) = y_hsu(:,3);
    v_hsu_ri_weighted(:,i) = v_hsu_ri(:,i)*GQMOM_model_w(i);
    c_hsu_ri_weighted(:,i) = c_hsu_ri(:,i)*GQMOM_model_w(i);
    t_hsu_ri_weighted(:,i) = t_hsu_ri(:,i)*GQMOM_model_w(i);
end

%Hsu modeling <v(x)> and <c(x)> using weights
v_hsu_x = sum(v_hsu_ri_weighted,2);
c_hsu_x = c_hsu_ri*volW;
t_hsu_x = t_hsu_ri*volW;


%Ruckenstein modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_ruck,y_ruck] = ode45(@(x,y) ODEsys_Ruckenstein(x,y,GQMOM_model_x(i),T_C,P,fH,fD), xspan, y0);
    v_ruck_ri(:,i) = y_ruck(:,1);
    c_ruck_ri(:,i) = y_ruck(:,2);
    t_ruck_ri(:,i) = y_ruck(:,3);
    v_ruck_ri_weighted(:,i) = v_ruck_ri(:,i)*GQMOM_model_w(i);
    c_ruck_ri_weighted(:,i) = c_ruck_ri(:,i)*GQMOM_model_w(i);
    t_ruck_ri_weighted(:,i) = t_ruck_ri(:,i)*GQMOM_model_w(i);
end

%Ruckenstein modeling <v(x)> and <c(x)> using weights
v_ruck_x = sum(v_ruck_ri_weighted,2);
c_ruck_x = c_ruck_ri*volW;
t_ruck_x = t_ruck_ri*volW;


%Amokrane modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_amokrane,y_amokrane] = ode45(@(x,y) ODEsys_Amokrane(x,y,GQMOM_model_x(i),T_C,P,fH,fD), xspan, y0);
    v_amokrane_ri(:,i) = y_amokrane(:,1);
    c_amokrane_ri(:,i) = y_amokrane(:,2);
    t_amokrane_ri(:,i) = y_amokrane(:,3);
    v_amokrane_ri_weighted(:,i) = v_amokrane_ri(:,i)*GQMOM_model_w(i);
    c_amokrane_ri_weighted(:,i) = c_amokrane_ri(:,i)*GQMOM_model_w(i);
    t_amokrane_ri_weighted(:,i) = t_amokrane_ri(:,i)*GQMOM_model_w(i);
end

%Amokrane modeling <v(x)> and <c(x)> using weights
v_amokrane_x = sum(v_amokrane_ri_weighted,2);
c_amokrane_x = c_amokrane_ri*volW;
t_amokrane_x = t_amokrane_ri*volW;

%system parameters 
tower_d = 9*0.0254+0.1; %tower diameter (m)
tower_h = max(xspan); %tower height (m)
tower_v_tot = pi*(tower_d/2)^2*tower_h; %tower volume (m3)

%droplet parameters
droplet_r = sum(GQMOM_model_x.*GQMOM_model_w); %average droplet radius (m)
droplet_v = sum(GQMOM_model_w(:).*(4/3)*pi.*GQMOM_model_x(:).^3)/sum(GQMOM_model_w); %mean droplet volume (m3)

%Liquid holdup parameters
droplet_v_tot_m3 = Q*t_angelo_x(end); %holdup over the full 1 m column (m3); use t_angelo_x(iH) for a measurement height
droplet_v_tot_L = 1000*droplet_v_tot_m3; %total droplet (LH) volume (L)
N_tot = droplet_v_tot_m3/droplet_v; %total number of droplets
droplet_density_tot = N_tot/tower_v_tot; %total droplet number density (droplets/m3)
liquid_holdup_tot = droplet_v_tot_m3/tower_v_tot; %total liquid holdup


%predicted volumetric transfer rates
v_dot_ruck(1) = 0;
v_dot_angelo(1) = 0;
v_dot_hsu(1) = 0;
v_dot_amokrane(1) = 0;

for i = 2:length(xspan)
    v_dot_ruck(i) = (c_ruck_x(i)-c_ruck_x(i-1))/(t_ruck_x(i)-t_ruck_x(i-1)); %g/L/s
    v_dot_angelo(i) = (c_angelo_x(i)-c_angelo_x(i-1))/(t_angelo_x(i)-t_angelo_x(i-1)); %g/L/s
    v_dot_hsu(i) = (c_hsu_x(i)-c_hsu_x(i-1))/(t_hsu_x(i)-t_hsu_x(i-1)); %g/L/s
    v_dot_amokrane(i) = (c_amokrane_x(i)-c_amokrane_x(i-1))/(t_amokrane_x(i)-t_amokrane_x(i-1)); %g/L/s
end

v_dot_ruck_overall = (max(c_ruck_x)-min(c_ruck_x))/(max(t_ruck_x)-min(t_ruck_x));
v_dot_angelo_overall = (max(c_angelo_x)-min(c_angelo_x))/(max(t_angelo_x)-min(t_angelo_x));
v_dot_hsu_overall = (max(c_hsu_x)-min(c_hsu_x))/(max(t_hsu_x)-min(t_hsu_x));
v_dot_amokrane_overall = (max(c_amokrane_x)-min(c_amokrane_x))/(max(t_amokrane_x)-min(t_amokrane_x));

Temp = repmat(T_C,size(x_amokrane));

%Creating a table of mass transfer model results 
resultsTable = array2table([Temp,x_amokrane,c_amokrane_x,c_angelo_x,c_hsu_x,c_ruck_x], ...
    "VariableNames",{'Temperature C','Height m','Amokrane','Angelo','Hsu','Ruckenstein'});

if ~isempty(csvFile)
       writetable(resultsTable,csvFile); %uncertainty runs pass '' and skip writing
end
end
