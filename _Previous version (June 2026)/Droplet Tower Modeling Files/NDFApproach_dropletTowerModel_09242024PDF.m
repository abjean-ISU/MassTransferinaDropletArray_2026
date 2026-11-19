clear all; close all; clc;

%import cleaned droplet data 
dropletData = readtable('C:\Users\anbarron\High Speed Camera Files\09242024\AllDropletData.csv');

%compute average diameter from x and y diameters 
dropletData.trueDiameter_mm = ((dropletData.xDiameter_mm.^2).*dropletData.yDiameter_mm).^(1/3);

%designate the average droplet diameter for each distinct droplet 
%identify each unique droplet identifier tag
[unique_droplet,~,identifier_index] = unique(dropletData.Droplet_Designation);
%take the mean of each uniqud droplet
average_diameter_mm = accumarray(identifier_index, dropletData.trueDiameter_mm, [length(unique_droplet),1], @mean); 
average_diameter_mm = sort(average_diameter_mm);

%export data as a csv
writematrix(average_diameter_mm,'C:\Users\anbarron\High Speed Camera Files\09242024\09242024_DropletDiameterData.csv');

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

%using the GQMOM routine to replicate the normalized histogram
[GQMOM_hist_w,GQMOM_hist_x,GQMOM_hist_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),20);

%plotting the overlayed histograms 
figure('color','white')
manual_hist = bar(hist_abscissae,hist_weights,'k');
manual_hist.BarWidth = 1;
manual_hist.FaceAlpha = 0.5; 

hold on 
GQMOM_hist = bar(GQMOM_hist_x,GQMOM_hist_w,'w');
GQMOM_hist.BarWidth = 1;
GQMOM_hist.FaceAlpha=0.8;

legend ('Raw Data','GQMOM Routine','Location','northeast')
xlabel('Normalized Droplet Radius (r/rmax)')
ylabel('Probability')
fontsize(32,'points')

hold off 

%Using GQMOM with 4 abscissae to determine odes 
[GQMOM_model_w,GQMOM_model_x,GQMOM_model_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),4);

%converting x to radii in m
GQMOM_model_x = GQMOM_model_x*diameter_max/2/1000; %m

%system parameters
d_manifold = 9*0.0254; %9" diameter to m
d_orifice = 0.8/1000; %0.8mm orifice diameter to m
A_manifold = pi*(d_manifold/2)^2; %cross sectional manifold area (m2)
Q = 3/1000/60; % flowrate 3 LPM to m3/s
N0 = 200; %initial number of droplets 
A_orifice = pi*(d_orifice/2)^2*N0; %cross sectional orifice area (m2)

%initial values for system solutions 
xspan = 0:0.01:0.98; %location (m)
v0 = Q/A_orifice; %initial velocity (m/s)
v0N0 = N0*v0; %constant value for computing <N(x|r)> and <N(x)>
c0 = 0; %initial concentration (g/L)
t0 = 0; %initial time (s)
y0 = [v0,c0,t0]; %initial conditions 


%Angelo modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_angelo,y_angelo] = ode45(@(x,y) ODEsys_Angelo(x,y,GQMOM_model_x(i)), xspan, y0);
    v_angelo_ri(:,i) = y_angelo(:,1);
    c_angelo_ri(:,i) = y_angelo(:,2);
    t_angelo_ri(:,i) = y_angelo(:,3);
    v_angelo_ri_weighted(:,i) = v_angelo_ri(:,i)*GQMOM_model_w(i);
    c_angelo_ri_weighted(:,i) = c_angelo_ri(:,i)*GQMOM_model_w(i);
    t_angelo_ri_weighted(:,i) = t_angelo_ri(:,i)*GQMOM_model_w(i);
end

%Angelo modeling <v(x)> and <c(x)> using weights
v_angelo_x = sum(v_angelo_ri_weighted,2);
c_angelo_x = sum(c_angelo_ri_weighted,2);
t_angelo_x = sum(t_angelo_ri_weighted,2);

%Hsu modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_hsu,y_hsu] = ode45(@(x,y) ODEsys_Hsu(x,y,GQMOM_model_x(i)), xspan, y0);
    v_hsu_ri(:,i) = y_hsu(:,1);
    c_hsu_ri(:,i) = y_hsu(:,2);
    t_hsu_ri(:,i) = y_hsu(:,3);
    v_hsu_ri_weighted(:,i) = v_hsu_ri(:,i)*GQMOM_model_w(i);
    c_hsu_ri_weighted(:,i) = c_hsu_ri(:,i)*GQMOM_model_w(i);
    t_hsu_ri_weighted(:,i) = t_hsu_ri(:,i)*GQMOM_model_w(i);
end

%Hsu modeling <v(x)> and <c(x)> using weights
v_hsu_x = sum(v_hsu_ri_weighted,2);
c_hsu_x = sum(c_hsu_ri_weighted,2);
t_hsu_x = sum(t_hsu_ri_weighted,2);


%Ruckenstein modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_ruck,y_ruck] = ode45(@(x,y) ODEsys_Ruckenstein(x,y,GQMOM_model_x(i)), xspan, y0);
    v_ruck_ri(:,i) = y_ruck(:,1);
    c_ruck_ri(:,i) = y_ruck(:,2);
    t_ruck_ri(:,i) = y_ruck(:,3);
    v_ruck_ri_weighted(:,i) = v_ruck_ri(:,i)*GQMOM_model_w(i);
    c_ruck_ri_weighted(:,i) = c_ruck_ri(:,i)*GQMOM_model_w(i);
    t_ruck_ri_weighted(:,i) = t_ruck_ri(:,i)*GQMOM_model_w(i);
end

%Ruckenstein modeling <v(x)> and <c(x)> using weights
v_ruck_x = sum(v_ruck_ri_weighted,2);
c_ruck_x = sum(c_ruck_ri_weighted,2);
t_ruck_x = sum(t_ruck_ri_weighted,2);


%Amokrane modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_amokrane,y_amokrane] = ode45(@(x,y) ODEsys_Amokrane(x,y,GQMOM_model_x(i)), xspan, y0);
    v_amokrane_ri(:,i) = y_amokrane(:,1);
    c_amokrane_ri(:,i) = y_amokrane(:,2);
    t_amokrane_ri(:,i) = y_amokrane(:,3);
    v_amokrane_ri_weighted(:,i) = v_amokrane_ri(:,i)*GQMOM_model_w(i);
    c_amokrane_ri_weighted(:,i) = c_amokrane_ri(:,i)*GQMOM_model_w(i);
    t_amokrane_ri_weighted(:,i) = t_amokrane_ri(:,i)*GQMOM_model_w(i);
end

%Amokrane modeling <v(x)> and <c(x)> using weights
v_amokrane_x = sum(v_amokrane_ri_weighted,2);
c_amokrane_x = sum(c_amokrane_ri_weighted,2);
t_amokrane_x = sum(t_amokrane_ri_weighted,2);


%generating plots 

%velocity plot 
figure('color','white')
p1 = plot(x_angelo,v_angelo_ri(:,1),'-','Color','#44AA99','LineWidth',2);
hold on
p2 = plot(x_angelo,v_angelo_ri(:,2),'--','Color','#CC6677','LineWidth',2);
p3 = plot(x_angelo,v_angelo_ri(:,3),'-.','Color','#332288','LineWidth',2);
p4 = plot(x_angelo,v_angelo_ri(:,4),':','Color','#117733','LineWidth',2);

legend (strcat('r = ',num2str(GQMOM_model_x(1)*1000,2),' mm'),strcat('r = ',num2str(GQMOM_model_x(2)*1000,2),' mm'), ...
    strcat('r = ',num2str(GQMOM_model_x(3)*1000,2),' mm'),strcat('r = ',num2str(GQMOM_model_x(4)*1000,2),' mm'),'Location','northwest')
xlabel('Fall Height (m)')
ylabel('Droplet Velocity (m/s)')
fontsize(32,'points')

hold off

%concentration plot 
figure('color','white')
p1 = plot(x_ruck,c_ruck_x,'-.','Color','#332288','LineWidth',2);
hold on
p2 = plot(x_angelo,c_angelo_x,'-','Color','#44AA99','LineWidth',2);
p3 = plot(x_hsu,c_hsu_x,'--','Color','#CC6677','LineWidth',2);
p4 = plot(x_amokrane,c_amokrane_x,':','Color','#117733','LineWidth',2);

legend ('Ruckenstein et al.','Angelo et al.','Hsu et al.','Amokrane et al.','Location','northwest')
xlabel('Fall Height (m)')
ylabel('Droplet Concentration (g/L)')
fontsize(32,'points')

hold off

%system parameters 
tower_d = 0.25; %tower diameter (m)
tower_h = max(xspan); %tower height (m)
tower_v_tot = pi*(tower_d/2)^2*tower_h; %tower volume (m3)

%droplet parameters
droplet_r = sum(GQMOM_model_x.*GQMOM_model_w); %average droplet radius (m)
droplet_v = 4/3*pi*droplet_r^3; %average droplet volume (m3)

%Liquid holdup parameters
droplet_v_tot_m3 = Q*max(t_angelo_x); %total droplet (LH) volume (m3)
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

%transfer rate plot 
figure('color','white')
p1 = plot(x_ruck,v_dot_ruck,'-.','Color','#332288','LineWidth',2);
hold on
p2 = plot(x_angelo,v_dot_angelo,'-','Color','#44AA99','LineWidth',2);
p3 = plot(x_hsu,v_dot_hsu,'--','Color','#CC6677','LineWidth',2);
p4 = plot(x_amokrane,v_dot_amokrane,':','Color','#117733','LineWidth',2);

legend ('Ruckenstein et al.','Angelo et al.','Hsu et al.','Amokrane et al.','Location','northeast')
xlabel('Fall Height (m)')
ylabel('Volumetric Transfer Rate (g/L/s)')
fontsize(32,'points')

hold off

Temp = repmat(19.35,size(x_amokrane));

%Creating a table of mass transfer model results 
resultsTable = array2table([Temp,x_amokrane,c_amokrane_x,c_angelo_x,c_hsu_x,c_ruck_x], ...
    "VariableNames",{'Temperature C','Height m','Amokrane','Angelo','Hsu','Ruckenstein'});

csvFile = 'MT_Models.csv';
writetable(resultsTable,csvFile);