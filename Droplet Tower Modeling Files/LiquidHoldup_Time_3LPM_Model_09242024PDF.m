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
Q = 3/1000/60; % flowrate 1.2 gpm to m3/s
N0 = 200; %initial number of droplets 
A_orifice = pi*(d_orifice/2)^2*N0; %cross sectional orifice area (m2)


%initial values for system solutions 
xspan = 0:0.01:1; %location (m)
v0 = Q/A_orifice; %initial velocity (m/s)
v0N0 = N0*v0; %constant value for computing <N(x|r)> and <N(x)>
t0 = 0; %initial time (s)
y0 = [v0,t0]; %initial conditions 


% modeling <v(x|r)> and <c(x|r)>
for i=1:length(GQMOM_model_x)
    [x_LH,y_LH] = ode45(@(x,y) ODEsys_LH(x,y,GQMOM_model_x(i)), xspan, y0);
    v_LH_ri(:,i) = y_LH(:,1);
    t_LH_ri(:,i) = y_LH(:,2);
    v_LH_ri_weighted(:,i) = v_LH_ri(:,i)*GQMOM_model_w(i);
    t_LH_ri_weighted(:,i) = t_LH_ri(:,i)*GQMOM_model_w(i);
end

%modeling <v(x)> and <t(x)> using weights
v_x = sum(v_LH_ri_weighted,2);
t_x = sum(t_LH_ri_weighted,2);

%velocity plot 
figure('color','white')
p1 = plot(x_LH,v_LH_ri(:,1),'-','Color','#44AA99','LineWidth',2);
hold on
p2 = plot(x_LH,v_LH_ri(:,2),'--','Color','#CC6677','LineWidth',2);
p3 = plot(x_LH,v_LH_ri(:,3),'-.','Color','#332288','LineWidth',2);
p4 = plot(x_LH,v_LH_ri(:,4),':','Color','#117733','LineWidth',2);

legend (strcat('r = ',num2str(GQMOM_model_x(1)*1000,2),' mm'),strcat('r = ',num2str(GQMOM_model_x(2)*1000,2),' mm'), ...
    strcat('r = ',num2str(GQMOM_model_x(3)*1000,2),' mm'),strcat('r = ',num2str(GQMOM_model_x(4)*1000,2),' mm'),'Location','northwest')
xlabel('Distance from Droplet Manifold (m)')
ylabel('Droplet Velocity (m/s)')
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
droplet_v_tot_m3 = Q*max(t_x); %total droplet (LH) volume (m3)
droplet_v_tot_L = 1000*droplet_v_tot_m3; %total droplet (LH) volume (L)
N_tot = droplet_v_tot_m3/droplet_v; %total number of droplets
droplet_density_tot = N_tot/tower_v_tot; %total droplet number density (droplets/m3)
liquid_holdup_tot = droplet_v_tot_m3/tower_v_tot; %total liquid holdup

%Important parameters at different heights
fall_heights = [0.53; 0.68; 0.98];
Indicies = [find(xspan>0.52 & xspan<0.54),find(xspan>0.67 & xspan<0.69),find(xspan>0.97 & xspan<0.99)];
fall_times = t_x(Indicies);
volume_holdup_m3 = Q*fall_times;
volume_holdup_L = 1000*volume_holdup_m3;
N = volume_holdup_m3/droplet_v; %total number of droplets
a_m = [3/droplet_r;3/droplet_r;3/droplet_r];

%creating table to export 
resultsTable = array2table([fall_heights,fall_times,volume_holdup_m3,volume_holdup_L,N,a_m], ...
    "VariableNames",{'fall.height.m','fall.time.s','volume.holdup.m3','volume.holdup.L','N','SA.to.V.1/m'});

csvFile = '3LPM_LH_t_model.csv';
writetable(resultsTable,csvFile);

