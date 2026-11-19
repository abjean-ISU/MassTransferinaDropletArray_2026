clear all; close all; clc;
dropletFile = '../Experimental Data/Stack Analysis/09242024/09242024_AllDropletData.csv';   %corrected 4.5 L/min distribution

%average diameter of each droplet (same steps as RunTrial)
dropletData = readtable(dropletFile);
dropletData.trueDiameter_mm = ((dropletData.xDiameter_mm.^2).*dropletData.yDiameter_mm).^(1/3);
[unique_droplet,~,identifier_index] = unique(dropletData.Droplet_Designation);
average_diameter_mm = accumarray(identifier_index, dropletData.trueDiameter_mm, [length(unique_droplet),1], @mean);
diameter_max = max(average_diameter_mm);
average_diameter_norm = average_diameter_mm/diameter_max;

%histogram and moments
[hist_weights,edges] = histcounts(average_diameter_norm,20,'Normalization','probability');
hist_abscissae = edges(1:end-1)+ diff(edges)/2;
for i = 0:10
    hist_moments(i+1) = sum(hist_weights.*(hist_abscissae.^i));
end

%20-node beta-GQMOM reconstruction for the figure
[GQMOM_hist_w,GQMOM_hist_x,GQMOM_hist_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),20);

%Figure 6: measured histogram and beta-GQMOM reconstruction
figure('color','white')
manual_hist = bar(hist_abscissae,hist_weights,'k');
manual_hist.BarWidth = 1;
manual_hist.FaceAlpha = 0.5;
hold on
GQMOM_hist = bar(GQMOM_hist_x,GQMOM_hist_w,'w');
GQMOM_hist.BarWidth = 1;
GQMOM_hist.FaceAlpha = 0.8;
legend('Raw Data','beta-GQMOM Routine','Location','northeast')
xlabel('Normalized Droplet Radius (r/rmax)')
ylabel('Probability')
fontsize(32,'points')
hold off

%model nodes for the Supplement (R2-8)
[GQMOM_model_w,GQMOM_model_x,GQMOM_model_werror]=PD_gqmom_0_1(hist_moments,floor(length(hist_moments)/2-1),4);
nodeTable = table(GQMOM_model_x(:)*diameter_max/2,GQMOM_model_w(:),'VariableNames',{'Radius_mm','Weight'})
fprintf('d_max = %.3f mm, n = %d droplets, realizability flag = %d\n',diameter_max,length(average_diameter_mm),any(GQMOM_model_werror(:)~=0));