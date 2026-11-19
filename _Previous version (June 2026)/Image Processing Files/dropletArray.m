clear all; close all; clc;

% specify the folder
myFolder = 'C:\Users\anbarron\Box\Research\Projects\Liquid Holdup and Mass Transfer\Image Processing and Modeling\Image Processing\TestImages';

%get list of files in the folder
filePattern = fullfile(myFolder,'*.tiff');
theFiles = dir(filePattern);

%viewing the b&w images
for i=1:length(theFiles)
    baselineFileName = theFiles(i).name;
    fullFileName=fullfile(theFiles(i).folder, baselineFileName);
    fprintf(1,'Now reading %s\n', fullFileName);

    %read each image and convert to grayscale
    originalImage = rgb2gray(imread(fullFileName));
    imageArray = im2double(originalImage);

    imshow(imageArray,[]);
    drawnow;

    %create an average image to use for scale setting
    if i==1
        sumImage = imageArray;
    else
        sumImage = sumImage + imageArray;
    end
end

%creating an average background-------------------------------------------
averageImage = sumImage/length(theFiles);
fig1 = figure('Name','1 - Average Background');
imshow(averageImage,[]);

%spacial calibration-----------------------------------------------------
[rows, columns, colorBands] = size(averageImage);

message = sprintf('Now beginning spacial calibration');
reply = questdlg(message, 'Calibrate spatially', 'OK', 'OK');

instructions = sprintf('Left click to anchor first endpoint of line.\nRight-click or double-left-click to anchor second endpoint of line.\n\nThen you will enter the real-world distance of the line.');
uiwait(msgbox(instructions));

[cx, cy, rgbValues, xi,yi] = improfile(1000);
rgbValues = squeeze(rgbValues);
distanceInPixels = sqrt( (xi(2)-xi(1)).^2 + (yi(2)-yi(1)).^2);
if length(xi) < 2
    return;
end

% Plot the line.
hold on;
lastDrawnHandle = plot(xi, yi, 'y-', 'LineWidth', 2);

userPrompt = {'Enter real world units (e.g. microns):', ...
              'Enter distance in those units:', ...
              'Enter desired scale bar length (in same units):'};
dialogTitle = 'Specify calibration information';
numberOfLines = 1;
def = {'microns', '500', '100'};
answer = inputdlg(userPrompt, dialogTitle, numberOfLines, def);
if isempty(answer)
    return;
end

units = string(answer(1));
distanceInUnits = str2double(answer(2));
distancePerPixel = distanceInUnits/distanceInPixels;
scaleBarLength_units = str2double(answer(3));
scaleBarLength_pixels = scaleBarLength_units / distancePerPixel;
scaleBarLabel = sprintf('%d %s', scaleBarLength_units, units);


message = sprintf('The distance you drew is %.2f pixels = %.0f %s.\nThe number of %s per pixel is %.2f.\nThe number of pixels per %s is %.2f',...
    distanceInPixels, distanceInUnits, units, ...
    units, distancePerPixel, ...
    units, 1/distancePerPixel);
uiwait(msgbox(message));

close all

% Choose which diagnostic images to display
showImg = struct();
showImg.backgroundSubtracted = false;
showImg.cannyEdge             = false;
showImg.noiseReduction        = false;
showImg.dropletSegmented      = fasle;
showImg.dropletOutline        = false;

% Droplet Tracking Parameters ---------------------------------------------
m=0; dropletNo = []; previousFrame = 0;
frameNo = [];
realCentroidX = [];
realCentroidY = [];
realXDiameter = [];
realYDiameter = [];
dropFrame = [];
Darkness = [];
Gradient = [];

for p=1:length(theFiles)

    baselineFileName = theFiles(p).name;
    fullFileName=fullfile(theFiles(p).folder, baselineFileName);
    fprintf(1,'Now reading %s\n', fullFileName);
    frameNo(p) = p;
    acceptCriteria = 0; %default to not ID a droplet

    %read each image and convert to grayscale
    originalImage = rgb2gray(imread(fullFileName));
    imageArray = im2double(originalImage);

    %Workflow to subtract the average background image ----------------------
    newArray = imsubtract(imageArray,averageImage);
    newArray = imcomplement(imfill(imcomplement(newArray),'holes'));
    
    if showImg.backgroundSubtracted
    fig2 = figure('Name', sprintf('2 - Background Subtracted - Frame %d', p));
    ax2 = axes(fig2);
    imshow(newArray, [], 'Parent', ax2);
    addScaleBar(ax2, scaleBarLength_pixels, scaleBarLabel, 'black');
    waitforbuttonpress; close(fig2);
    drawnow
    end
    %
    %detect the droplet with Canny edge detection
    dropletDetection = edge(newArray,'Canny',max(newArray(:))*0.75);
    %fill all the blobs
    dropletDetection = bwmorph(dropletDetection,'close');
    dropletDetection = imfill(dropletDetection,'holes');
    %remove any droplets touching the border
    dropletDetection = imclearborder(dropletDetection);
    
    if showImg.cannyEdge
    fig3 = figure('Name', sprintf('3 - Canny Edge Detection - Frame %d', p));
    ax3 = axes(fig3);
    imshow(dropletDetection, [], 'Parent', ax3);
    addScaleBar(ax3, scaleBarLength_pixels, scaleBarLabel);
    waitforbuttonpress; close(fig3);
    end

    %
    %Labeling the droplets in the image -------------------------------------
    [labeledDroplets,numberOfDrops] = bwlabel(dropletDetection,8);

    %finding relevent parameters --------------------------------------------
    dropletProps = regionprops(labeledDroplets,'all');

    %Remove non-circular droplets -------------------------------------------
    keepDroplets = [dropletProps.Extent] > 0.75 & [dropletProps.Area]>500 & [dropletProps.Circularity]>0.8 & [dropletProps.Circularity]<1.25;
    notKeepDroplets = ~keepDroplets;
    acceptDroplets = ismember(labeledDroplets,find(keepDroplets));
    
    if showImg.noiseReduction
    fig4 = figure('Name', sprintf('4 - After Noise Reduction - Frame %d', p));
    ax4 = axes(fig4);
    imshow(acceptDroplets, [], 'Parent', ax4);
    addScaleBar(ax4, scaleBarLength_pixels, scaleBarLabel);
    waitforbuttonpress; close(fig4);
    end

    acceptDropletsNo = sum(keepDroplets);

    k=[];
    acceptCriteria=[];
    meanDropletIntensity=[];
    meanGrad = [];
    
    if acceptDropletsNo ~= 0
        %generate new parameters for relevent droplets ---------------------------
        [labeledDroplets,numberOfDrops] = bwlabel(acceptDroplets,8);
        dropletProps = regionprops(labeledDroplets,'all');

        %Image Segmentation -> isolate individual droplets ------------------------
        maskedImage =originalImage;
        imshow(maskedImage,[])

        %Pulling our individual droplet images
        for k = 1:acceptDropletsNo
            dropletCenter = dropletProps(k).Centroid;
            boxRadius = dropletProps(k).MajorAxisLength/2*1.5;
            rect = [dropletCenter-boxRadius,boxRadius*2,boxRadius*2];
            subImage = imcrop(maskedImage,rect);
            subImage = imcomplement(imfill(imcomplement(subImage),'holes'));
            
            if showImg.dropletSegmented
            fig5 = figure('Name', sprintf('5 - Droplet %d Frame %d', k, p));
            ax5 = axes(fig5);
            imshow(subImage, [], 'Parent', ax5);
            addScaleBar(ax5, scaleBarLength_pixels, scaleBarLabel,'black');
            waitforbuttonpress; close(fig5);
            end

            %determining edge of Droplet
            [Gmag,Gdir]=imgradient(subImage,'sobel');
            CannyEdge = edge(subImage,'Canny',max(im2double(subImage(:)))*.75);
            CannyDroplet = bwconvhull(CannyEdge,'union');
            CannyDroplet = imfill(CannyEdge,'holes');
            
            if showImg.dropletOutline
            fig6 = figure('Name', sprintf('6 - Droplet Outline %d Frame %d', k, p));
            ax6 = axes(fig6);
            imshow(CannyEdge, [], 'Parent', ax6);
            addScaleBar(ax6, scaleBarLength_pixels, scaleBarLabel);
            waitforbuttonpress; close(fig6);
            end

            %determining droplet edge gradiant
            [yEdge,xEdge]=find(CannyEdge == 1);
            GradMag=[];
            i=[];
            for i = 1:length(yEdge)
                GradMag(i) = Gmag(yEdge(i),xEdge(i));
            end
            meanGrad(k) = mean(GradMag);
            
            %determining background and droplet intensity
            [yDroplet,xDroplet]=find(CannyDroplet == 1);
            [yBackground,xBackground]=find(CannyDroplet ==0);

            dropletPixels = subImage;
            i=[];
            for i=1:length(xBackground)
                dropletPixels(yBackground(i),xBackground(i))=255;
            end
            dropletPixels = medfilt2(dropletPixels);

            i=[];
            dropletIntensity=[];
            for i=1:length(xDroplet)
                dropletIntensity(i) = dropletPixels(yDroplet(i),xDroplet(i));
            end
            meanDropletIntensity(k) = mean(dropletIntensity);

        %Select and Size Droplets with In-Focus Criteria ----------------------
        if meanDropletIntensity(k) <65 & meanGrad(k) > 250;
            acceptCriteria(k) = 1 %droplet is dark enough and in focus
            %figure
            %imshow(subImage,[])
        else
            acceptCriteria(k) = 0
            continue
        end
        end
    else
        k=[1];
        acceptCriteria(k) = 0;
        continue
    end

    indicies = [];
    b=[];

    for b=1:k
        if acceptCriteria(b) == 1 %recording values for acceptable droplets
            m=m+1;

            centroids = vertcat(dropletProps(b).Centroid);
            centroidX = centroids(:,1);
            centroidY = centroids(:,2);

            boxes = vertcat(dropletProps(b).BoundingBox);
            xDiameter = boxes(:,3);
            yDiameter = boxes(:,4);

            dropFrame(m) = frameNo(p);
            realCentroidX(m) = centroidX;
            realCentroidY(m) = centroidY;
            realXDiameter(m) = xDiameter;
            realYDiameter(m) = yDiameter;
            Gradient(m) = meanGrad(b);
            Darkness(m) = meanDropletIntensity(b);
            

            dropletNo(m) = 0;
            if frameNo(p) == 1
                dropletNo(m) = m; %these are new droplets by definition
            elseif m == 1
                dropletNo(m) = m;
            else
                for f = max(dropletNo):-1:1
                    %determining if a droplet is the same as one in any previous frames 
                    indicies = max(find(dropletNo == f)); %identify the most recent occurrence of that droplet
                    if (frameNo(p)-dropFrame(indicies))<=50 %make sure it is within 50 slides 
                        if abs(realCentroidX(m) - realCentroidX(indicies)) < 1 %in the same x location 
                            if (realCentroidY(m) - realCentroidY(indicies)) > 0 %is below where it was previously
                                dropletNo(m) = dropletNo(indicies); %is the same droplet 
                                break
                            else
                                continue
                            end
                        else 
                            continue
                        end
                    else 
                        continue
                    end
                end
                %identifying a new droplet
                if dropletNo(m)==0
                    dropletNo(m) = max(dropletNo)+1;
                end
            end
        else
            continue
        end
    end
end

%organizing the data
dropFrame = dropFrame';
dropletNo = dropletNo';
realCentroidX = realCentroidX';
realCentroidY = realCentroidY';
realXDiameter = realXDiameter';
realYDiameter = realYDiameter';
Gradient = Gradient';
Darkness = Darkness';

%converting pixels to mm
realCentroidXmm = realCentroidX/distanceInPixels;
realCentroidYmm = realCentroidY/distanceInPixels;
realXDiametermm = realXDiameter/distanceInPixels;
realYDiametermm = realYDiameter/distanceInPixels;

%creating a table
tableNames = ["Frame No","Droplet No","X Centroid (Pixels)","X Centroid (mm)",...
    "Y Centroid (Pixels)","Y Centroid (mm)","X Diameter (Pixels)","X Diameter (mm)",...
    "Y Diameter (Pixels)","Y Diameter (mm)",'Gradient','Darkness'];
resultsTable = table(dropFrame, dropletNo, realCentroidX, realCentroidXmm,...
    realCentroidY, realCentroidYmm, realXDiameter, realXDiametermm,...
    realYDiameter, realYDiametermm, Gradient, Darkness, 'VariableNames',tableNames);

%exporting results
saveFile = fullfile(myFolder,'dropData.csv');
writetable(resultsTable,saveFile);
disp(resultsTable);

function addScaleBar(ax, scaleBarLength_pixels, scaleBarLabel, color)
    if nargin < 4
        color = 'white'; % default if no color specified
    end
    
    xlims = xlim(ax);
    ylims = ylim(ax);
    imgWidth = xlims(2) - xlims(1);
    imgHeight = ylims(2) - ylims(1);
    
    margin = 0.05;
    x_end = xlims(1) + imgWidth * (1 - margin);
    x_start = x_end - scaleBarLength_pixels;
    y_pos = ylims(1) + imgHeight * (1 - margin);
    
    hold(ax, 'on');
    plot(ax, [x_start, x_end], [y_pos, y_pos], '-', 'Color', color, 'LineWidth', 4);
    text(ax, (x_start + x_end)/2, y_pos - imgHeight*0.03, scaleBarLabel, ...
        'Color', color, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
end