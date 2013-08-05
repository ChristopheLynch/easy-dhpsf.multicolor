function [] = Combine_MulticolorSMACM()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

scrsz = get(0,'ScreenSize');


%% Ask user for relevant datafiles

[tformFile tformPath] = uigetfile({'*.mat';'*.*'},'Open 3D_Transform.mat');
if isequal(tformFile,0)
   error('User cancelled the program');
end
load([tformPath tformFile]);
% Prepare for data interpolation
x = matched_cp_reflected(:,5);
y = matched_cp_reflected(:,6);
z = matched_cp_reflected(:,7);
F_FRE = TriScatteredInterp(x,y,z,FRE_full(:,1), 'natural');
F_FRE_X = TriScatteredInterp(x,y,z,FRE_full(:,2), 'natural');
F_FRE_Y = TriScatteredInterp(x,y,z,FRE_full(:,3), 'natural');
F_FRE_Z = TriScatteredInterp(x,y,z,FRE_full(:,4), 'natural');
F_TRE = TriScatteredInterp(x,y,z,TRE_full(:,1), 'natural');
F_TRE_X = TriScatteredInterp(x,y,z,TRE_full(:,2), 'natural');
F_TRE_Y = TriScatteredInterp(x,y,z,TRE_full(:,3), 'natural');
F_TRE_Z = TriScatteredInterp(x,y,z,TRE_full(:,4), 'natural');

[whiteLightFile whiteLightPath] = uigetfile({'*.tif';'*.*'},'Open image stack with white light image');

[LocFile LocPath] = uigetfile({'*.mat';'*.*'},'Open data file file #1 with filtered molecule fits');
if isequal(LocFile,0)
   error('User cancelled the program');
end
load([LocPath LocFile]);

fileNum = 1;
LocFiles = {};
dataSets = [];

while ~isequal(LocFile,0)
    
    LocFiles = [LocFiles; {[LocPath LocFile]}];
    % load data
    load([LocPath LocFile]);
    
    %    Assemble the structure for this dataset
    dataSet.frameNum = frameNum;
    dataSet.xLoc = xLoc;
    dataSet.yLoc = yLoc;
    dataSet.zLoc = zLoc;    
    dataSet.sigmaX = sigmaX;
    dataSet.sigmaY = sigmaY;
    dataSet.sigmaZ = sigmaZ;
    dataSet.numPhotons = numPhotons;
    dataSet.meanBkgnd = meanBkgnd;
    dataSet.fidTrackX = fidTrackX;
    dataSet.fidTrackY = fidTrackY;
    dataSet.fidTrackZ = fidTrackZ;
            
    dlg_title = 'Transform dataset';
    prompt = {'Do you want to transform this dataset?'};
    def =       { 'Yes' };
    questiondialog = questdlg(prompt,dlg_title, def);
    % Handle response
    switch questiondialog
        case 'Yes'
            
            dataSet.transformedDataset = true;
            transformedData = transformData([xLoc, yLoc, zLoc],tform);
            dataSet.xLoc_transformed = transformedData(:,1);
            dataSet.yLoc_transformed = transformedData(:,2);
            dataSet.zLoc_transformed = transformedData(:,3);
            
            
            dataSet.fidTrackX_transformed = NaN(length(fidTrackX),1);
            dataSet.fidTrackY_transformed = NaN(length(fidTrackX),1);
            dataSet.fidTrackZ_transformed = NaN(length(fidTrackX),1);
            transformedData = transformData([fidTrackX(~isnan(fidTrackX)),...
                fidTrackY(~isnan(fidTrackY)),...
                fidTrackZ(~isnan(fidTrackZ))],tform);
            dataSet.fidTrackX_transformed(~isnan(fidTrackX)) = transformedData(:,1);
            dataSet.fidTrackY_transformed(~isnan(fidTrackY)) = transformedData(:,2);
            dataSet.fidTrackZ_transformed(~isnan(fidTrackZ)) = transformedData(:,3);
            
            dataSet.fidTrack_interpolated_FRE = NaN(length(fidTrackX),4);
            dataSet.fidTrack_interpolated_TRE = NaN(length(fidTrackX),4);
            dataSet.fidTrack_interpolated_FRE(~isnan(fidTrackX),:) = [...
                F_FRE(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                F_FRE_X(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                F_FRE_Y(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                F_FRE_Z(transformedData(:,1),transformedData(:,2),transformedData(:,3)) ];
            dataSet.fidTrack_interpolated_TRE(~isnan(fidTrackX),:) = [...
                F_TRE(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                F_TRE_X(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                F_TRE_Y(transformedData(:,1),transformedData(:,2),transformedData(:,3)),...
                F_TRE_Z(transformedData(:,1),transformedData(:,2),transformedData(:,3)) ];
            
        case 'No'
            
            dataSet.transformedDataset = false;
            dataSet.xLoc_transformed = NaN(length(xLoc),1);
            dataSet.yLoc_transformed = NaN(length(xLoc),1);
            dataSet.zLoc_transformed = NaN(length(xLoc),1);
            dataSet.fidTrackX_transformed = NaN(length(xLoc),1);
            dataSet.fidTrackY_transformed = NaN(length(xLoc),1);
            dataSet.fidTrackZ_transformed = NaN(length(xLoc),1);
            dataSet.fidTrack_interpolated_FRE = NaN(length(fidTrackX),4);
            dataSet.fidTrack_interpolated_TRE = NaN(length(fidTrackX),4);
            
        case 'Cancel'
            error('User cancelled the program');
    end
    
    dataSet.LocFile = LocFile;
    dataSet.LocPath = LocPath;
    
    dataSet.frameRange = frameRange;
    dataSet.zRange = zRange;
    dataSet.ampRatioLimit = ampRatioLimit;
    dataSet.fitErrorRange = fitErrorRange;
    dataSet.lobeDistBounds = lobeDistBounds;
    dataSet.sigmaBounds = sigmaBounds;
    dataSet.sigmaRatioLimit = sigmaRatioLimit;
    dataSet.numPhotonRange = numPhotonRange;

    
    % Append the structure to the previous structures in the array
    dataSets = [dataSets, dataSet];
    
    fileNum = fileNum+1;
    [LocFile LocPath] = uigetfile({'*.mat';'*.*'},...
        ['Open data file #' num2str(fileNum) ' with filtered molecule fits']);
    
end

clear LocFile LocPath ampRatioLimit dataSet def dlg_title fileNum
clear fitErrorRange frameNum frameRange lobeDistBounds meanBkgnd 
clear numPhotonRange numPhotons prompt questiondialog sigmaBounds
clear sigmaRatioLimit sigmaX sigmaY sigmaZ transformedData
clear xLoc yLoc zLoc zRange

% save('workspace.mat');
%% Show the fiducial tracks in the same coordinate system


display = 1:length(dataSets);    
color = {'green', 'red', 'blue', 'cyan', 'yellow' };

figure_h_a = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
    'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
    'DefaultAxesLineWidth',1.2);


for i = display
    if dataSets(i).transformedDataset
        scatter3(dataSets(i).fidTrackX_transformed,...
            dataSets(i).fidTrackY_transformed,...
            dataSets(i).fidTrackZ_transformed,...
            5,'filled', color{i});
%         xlim([min(dataSets(i).xLoc_transformed) max(dataSets(i).xLoc_transformed)]);
%         ylim([min(dataSets(i).yLoc_transformed) max(dataSets(i).yLoc_transformed)]);
    else
        scatter3(dataSets(i).fidTrackX,...
            dataSets(i).fidTrackY,...
            dataSets(i).fidTrackZ,...
            5,'filled', color{i});
%         xlim([min(dataSets(i).xLoc) max(dataSets(i).xLoc)]);
%         ylim([min(dataSets(i).yLoc) max(dataSets(i).yLoc)]);
    end
    
    xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
    axis vis3d equal;
    hold on
    
end
hold off

saveas(figure_h_a,'3DFidCorrelation.fig');
saveas(figure_h_a,'3DFidCorrelation.png');
%%
figure_h_b = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
    'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
    'DefaultAxesLineWidth',1.2);

for i = display
    if dataSets(i).transformedDataset
        subplot(3,1,1)
        plot(dataSets(i).fidTrackX_transformed, color{i});
        xlabel('frame number');ylabel('X position (nm)');
        hold on
        subplot(3,1,2)
        plot(dataSets(i).fidTrackY_transformed, color{i});
        xlabel('frame number');ylabel('Y position (nm)');
        hold on
        subplot(3,1,3)
        plot(dataSets(i).fidTrackZ_transformed, color{i});
        xlabel('frame number');ylabel('Z position (nm)');
        hold on
    else
        subplot(3,1,1)
        plot(dataSets(i).fidTrackX, color{i});
        xlabel('frame number');ylabel('X position (nm)');
        hold on
        subplot(3,1,2)
        plot(dataSets(i).fidTrackY, color{i});
        xlabel('frame number');ylabel('Y position (nm)');
        hold on
        subplot(3,1,3)
        plot(dataSets(i).fidTrackZ, color{i});
        xlabel('frame number');ylabel('Z position (nm)');
        hold on
    end
    hold on
    
end
hold off

saveas(figure_h_b,'XYZFidTracks.fig');
saveas(figure_h_b,'XYZFidTracks.png');
%%
figure_h_c = figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w','renderer','painters');
set(gcf,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial',...
    'DefaultTextFontSize',12,'DefaultAxesFontSize',12,...
    'DefaultAxesTickLength',[0.01 0.01],'DefaultAxesTickDir','out',...
    'DefaultAxesLineWidth',1.2);

subplot(4,1,1)
plot(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
xlabel('frame number');ylabel('X position shift (nm)');
mean = nanmean(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
stdev = nanstd(dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed);
legend(['offset = ' num2str(mean),...
    ' +/- ' num2str(stdev)])
ylim([mean-5*stdev mean+7*stdev])
hold on 
plot([0;length(dataSets(1).fidTrackX)],[mean; mean],'red','LineWidth' ,2)
hold off

subplot(4,1,2)
plot(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
xlabel('frame number');ylabel('Y position shift (nm)');
mean = nanmean(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
stdev = nanstd(dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed);
legend(['offset = ' num2str(mean),...
    ' +/- ' num2str(stdev)])
ylim([mean-5*stdev mean+7*stdev])
hold on 
plot([0;length(dataSets(1).fidTrackX)],[mean; mean],'red','LineWidth' ,2)
hold off

subplot(4,1,3)
plot(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
xlabel('frame number');ylabel('Z position shift (nm)');
mean = nanmean(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
stdev = nanstd(dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed);
legend(['offset = ' num2str(mean),...
    ' +/- ' num2str(stdev)])
ylim([mean-5*stdev mean+7*stdev])
hold on 
plot([0;length(dataSets(1).fidTrackX)],[mean; mean],'red','LineWidth' ,2)
hold off

subplot(4,1,4)
euclid_Dist = sqrt(((dataSets(1).fidTrackX-dataSets(2).fidTrackX_transformed).^2)+...
    ((dataSets(1).fidTrackY-dataSets(2).fidTrackY_transformed).^2+...
    ((dataSets(1).fidTrackZ-dataSets(2).fidTrackZ_transformed).^2)));
plot(euclid_Dist);
xlabel('frame number');ylabel('3D shift (nm)');
mean = nanmean(euclid_Dist);
stdev = nanstd(euclid_Dist);
legend(['offset  = ' num2str(mean),...
    ' +/- ' num2str(stdev)])
ylim([mean-5*stdev mean+7*stdev])
hold on 
plot([0;length(dataSets(1).fidTrackX)],[mean; mean],'red','LineWidth' ,2)
hold off


saveas(figure_h_c,'XYZFidMisregistration.fig');
saveas(figure_h_c,'XYZFidMisregistration.png');

%% Display the results
% load('workspace.mat');

useTimeColors = 0;
frameRange = [1 100000];
numPhotonRange = [0 100000];

dlg_title = 'Please Input Parameters';
prompt = {  'Pixel size (in nm)',...
    'Size of Points in reconstruction',...
    'White Light Shift X (in nm)',...
    'White Light Shift Y (in nm)',...
    'Laser Power at objective (in mW)'...
    };
def = {    '125.78', ... 
    '30', ...
    '0', ...
    '0', ...
    '9.5' ...
    };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);

nmPerPixel = str2double(inputdialog{1});
scatterSize = str2double(inputdialog{2});
wlShiftX = str2double(inputdialog{3});
wlShiftY = str2double(inputdialog{4});
powerAtObjective = str2double(inputdialog{5})/1000;

pass = 1;
anotherpass = true;

while anotherpass == true
    close all
    
    %% Plot the white light image if specified
    if whiteLightFile ~= 0
        whiteLightInfo = imfinfo([whiteLightPath whiteLightFile]);
        whiteLight = zeros(whiteLightInfo(1).Height, whiteLightInfo(1).Width);
        % average white light images together to get better SNR
        for a = 1:length(whiteLightInfo)
            whiteLight = whiteLight + double(imread([whiteLightPath whiteLightFile], ...
                'Info', whiteLightInfo));
        end
        % resize white light to the size of the ROI of the single molecule fits
        ROI_initial = [1, 1, 270, 270];
%         ROI_initial = [1, 1, 200, 200];
        whiteLight = whiteLight(ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1,ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1);
        % rescale white light image to vary from 0 to 1
        whiteLight = (whiteLight-min(whiteLight(:)))/(max(whiteLight(:))-min(whiteLight(:)));
        [xWL yWL] = meshgrid((ROI_initial(1):ROI_initial(1)+ROI_initial(3)-1) * nmPerPixel + wlShiftX, ...
            (ROI_initial(2):ROI_initial(2)+ROI_initial(4)-1) * nmPerPixel + wlShiftY);
        %         [xWL yWL] = meshgrid((1:(whiteLightInfo(1).Width)) * nmPerPixel + wlShiftX, ...
        %         (1:(whiteLightInfo(1).Height)) * nmPerPixel + wlShiftY);
    end
    
    
    %% Chose a desired parameter set for reconstruction
    
    dlg_title = 'Please Input Parameters';
    prompt = {  'Size of points in reconstruction',...
        'Temporal Color Coding',...
        'White light shift X (in nm)',...
        'White light shift Y (in nm)',...
        'First frame',...
        'Last frame',...
        'Number of photons lower bound',...
        'Number of photons upper bound',...
        };
    def = { ...
        num2str(scatterSize), ...
        num2str(useTimeColors), ...
        num2str(wlShiftX), ...
        num2str(wlShiftY), ...
        num2str(frameRange(1)), ...
        num2str(frameRange(2)), ...
        num2str(numPhotonRange(1)), ...
        num2str(numPhotonRange(2)), ...
        };
    num_lines = 1;
    inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
    
    
    scatterSize = str2double(inputdialog{1});
    useTimeColors = str2double(inputdialog{2});
    wlShiftX = str2double(inputdialog{3});
    wlShiftY = str2double(inputdialog{4});
    frameRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
    numPhotonRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];
  
    
    %% load valid xyz locations
    
%     goodFits = totalPSFfits(:,17) > 0; % totalPSFfits(:,17) > -inf;
%     badFits = totalPSFfits(:,17) < 0;
%     goodFits = goodFits & totalPSFfits(:,27) >= zRange(1) & totalPSFfits(:,27) <= zRange(2);
%     goodFits = goodFits & totalPSFfits(:,1) >= frameRange(1) & totalPSFfits(:,1) <= frameRange(2);
%     goodFits = goodFits & totalPSFfits(:,numPhotonCol) >= numPhotonRange(1) & totalPSFfits(:,numPhotonCol) <= numPhotonRange(2);
    
%     if transformedData
%         xLoc = totalPSFfits(goodFits,31);
%         yLoc = totalPSFfits(goodFits,32);
%         zLoc = totalPSFfits(goodFits,33);
%         xLoc_bad = totalPSFfits(badFits,31);
%         yLoc_bad = totalPSFfits(badFits,32);
%         zLoc_bad = totalPSFfits(badFits,33);
%     else
%         if useFidCorrections
%             xLoc = xFidCorrected(goodFits);
%             yLoc = yFidCorrected(goodFits);
%             zLoc = zFidCorrected(goodFits);
%             xLoc_bad = xFidCorrected(badFits);
%             yLoc_bad = yFidCorrected(badFits);
%             zLoc_bad = zFidCorrected(badFits);
%         else
%             xLoc = totalPSFfits(goodFits,25);
%             yLoc = totalPSFfits(goodFits,26);
%             zLoc = totalPSFfits(goodFits,27);
%             xLoc_bad = totalPSFfits(badFits,25);
%             yLoc_bad = totalPSFfits(badFits,26);
%             zLoc_bad = totalPSFfits(badFits,27);
%         end
%     end
%     zLoc = zLoc * nSample/nOil;
%     numPhotons = totalPSFfits(goodFits,21);
%     meanBkgnd = totalPSFfits(goodFits,15)*conversionFactor;
%     frameNum = totalPSFfits(goodFits,1);
%     PSFfits_bad = totalPSFfits(badFits,:);
    
    %% ask user what region to plot in superresolution image
    
    if whiteLightFile ~= 0
        xRange = xWL(1,:);
        yRange = yWL(:,1);
        % pick region that contains background
        figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');
        imagesc(xRange,yRange,whiteLight);axis image;colormap gray;
        hold on;
    end
    %scatter3(xLoc,yLoc,zLoc,1,'filled');

% display = 1:length(dataSets);    
% color = {'green', 'red', 'blue', 'cyan', 'yellow' };

for i = display
    if dataSets(i).transformedDataset
        scatter(dataSets(i).xLoc_transformed,dataSets(i).yLoc_transformed,1,'filled', color{i});
        xlim([min(dataSets(i).xLoc_transformed) max(dataSets(i).xLoc_transformed)]);
        ylim([min(dataSets(i).yLoc_transformed) max(dataSets(i).yLoc_transformed)]);
%         xMean = mean(dataSets(i).xLoc_transformed);
%         yMean = mean(dataSets(i).yLoc_transformed);
    else
        scatter(dataSets(i).xLoc,dataSets(i).yLoc,1,'filled', color{i});
        xlim([min(dataSets(i).xLoc) max(dataSets(i).xLoc)]);
        ylim([min(dataSets(i).yLoc) max(dataSets(i).yLoc)]);
%         xMean = mean(dataSets(i).xLoc);
%         yMean = mean(dataSets(i).yLoc);
    end
    
    xlabel('x (nm)');ylabel('y (nm)');
    axis ij;
    hold on
    
end

[ROI, xi, yi] = roipoly;
plot(xi, yi, 'Color','black', 'LineWidth',2)

hold off

    
    %% Display the results in 3D
    
    f = figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','k','renderer','opengl', 'Toolbar', 'figure');
    if whiteLightFile~=0
        %imagesc(xRange,yRange,whiteLight);axis image;colormap gray;hold on;
        [x,y,z] = meshgrid(xRange,yRange,[-2000 2000]);
        xslice = []; yslice = []; zslice = -600;
        h=slice(x,y,z,repmat(whiteLight,[1 1 2]),xslice,yslice,zslice,'nearest');
        set(h,'EdgeColor','none','FaceAlpha',0.75);
        colormap gray; hold on;
    end   
    
    
    %% filter out localizations outside of ROI
    
    croppedDataSets = [];
    interpolated_FREs = [];
    interpolated_TREs = [];
    plotRange = [];
    
    for i = display
        
        if dataSets(i).transformedDataset
            xLoc = dataSets(i).xLoc_transformed;
            yLoc = dataSets(i).yLoc_transformed;
            zLoc = dataSets(i).zLoc_transformed;
        else
            xLoc = dataSets(i).xLoc;
            yLoc = dataSets(i).yLoc;
            zLoc = dataSets(i).zLoc;
        end
        
        validPoints = inpolygon(xLoc,yLoc,xi, yi);
        
        xLoc = xLoc(validPoints);
        yLoc = yLoc(validPoints);
        zLoc = zLoc(validPoints);
        
        %% Assemble the croppedDataSet structure
        
        croppedDataSet.frameNum =  dataSets(i).frameNum(validPoints);
        croppedDataSet.xLoc = xLoc;
        croppedDataSet.yLoc = yLoc;
        croppedDataSet.zLoc = zLoc;
        croppedDataSet.sigmaX = dataSets(i).sigmaX(validPoints);
        croppedDataSet.sigmaY = dataSets(i).sigmaY(validPoints);
        croppedDataSet.sigmaZ = dataSets(i).sigmaZ(validPoints);
        croppedDataSet.numPhotons = dataSets(i).numPhotons(validPoints);
        croppedDataSet.meanBkgnd = dataSets(i).meanBkgnd(validPoints);
        if dataSets(i).transformedDataset == 1
            interpolated_FRE = F_FRE(xLoc,yLoc,zLoc);
            interpolated_TRE = F_TRE(xLoc,yLoc,zLoc);
            croppedDataSet.interpolated_FRE = interpolated_FRE;
            croppedDataSet.interpolated_TRE = interpolated_TRE;
            interpolated_FREs = [interpolated_FREs; interpolated_FRE];
            interpolated_TREs = [interpolated_TREs; interpolated_TRE];
        else
            croppedDataSet.interpolated_FRE = nan;
            croppedDataSet.interpolated_TRE = nan;
        end

        
        
        %% plot 3D scatterplot of localizations with white light
          
        
        scatter3(xLoc,yLoc,zLoc,scatterSize,'filled',color{i});
        axis vis3d equal;
        
        plotRange = [plotRange; [min(xLoc) max(xLoc) min(yLoc) max(yLoc) min(zLoc) max(zLoc)]];
        plotRange = max(plotRange,[],1);
        
        %         xlim([min(xLoc) max(xLoc)]);
        %         ylim([min(yLoc) max(yLoc)]);
        %         zlim([min(zLoc) max(zLoc)]);
        xlim([plotRange(1) plotRange(2)]);
        ylim([plotRange(3) plotRange(4)]);
        zlim([plotRange(5) plotRange(6)]);
        xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
        title({[num2str(length(interpolated_FREs)) ' transformed localizations'];...
            ['Mean FRE = ' num2str(mean(interpolated_FREs)) ' nm'];...
            ['Mean TRE = ' num2str(mean(interpolated_TREs)) ' nm']},...
            'color','w');
        
        set(gca,'color','k');
        set(gca,'xcolor','w');set(gca,'ycolor','w');set(gca,'zcolor','w');
     

        
        croppedDataSets = [croppedDataSets, croppedDataSet];
        
    end
    
    hold off
    
    %% Construct a questdlg with three options
    
    % f = figure;
    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    % disp('This will print immediately');
    uiwait(gcf);
    % disp('This will print after you click Continue');
%     close(f);
    
    dlg_title = 'Replot';
    prompt = {'Would you like to replot with a different parameter set?'};
    def =       { 'Yes'  };
    questiondialog = questdlg(prompt,dlg_title, def);
    % Handle response
    switch questiondialog
        case 'Yes'
            pass = pass + 1;
        case 'No'
            anotherpass = false;
        case 'Cancel'
            error('User cancelled the program');
    end
    
%     z
    
end


%% prompt to save data
[saveFile, savePath] = uiputfile({'*.mat';'*.*'},'Enter a filename to save this ROI. Otherwise, click cancel.');
if ~isequal(saveFile,0)
    save([savePath saveFile(1:length(saveFile)-4) '_multicolorSMACM.mat'],'croppedDataSets','yLoc','LocFiles','ROI','dataSets','tformFile','tformPath','validPoints',...
        'whiteLightFile','whiteLightPath','whiteLight','wlShiftX','wlShiftY', 'tform');
end
%%
% % output excel spreadsheet
% textHeader = {'frame number' ...
%     'fiduciary corrected x location (nm)' ...
%     'fiduciary corrected y location (nm)' ...
%     'fiduciary corrected z location (nm)' ...
%     'sigma x (nm)' ...
%     'sigma y (nm)' ...
%     'sigma z (nm)' ...
%     'number of photons' ...
%     'mean background photons' }
% output = [frameNum, xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons, meanBkgnd]
% xlswrite([savePath saveFile(1:length(saveFile)-4) '.xlsx'], [textHeader; ...
%     num2cell(output)], ...
%     'valid PSF fits');


%% Save figures  
saveas(gcf,[savePath saveFile(1:length(saveFile)-4) '_multicolorSMACM_3D.fig']);
close
saveas(gcf,[savePath saveFile(1:length(saveFile)-4) '_multicolorSMACM_2D.fig']);
close all



end

