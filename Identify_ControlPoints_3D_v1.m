function [] = Identify_ControlPoints_3D_v1()
% This function pairs up corresponding localization in the two channels 

% Instrument specific parameters
nmPerPixel = 125.78;
roiSize = 270;


%% Ask user for relevant datafiles
% these are the localizations for control point candidates in the two
% channels
[reflectedFile reflectedPath] = uigetfile({'*.mat';'*.*'},'Open data file with raw molecule fits in reflected channel');
if isequal(reflectedFile,0)
   error('User cancelled the program');
end
load([reflectedPath reflectedFile]);
totalPSFfits_reflected = [frameNum, xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons meanBkgnd];

[transmittedFile transmittedPath] = uigetfile({'*.mat';'*.*'},'Open data file with raw molecule fits in transmitted channel');
if isequal(transmittedFile,0)
   error('User cancelled the program');
end
load([transmittedPath transmittedFile]);
totalPSFfits_transmitted = [frameNum, xLoc, yLoc, zLoc, sigmaX, sigmaY, sigmaZ, numPhotons meanBkgnd];
clear frameNum xLoc yLoc zLoc sigmaX sigmaY sigmaZ numPhotons meanBkgnd

% this is the logfile containing information when the frames correspond to
% stationary or moving z-positions
[logFile logPath] = uigetfile({'*.dat';'*.*'},'Open sif log file');
if isequal(logFile,0)
   error('User cancelled the program');
end

%% Find control point candidates

[ PSFfits_reflected, PSFfits_transmitted, validFrames, maxNumMeasurement] = findCPCandidates(...
    logPath, logFile, totalPSFfits_reflected, totalPSFfits_transmitted );

save('Identify_ControlPoints_3D_workspace.mat');

%% Calulate the average positions of beads at each stationary point

load('Identify_ControlPoints_3D_workspace.mat');

[ Locs_reflected, Locs_transmitted ] = avgBeadPos(...
    PSFfits_reflected, PSFfits_transmitted, validFrames, maxNumMeasurement);


outputFilePrefix = [reflectedPath '\FilteredLocalizations_std_' ...
    num2str(max(Locs_reflected(:,8))) '_' num2str(max(Locs_reflected(:,9))) '_' num2str(max(Locs_reflected(:,10))) '\'];
mkdir(outputFilePrefix);

% show the averaged bead positions
figure
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
scatter3(Locs_reflected(:,5), Locs_reflected(:,6), Locs_reflected(:,7), 'filled')
title('Inspect averaged bead positions in reflected channel. Press any key to continue')
pause
saveas(gcf,[outputFilePrefix 'Locs_reflected_3D.fig']);
figure
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
scatter3(Locs_transmitted(:,5), Locs_transmitted(:,6), Locs_transmitted(:,7), 'filled')
title('Inspect averaged bead positions in transmitted channel. Press any key to continue')
pause
saveas(gcf,[outputFilePrefix 'Locs_transmitted_3D.fig']);
close all

save('Identify_ControlPoints_3D_workspace.mat');
save([outputFilePrefix 'Identify_ControlPoints_3D_output.mat']);

%% Find control points in three different z slices
load('Identify_ControlPoints_3D_workspace.mat');

matched_cpLocs_reflected = [];
matched_cpLocs_transmitted = [];
zRanges = [ min(Locs_reflected(:,7)),       max(Locs_reflected(:,7)) ;...
            min(Locs_transmitted(:,7)),     max(Locs_transmitted(:,7)) ];
centerZ = (min(zRanges(:,2))+max(zRanges(:,1)))/2;
zStep = 500;
zSliceLimit = 150;

for zSlice = centerZ-zStep:zStep:centerZ+zStep

    % for 2D identification of control points limit the tested range to +- zSliceLimit nm

    cpLocs_reflected = Locs_reflected(Locs_reflected(:,7) < (zSlice+zSliceLimit) ...
        & Locs_reflected(:,7) > (zSlice-zSliceLimit),:);

    cpLocs_transmitted = Locs_transmitted(Locs_transmitted(:,7) < (zSlice+zSliceLimit) ...
        & Locs_transmitted(:,7) > (zSlice-zSliceLimit),:);
    cpSteps = unique(cpLocs_transmitted(:,1));

    cpPSFfits_transmitted = PSFfits_transmitted(...
        PSFfits_transmitted(:,4) < (zSlice+zSliceLimit) ...
        & PSFfits_transmitted(:,4) > (zSlice-zSliceLimit),:);
    cpFrames = unique(cpPSFfits_transmitted(:,1));

    % The first control points need to be hand picked
    if ~exist('cpChannel1_approx', 'var')
        [ cpChannel1_approx, cpChannel2_approx, selectedFrame] = ...
            handpickCPs( cpFrames,PSFfits_reflected,PSFfits_transmitted,nmPerPixel,roiSize );

        % Take hand-selected pairs and find nearest fit in the respective frame
        [cp_channel1, cp_channel2] = ...
            nearestFit(cpChannel1_approx, cpChannel2_approx,selectedFrame, PSFfits_reflected,PSFfits_transmitted);

        % Calculate the (preliminary) transform function based on the
        % handpicked control points using the MATLAB built-in
        % transformation module.
        tform = cp2tform(cp_channel1,cp_channel2,'lwm');
    end
    
    % Use this transform to transform all the average x,y locations in the
    % transmitted channel (Channel 2) to their corresponding location in
    % Channel 1.  (They may not be at the same z in channel 1).
    % Assemble the subset of control points for this slice
    [ matched_cpLocs_reflected_temp, matched_cpLocs_transmitted_temp, matchedCP ] = nearestFitCP(...
        cpLocs_reflected, cpLocs_transmitted, cpSteps, tform );

    % This function is evaluated to tell the user how well the 
    % control points can be registered within each z-Slice.
    [tform, FRE, TRE, FRE_full, TRE_full] = matlab_transformation(...
        matched_cpLocs_reflected_temp(:,5:6), matched_cpLocs_transmitted_temp(:,5:6), 'affine')

    matched_cpLocs_reflected_temp = sortrows(matched_cpLocs_reflected_temp,13);
    matched_cpLocs_reflected_temp(:,13) = matched_cpLocs_reflected_temp(:,13) + size(matched_cpLocs_reflected,1)
    matched_cpLocs_reflected = [matched_cpLocs_reflected ; matched_cpLocs_reflected_temp];

    matched_cpLocs_transmitted_temp = sortrows(matched_cpLocs_transmitted_temp,13);
    matched_cpLocs_transmitted_temp(:,13) = matched_cpLocs_transmitted_temp(:,13) + size(matched_cpLocs_transmitted,1)
    matched_cpLocs_transmitted = [matched_cpLocs_transmitted ; matched_cpLocs_transmitted_temp];

end

clear FRE FRE_full TRE TRE_full cp_channel1 cp_channel2 zSlice zSliceLimit
clear matchedCP matched_cpLocs_reflected_temp matched_cpLocs_transmitted_temp selectedFrame
save('Identify_ControlPoints_3D_workspace.mat');
save([outputFilePrefix 'Identify_ControlPoints_3D_output.mat']);

%% Evaluate a preliminary 3D transformation
load('Identify_ControlPoints_3D_workspace.mat');

% Assemble the full set of control points
% The parameter fed to this function are chosen empirically, based on
% previous results.
temp=parcluster;
matlabpool(temp,temp.NumWorkers-1);
clear temp;
[tform, FRE, TRE, FRE_full, TRE_full] = custom_transformation(...
               matched_cpLocs_reflected(:,5:7),matched_cpLocs_transmitted(:,5:7),'lwquadratic',60,'Gaussian',7,1, true);
matlabpool close

cpSteps = unique(Locs_transmitted(:,1));
[ matched_cp_reflected, matched_cp_transmitted, matchedCP ] = nearestFitCP_3D(...
    Locs_reflected, Locs_transmitted, cpSteps, tform );

matched_cp_reflected = sortrows(matched_cp_reflected,13);
matched_cp_transmitted = sortrows(matched_cp_transmitted,13);

save([outputFilePrefix 'Identify_ControlPoints_3D_output.mat']);

end

%% ---------------------------------------------------------------------------------------------
function [ PSFfits_reflected, PSFfits_transmitted, validFrames, maxNumMeasurement] = findCPCandidates(...
    logPath, logFile, totalPSFfits_reflected, totalPSFfits_transmitted )
% This function isolates frames/localizations when there was no xyz motion
% These are candidates for control point localizations

%% Ask for user input
dlg_title = 'Please Input Parameters';
prompt = {  'How many stationary frames for each position?',...
        };
def = {    '20', ... 
        };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
maxNumMeasurement = str2double(inputdialog{1});

%% Analize sif log file
sifLogData =  importdata([logPath logFile]);
motionFrames = find(sifLogData(:,1)==-1);               % This finds -1 entries in the shutters correspoding to moving frames
validPeriod = diff(motionFrames)==maxNumMeasurement+1;
validFrames = [];
validMotionFrames = motionFrames(validPeriod);

for i = 1:length(validMotionFrames)
    frame = validMotionFrames(i);
    temp = (frame+3:frame+maxNumMeasurement)';
    validFrames = [validFrames; temp];
end

numFrames = size(sifLogData,1);
PSFfits_reflected = [];
PSFfits_transmitted = [];

%% only keep localization when there was no xyz motion
for frame = 1:numFrames
    if logical(sum(frame == validFrames))
 
        temp_reflected = totalPSFfits_reflected(totalPSFfits_reflected(:,1)==frame,:);
        numBeads = (1:size(temp_reflected,1))';
        temp_reflected = [temp_reflected, numBeads];
        PSFfits_reflected = [PSFfits_reflected; temp_reflected];
        
        temp_transmitted = totalPSFfits_transmitted(totalPSFfits_transmitted(:,1)==frame,:);
        numBeads = (1:size(temp_transmitted,1))';
        temp_transmitted = [temp_transmitted, numBeads];
        PSFfits_transmitted = [PSFfits_transmitted; temp_transmitted];
        frame
    end

end

end

%% --------------------------------------------------------------------------------------------- 
function [ Locs_reflected, Locs_transmitted ] = avgBeadPos(...
    PSFfits_reflected, PSFfits_transmitted, validFrames, framesToAverage )
% This function calulates the average positions of beads at each stationary point.
% Control point localizations that are too close together ( currently 100 nm) are discarded.
% Control point localizations can be filtered according to localization
% precision and number of measurements of their localization.
% The function is specific to the current format/sequence of the log file (might need to be changed). 

transitionFrames = validFrames([1;find(diff(validFrames)-1)+1]);

% Reflected Channel
Locs_reflected = [];

for i = 1:length(transitionFrames)
    
    frameRange = [transitionFrames(i), transitionFrames(i)+(framesToAverage-3)];
    
    temp_reflected = PSFfits_reflected(PSFfits_reflected(:,1)>=frameRange(1) ...
        & PSFfits_reflected(:,1)<=frameRange(2),:);
    
    if size(temp_reflected,1) <= 1
        continue
    else
        
        numBeads = max(temp_reflected(:,10));
        meanX = zeros(numBeads,1);
        meanY = zeros(numBeads,1);
        meanZ = zeros(numBeads,1);
        stdX = zeros(numBeads,1);
        stdY = zeros(numBeads,1);
        stdZ = zeros(numBeads,1);
        meanPhotons = zeros(numBeads,1);
        numMeasurements = zeros(numBeads,1);
        
        % identify the frames that have the maximum number of beads
        maxBeadFrame = temp_reflected(temp_reflected(:,10)==numBeads,1);
        buffer = 100;                                       % this parameter might need to be adjustable 
                                                            % depending on the precision of the measurement 
                                                            % and the density of control points
        beads = 0;
        
        
        for j = 1:numBeads
            
            x = temp_reflected(temp_reflected(:,1)==maxBeadFrame(1) & temp_reflected(:,10)==j,2);
            y = temp_reflected(temp_reflected(:,1)==maxBeadFrame(1) & temp_reflected(:,10)==j,3);

            if length(x) == 0
                continue
            else
                % determine how many times each bead was measured at this position
                measurements = temp_reflected(:,2)<=x+buffer ...
                    & temp_reflected(:,2)>=x-buffer ...
                    & temp_reflected(:,3)<=y+buffer ...
                    & temp_reflected(:,3)>=y-buffer;
                numMeasurements(j) = sum(measurements);
                if numMeasurements(j) > framesToAverage-3+1  % two overlapping beads, 
                                                             % these parameters are specific to the 
                                                             % current log file format
                    continue
                end
                
                % determine statistical parameters for this bead
                meanX(j) = mean(temp_reflected(measurements,2));
                meanY(j) = mean(temp_reflected(measurements,3));
                meanZ(j) = mean(temp_reflected(measurements,4));
                
                stdX(j) = std(temp_reflected(measurements,2));
                stdY(j) = std(temp_reflected(measurements,3));
                stdZ(j) = std(temp_reflected(measurements,4));
                
                meanPhotons(j) = mean(temp_reflected(measurements,8));
                
                beads = beads +1;
            end
            
        end
        goodLoc = find(meanX);
        tempArray = [i*ones(beads,1), (1:beads)',...
            frameRange(1)*ones(beads,1), frameRange(2)*ones(beads,1),...
            meanX(goodLoc), meanY(goodLoc), meanZ(goodLoc),...
            stdX(goodLoc), stdY(goodLoc), stdZ(goodLoc),...
            meanPhotons(goodLoc), numMeasurements(goodLoc)];
        
        %throw away duplicate entries that may occur due to proximity of
        %two beads
        if ~isempty(tempArray)
            tempArray = sortrows(tempArray,5);
            temp = diff(tempArray(:,5))==0;
            badFit = logical(zeros(size(tempArray,1),1));
            badFit(1:size(temp)) = temp;
            badFit(1+1:1+size(temp)) = temp;
            tempArray = tempArray(~badFit,:);
            tempArray = sortrows(tempArray,2);
            clear temp badFit
            Locs_reflected = cat(1,Locs_reflected, tempArray);
        end

    end
    
end


% Transmitted Channel
Locs_transmitted = [];

for i = 1:length(transitionFrames)
    
    frameRange = [transitionFrames(i), transitionFrames(i)+(framesToAverage-3)];
    
    temp_transmitted = PSFfits_transmitted(PSFfits_transmitted(:,1)>=frameRange(1) ...
        & PSFfits_transmitted(:,1)<=frameRange(2),:);
    
    if size(temp_transmitted,1) <= 1
        continue
    else
        
        numBeads = max(temp_transmitted(:,10));
        meanX = zeros(numBeads,1);
        meanY = zeros(numBeads,1);
        meanZ = zeros(numBeads,1);
        stdX = zeros(numBeads,1);
        stdY = zeros(numBeads,1);
        stdZ = zeros(numBeads,1);
        meanPhotons = zeros(numBeads,1);
        numMeasurements = zeros(numBeads,1);
        
        % identify the frames that have the maximum number of beads
        maxBeadFrame = temp_transmitted(temp_transmitted(:,10)==numBeads,1);      
        beads = 0;
        
        for j = 1:numBeads
            
            x = temp_transmitted(temp_transmitted(:,1)==maxBeadFrame(1) & temp_transmitted(:,10)==j,2);
            y = temp_transmitted(temp_transmitted(:,1)==maxBeadFrame(1) & temp_transmitted(:,10)==j,3);
            
            if length(x) == 0
                continue
            else
                measurements = temp_transmitted(:,2)<=x+buffer ...
                    & temp_transmitted(:,2)>=x-buffer ...
                    & temp_transmitted(:,3)<=y+buffer ...
                    & temp_transmitted(:,3)>=y-buffer;
                numMeasurements(j) = sum(measurements);
                if numMeasurements(j) > framesToAverage-3+1  % two overlapping beads, 
                                                             % these parameters are specific to the 
                                                             % current log file format
                    continue
                end
                
                meanX(j) = mean(temp_transmitted(measurements,2));
                meanY(j) = mean(temp_transmitted(measurements,3));
                meanZ(j) = mean(temp_transmitted(measurements,4));
                
                stdX(j) = std(temp_transmitted(measurements,2));
                stdY(j) = std(temp_transmitted(measurements,3));
                stdZ(j) = std(temp_transmitted(measurements,4));
                
                meanPhotons(j) = mean(temp_transmitted(measurements,8));
                
                beads = beads +1;
            end
            
        end
        goodLoc = find(meanX);
        tempArray = [i*ones(beads,1), (1:beads)',...
            frameRange(1)*ones(beads,1), frameRange(2)*ones(beads,1),...
            meanX(goodLoc), meanY(goodLoc), meanZ(goodLoc),...
            stdX(goodLoc), stdY(goodLoc), stdZ(goodLoc),...
            meanPhotons(goodLoc), numMeasurements(goodLoc)];
        
        %throw away duplicate entries that may occur due to proximity of
        %two beads
        if ~isempty(tempArray) 
            tempArray = sortrows(tempArray,5);
            temp = diff(tempArray(:,5))==0;
            badFit = logical(zeros(size(tempArray,1),1));
            badFit(1:size(temp)) = temp;
            badFit(1+1:1+size(temp)) = temp;
            tempArray = tempArray(~badFit,:);
            tempArray = sortrows(tempArray,2);
            clear temp badFit
            Locs_transmitted = cat(1,Locs_transmitted, tempArray);
        end

    end
    
end

%% Ask for user input

scrsz = get(0,'ScreenSize');
h=figure('Position',[(scrsz(3)-1280)/2 (scrsz(4)-720)/2 1280 720],'color','w');

subplot(2,4,1)
hist(Locs_reflected(:,8),300) 
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('X Localization Precision');
subplot(2,4,2)
hist(Locs_reflected(:,9),300)
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Y Localization Precision');
subplot(2,4,3)
hist(Locs_reflected(:,10),2400)
xlim([0 20])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Z Localization Precision');
subplot(2,4,4)
hist(Locs_reflected(:,12), 30)
xlim([0 framesToAverage])
xlabel('Number of Measurements');
ylabel('Frequency');
title('Number of Measurements');

subplot(2,4,5)
hist(Locs_transmitted(:,8),300) 
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('X Localization Precision');
subplot(2,4,6)
hist(Locs_transmitted(:,9),300)
xlim([0 10])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Y Localization Precision');
subplot(2,4,7)
hist(Locs_transmitted(:,10),2400)
xlim([0 20])
xlabel('Distance (nm)');
ylabel('Frequency');
title('Z Localization Precision');
subplot(2,4,8)
hist(Locs_transmitted(:,12), 30)
xlim([0 framesToAverage])
xlabel('Number of Measurements');
ylabel('Frequency');
title('Number of Measurements');

saveas(gcf,['AveragedCPLocs.fig']);
saveas(gcf,['AveragedCPLocs.png']);

% Restricting the range of localization precisions for the control points
% localizations can ensure that only high quality data points are used
% later on, but it also limits the control point density
dlg_title = 'Please Input Parameters';
prompt = {  'Standard deviation X lower bound',...
    'Standard deviation X upper bound',...
    'Standard deviation Y lower bound',...
    'Standard deviation Y upper bound',...
    'Standard deviation Z lower bound',...
    'Standard deviation z upper bound',...
    'Number of Measurements lower bound',...
    'Number of Measurements upper bound'...
        };
def = {    '0', ... 
    '4', ... 
    '0', ... 
    '4', ... 
    '0', ...  
    '7', ... 
    num2str(framesToAverage-3+1-3),...
    num2str(framesToAverage-3+1)...
        };
num_lines = 1;
inputdialog = inputdlg(prompt,dlg_title,num_lines,def);
maxNumMeasurement = str2double(inputdialog{1});


stdXRange = [str2double(inputdialog{1}) str2double(inputdialog{2})];
stdYRange = [str2double(inputdialog{3}) str2double(inputdialog{4})];
stdZRange = [str2double(inputdialog{5}) str2double(inputdialog{6})];
numMeasurementsRange = [str2double(inputdialog{7}) str2double(inputdialog{8})];

% Apply the chosen Filters 
Locs_reflected = Locs_reflected(...
    Locs_reflected(:,8)>stdXRange(1) & ...
    Locs_reflected(:,8)<=stdXRange(2) & ...
    Locs_reflected(:,9)>stdYRange(1) & ...
    Locs_reflected(:,9)<=stdYRange(2) & ...
    Locs_reflected(:,10)>stdZRange(1) & ...
    Locs_reflected(:,10)<=stdZRange(2) & ...
    Locs_reflected(:,12)>=numMeasurementsRange(1) & ...
    Locs_reflected(:,12)<=numMeasurementsRange(2) ,...
    :);

 Locs_transmitted = Locs_transmitted(...
    Locs_transmitted(:,8)>stdXRange(1) & ...
    Locs_transmitted(:,8)<=stdXRange(2) & ...
    Locs_transmitted(:,9)>stdYRange(1) & ...
    Locs_transmitted(:,9)<=stdYRange(2) & ...
    Locs_transmitted(:,10)>stdZRange(1) & ...
    Locs_transmitted(:,10)<=stdZRange(2) & ...
    Locs_transmitted(:,12)>=numMeasurementsRange(1) & ...
    Locs_transmitted(:,12)<=numMeasurementsRange(2) ,...
    :);   


end

%% --------------------------------------------------------------------------------------------- 
function [ cpChannel1_approx, cpChannel2_approx, selectedFrame ] = handpickCPs(...
    cpFrames,fits_reflected,fits_transmitted,nmPerPixel,roiSize  )
% This function shows frames side by side to let the user handpick control points
% The GUI interface and how many control points need to be selected could be optimized
% Also the user should be allowed to terminate the handpicking himself
%%% TODO: need to make sure that fits are NOT used when canceling:
%%% numAssigned seems to keep increasing, which reduces total pairs chosen
%%% before that routine quits if you change your mind about a pair
scrsz = get(0,'ScreenSize');

% Preallocate for vectors
number_handpicked = 15;
cpChannel1_approx = zeros(number_handpicked,2);
cpChannel2_approx = zeros(number_handpicked,2);
framechannels = zeros(number_handpicked,1); numAssigned = 0;  isOKAY = 0;
selectedFrame = zeros(number_handpicked,1);


for i = 1:10:length(cpFrames)
    goodFitFrame = 0;
    grayBox = zeros(1,4);
    close all
    figure('Position',[(scrsz(3)-1280)/2+1 (scrsz(4)-720)/2 1280 720],'color','w',...
        'Toolbar','figure');
    set(gcf,'DefaultTextFontSize',12,'DefaultAxesFontSize',12);
    showFrameAndFits(fits_reflected,fits_transmitted,cpFrames,i,grayBox,roiSize,nmPerPixel);
    subplot(1,2,1)
    title({'Is this a good frame?';'Hit enter if yes, click if no'})
    subplot(1,2,2)
    title({['Assigned Control Points so far: ' num2str(numAssigned)]})
    
    goodFitFrame = waitforbuttonpress;
    if goodFitFrame == 1
        isOKAY = 1;
        while isOKAY == 1
            [x,y,isOKAY] = clickPairs(numAssigned,isOKAY);
            if isOKAY == 0
                break
            end
            title('Pick some more!')
            cpChannel1_approx(numAssigned+1,:) = [x(1),y(1)];
            cpChannel2_approx(numAssigned+1,:) = [x(2),y(2)];
            framechannels(numAssigned+1,:) = cpFrames(i);
            selectedFrame(numAssigned+1) = cpFrames(i);
            grayBox(numAssigned+1,:) = [x(1),y(1),x(2),y(2)];
            numAssigned = numAssigned + 1;
            showFrameAndFits(fits_reflected,fits_transmitted,cpFrames,i,grayBox,roiSize,nmPerPixel);
            subplot(1,2,2)
            title({['Assigned Control Points so far: ' num2str(numAssigned)]})
        end
    elseif numAssigned >= number_handpicked
        break
    else
        continue
    end
end

end

%% --------------------------------------------------------------------------------------------- 
function showFrameAndFits(...
    PSFfits_reflected,PSFfits_transmitted, cpFrames, frame, grayBox, roiSize, nmPerPixel)
% This function shows the control point localizations  frame-by-frame next to each other 

if frame<1 || frame>length(cpFrames)
    return;
end

CCDChipSize = 512;
frameCol = 1;
nPhotonsCol = 8;
xCenterCol = 2;
yCenterCol = 3;

numFits_reflected = size(PSFfits_reflected(PSFfits_reflected(:,1)==cpFrames(frame),:),1);
numFits_transmitted = size(PSFfits_transmitted(PSFfits_transmitted(:,1)==cpFrames(frame),:),1);
numFits = max(numFits_reflected, numFits_transmitted);
markerColors = jet(numFits);


% First subplot (Channel1=Reflected)
temp = sortrows(PSFfits_reflected(PSFfits_reflected(:,frameCol)==cpFrames(frame),:),-nPhotonsCol);
subplot(1,2,1)
hold on
for a = 1:numFits_reflected
    scatter(temp(a,xCenterCol),temp(a,yCenterCol),temp(a,nPhotonsCol)/temp(1,nPhotonsCol)*100,...
        'MarkerFaceColor', markerColors(numFits+1-a,:),...
        'MarkerEdgeColor', markerColors(numFits+1-a,:))
end
xlim([0 roiSize*nmPerPixel])
ylim([0 roiSize*nmPerPixel])
axis square
xlabel({'Reflected Channel'; ['Frame: ' num2str(cpFrames(frame))];[num2str(numFits_reflected), ' SM localizations']})
plot([0 roiSize*nmPerPixel], [roiSize*nmPerPixel 0])
plot([0 roiSize*nmPerPixel], [0 roiSize*nmPerPixel])

if sum(grayBox)>0
    for k=1:size(grayBox,1)
        plot(grayBox(k,1),grayBox(k,2),'square','LineWidth',2,'MarkerSize',20)
    end
end
hold off

% Second subplot (Channel2=Transmitted)
temp = sortrows(PSFfits_transmitted(PSFfits_transmitted(:,frameCol)==cpFrames(frame),:),-nPhotonsCol);
subplot(1,2,2)
hold on
for a = 1:numFits_transmitted
    scatter(temp(a,xCenterCol),temp(a,yCenterCol),temp(a,nPhotonsCol)/temp(1,nPhotonsCol)*100,...
        'MarkerFaceColor', markerColors(numFits+1-a,:),...
        'MarkerEdgeColor', markerColors(numFits+1-a,:))
end

xlim([(CCDChipSize-roiSize)*nmPerPixel CCDChipSize*nmPerPixel])
ylim([(CCDChipSize-roiSize)*nmPerPixel CCDChipSize*nmPerPixel])
axis square
xlabel({'Transmitted Channel'; ['Frame: ' num2str(cpFrames(frame))];[num2str(numFits_transmitted), ' SM localizations']})
plot([(512-roiSize)*nmPerPixel CCDChipSize*nmPerPixel],...
    [CCDChipSize*nmPerPixel (512-roiSize)*nmPerPixel])
plot([(512-roiSize)*nmPerPixel CCDChipSize*nmPerPixel],...
    [(512-roiSize)*nmPerPixel CCDChipSize*nmPerPixel])

if sum(grayBox)>0
    for k=1:size(grayBox,1)
        plot(grayBox(k,3),grayBox(k,4),'square','LineWidth',2,'MarkerSize',20)
    end
end
hold off

end

%% ---------------------------------------------------------------------------------------------
function [x, y, isOKAY] = clickPairs(numAssigned, isOKAY)
% This function lets the user click on control point pairs (Left-->Right) 

clear x y
subplot(1,2,1)
title({'Pick two molecules';'(Left plot, then Right plot)'})
[x,y] = ginput(2);
subplot(1,2,1)
text(x(1),y(1),['CP',num2str(numAssigned)],'Color',[0,0,0])
subplot(1,2,2)
text(x(2),y(2),['CP',num2str(numAssigned)],'Color',[0,0,0])
subplot(1,2,1)
title({'Does the following look okay?';'Hit enter if OK, click if NO'})
isOKAY = waitforbuttonpress();
end

%% --------------------------------------------------------------------------------------------- 
function [cp_channel1, cp_channel2] = nearestFit(cp_channel1_approx, cp_channel2_approx, selectedFrame, c1_allfits, c2_allfits)
% This function takes hand-selected pairs and finds the position of the nearest localization in the respective frame

frameCol = 1;
xCenterCol = 2;
yCenterCol = 3;

% do channel 1
[numMatch,dim] = size(cp_channel1_approx);
cp_channel1 = zeros(numMatch);

for i=1:numMatch
    frameFits1 = c1_allfits((c1_allfits(:,frameCol) == selectedFrame(i)),:);
    totalDifference = sqrt((cp_channel1_approx(i,1)-frameFits1(:,xCenterCol)).^2 + ...
        (cp_channel1_approx(i,2)-frameFits1(:,yCenterCol)).^2);
    [~,matchedCoord] = min(totalDifference);
    cp_channel1(i,:) = [frameFits1(matchedCoord,xCenterCol),frameFits1(matchedCoord,yCenterCol)];
    clear totalDifference frameFits1 matchedCoord
end

% do channel 2
[numMatch,dim] = size(cp_channel2_approx);
cp_channel2 = zeros(numMatch);

for i=1:numMatch
    frameFits2 = c2_allfits((c2_allfits(:,frameCol) == selectedFrame(i)),:);
    totalDifference = sqrt((cp_channel2_approx(i,1)-frameFits2(:,xCenterCol)).^2 + ...
        (cp_channel2_approx(i,2)-frameFits2(:,yCenterCol)).^2);
    [~,matchedCoord] = min(totalDifference);
    cp_channel2(i,:) = [frameFits2(matchedCoord,xCenterCol),frameFits2(matchedCoord,yCenterCol)];
    clear totalDifference frameFits1 matchedCoord
end

end

%% ---------------------------------------------------------------------------------------------
function [ matched_cpLocs_reflected, matched_cpLocs_transmitted, matchedCP ] = nearestFitCP(...
    cpLocs_reflected, cpLocs_transmitted, cpSteps, tform )
% This function identifies the remaining control point pairs based on the
% structure tform (2D transformation)
% If the target localization is within 60 nm of the transformed
% localization the control point pair is kept. 

cpLocs_reflected = [cpLocs_reflected, nan(length(cpLocs_reflected),1)];
cpLocs_transmitted = [cpLocs_transmitted, nan(length(cpLocs_transmitted),1)];
matchedCP = 0;

for i = 1:length(cpSteps)
    i
    temp = cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i),5:6);
    % Use tform to transform all the average x,y locations in the
    % transmitted channel (Channel 2) to their corresponding location in
    % Channel 1.  (They may not be at the same z in channel 1).
    trans_cpLocs_transmitted_temp = tforminv(tform, temp);
    tempX = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),5);
    tempY = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),6);
    
    for j = 1:size(trans_cpLocs_transmitted_temp,1)
        
        residual1 = tempX - trans_cpLocs_transmitted_temp(j,1);
        residual2 = tempY - trans_cpLocs_transmitted_temp(j,2);
        totalDifference = sqrt((residual1).^2 + (residual2).^2);
        [value,matchedCoord] = min(totalDifference);
        
        if value <= 60   % in units of nm, this parameter might need to be made adjustable.
            matchedCP = matchedCP + 1;
            if isnan(cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13))
                cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13)= matchedCP;
                cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i) & ...
                    cpLocs_transmitted(:,5)==temp(j,1),13)= matchedCP;
            else
                error('Attempting to assign the same control point twice.  Aborting Execution')
            end
            
        end
    end
end

% Align the matched control points between the two channels
matched_cpLocs_reflected = cpLocs_reflected(~isnan(cpLocs_reflected(:,13)),:);
matched_cpLocs_reflected = sortrows(matched_cpLocs_reflected,13);

matched_cpLocs_transmitted = cpLocs_transmitted(~isnan(cpLocs_transmitted(:,13)),:);
matched_cpLocs_transmitted = sortrows(matched_cpLocs_transmitted,13);

if size(matched_cpLocs_reflected,1) ~= size(matched_cpLocs_transmitted,1)
    error('Unequal number of matched control points')
end

end

%% ---------------------------------------------------------------------------------------------
function [ matched_cpLocs_reflected, matched_cpLocs_transmitted, matchedCP ] = nearestFitCP_3D(...
    cpLocs_reflected, cpLocs_transmitted, cpSteps, tform )
% This function identifies the remaining control point pairs in the entire image domain
% based on the structure tform (3D transformation)
% If the target localization is within 50 nm of the transformed
% localization the control point pair is kept. 

cpLocs_reflected = [cpLocs_reflected, nan(length(cpLocs_reflected),1)];
cpLocs_transmitted = [cpLocs_transmitted, nan(length(cpLocs_transmitted),1)];
matchedCP = 0;

for i = 1:length(cpSteps)
    cpSteps(i)
    temp = cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i),5:7);
    trans_cpLocs_transmitted_temp = transformData(temp,tform);
    tempX = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),5);
    tempY = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),6);
    tempZ = cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i),7);
    
    for j = 1:size(trans_cpLocs_transmitted_temp,1)
        
        residual1 = tempX - trans_cpLocs_transmitted_temp(j,1);
        residual2 = tempY - trans_cpLocs_transmitted_temp(j,2);
        residual3 = tempZ - trans_cpLocs_transmitted_temp(j,3);
        totalDifference = sqrt((residual1).^2 + (residual2).^2 + (residual3).^2);
        [value,matchedCoord] = min(totalDifference);
        
        if value <= 50   % in units of nm
            matchedCP = matchedCP + 1;
            if isnan(cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13))
                
                cpLocs_reflected(cpLocs_reflected(:,1)==cpSteps(i) & ...
                    cpLocs_reflected(:,5)==tempX(matchedCoord),13)= matchedCP;
                cpLocs_transmitted(cpLocs_transmitted(:,1)==cpSteps(i) & ...
                    cpLocs_transmitted(:,5)==temp(j,1),13)= matchedCP;
                
            else
                error('Attempting to assign the same control point twice.  Aborting Execution')
            end
            
        end
    end
end

% Align the matched control points between the two channels
matched_cpLocs_reflected = cpLocs_reflected(~isnan(cpLocs_reflected(:,13)),:);
matched_cpLocs_reflected = sortrows(matched_cpLocs_reflected,13);

matched_cpLocs_transmitted = cpLocs_transmitted(~isnan(cpLocs_transmitted(:,13)),:);
matched_cpLocs_transmitted = sortrows(matched_cpLocs_transmitted,13);

if size(matched_cpLocs_reflected,1) ~= size(matched_cpLocs_transmitted,1)
    error('Unequal number of matched control points')
end

end

%% --------------------------------------------------------------------------------------------- 
function [tform, FRE, TRE, FRE_full, TRE_full] = matlab_transformation(...
    cp_channel1, cp_channel2 , tform_mode)
% This function computes a 2D transformation using the given set of control points
% and also computes the associated FRE and TRE values.
% The code is adapted from "Single-Molecule High Resolution Colocalization of Single
% Probes" by L. Stirling Churchman and James A. Spudich in Cold Spring
% Harbor Protocols (2012), doi: 10.1101/pdb.prot067926
% Modified Definitions of FRE and TRE - Andreas Gahlmann, 20110511

% Transform the data using the cp2tform command
tform = cp2tform(cp_channel1,cp_channel2,tform_mode);

% Calculate the metrices to estimate the error associated with this
% Calculate the fiducial registration error (FRE)
trans_cp_channel2 = tforminv(cp_channel2,tform);

% FRE = sqrt(sum(sum((cp_channel1-trans_cp_channel2).^2))/(length(cp_channel1)));
% FRE_full = ((cp_channel1-trans_cp_channel2).^2)/(length(cp_channel1));

FRE_full = sqrt(sum((cp_channel1-trans_cp_channel2).^2,2));
FRE = mean(FRE_full)

% Calculate the target registration error (TRE)
number_cp = length(cp_channel1); % find the number of control points
% Loop through the control points
for i=1:number_cp
    i
    remain_cp = [1:i-1 i+1:number_cp]; % take out that control point
    
    % Calculate the transformation without the ith control point
    tform = cp2tform(cp_channel1(remain_cp,:),cp_channel2(remain_cp,:),tform_mode);
    
    % Transform left out control point with above found transformation
    trans_cp_channel2(i,:) = tforminv(cp_channel2(i,:),tform);
    
end

% TRE = sqrt(sum(sum((cp_channel1 - trans_cp_channel2).^2))/(length(cp_channel1)));
% TRE_full = ((cp_channel1 - trans_cp_channel2).^2)/(length(cp_channel1));
TRE_full = sqrt(sum((cp_channel1-trans_cp_channel2).^2,2));
TRE = mean(TRE_full(TRE_full<100000));

% Restore the full transform function again
tform = cp2tform(cp_channel1,cp_channel2,tform_mode);
channel2_trans = tforminv(cp_channel2,tform);

% show the results
figure
distlimit = 30;

subplot(2,2,1)
scatter(cp_channel1(:,1), cp_channel1(:,2))
title({'Reflected Channel';'Channel 1'})
hold on
scatter(channel2_trans(:,1), channel2_trans(:,2), 10, 'filled')
hold off
axis square
subplot(2,2,2)
scatter(cp_channel2(:,1), cp_channel2(:,2))
title({'Transmitted Channel';'Channel 2'})
axis square

subplot(2,2,3)
hist(FRE_full(FRE_full<=distlimit), 20)
title({['Target Registration Error']; [tform_mode ' Transformation']});
xlabel('Distance (nm)');
ylabel('Frequency');
xlim([0 distlimit]);
legend(['Mean = ' num2str(FRE, 3) ' nm']);

subplot(2,2,4)
hist(TRE_full(TRE_full<distlimit),20)
title({['Fiducial Registration Error']; [tform_mode ' Transformation']});
xlabel('Distance (nm)');
ylabel('Frequency');
xlim([0 distlimit]);
legend(['Mean = ' num2str(TRE, 3) ' nm']);

% saveas(gcf,['z0nm_stats_' tform_mode '.fig']);
% saveas(gcf,['z0nm_stats_' tform_mode '.png']);

end
