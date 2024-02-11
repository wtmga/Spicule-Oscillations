%% MAIN

close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read video

% CHANGE "spic.mov" TO THE NAME OF YOUR VIDEO TO BE ANALYSED
videoName = 'spic.mov';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vid = VideoReader(videoName); %reads video and gives information
totalFrames = vid.FrameRate()*vid.Duration;

vid.CurrentTime = 0;  %resets video to beginning
vidFrame = vid.readFrame();  %defines vidFrames to use read frame function


%BEGINNING OF PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 0 SETUP

numFrames = 25; %set frames to be read, numFrames = 1 for testing

framesPerSecond = 1/16;
realFrameWidth = 100*10^6;

newCrop = 0; %0 uses saved crop, 1 uses new crop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 1 ENHANCING VIDEO FRAMES (stores in temp1)

enhanceMethod = 3; %0 orginal, 1 imadjust, 2 histeq, 3 adapthishiteq (recommended)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 2 GRAY THRESHOLDING (stores in temp2)

threshMethod = 0; %0 no gray thresholding (recommended for edge detection technique), 1 manual, 2 automatic


%Method 1 parameter

grayThreshold = 110; %custom threshold value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 3 SPICULE DETECTION (stores in temp3, and temp3Doub for method 2)

detectionMethod = 2; %0 no detection, 1 edge detection, 2 curvilinear detection


%Method 1 parameters

edgeMethod = "canny"; %detect edges using sobel, prewitt, roberts, log, zerocross, canny, or approxcanny
manualEdge = 1; %0 automatic parameters (below not required), 1 manual parameters (adjust below)

edgeThresh = 0.35; %between 0 and 1, lower values increase sensitivity, [for all methods]
sobPrewDirect = "vertical"; %vertical, horizontal, or both,  [for sobel and prewitt methods]
logCanSigma = 2.5; %default sigma is 2, lower values increase separation, [for log and canny methods]


%Method 2 parameters of Step Filter

StepFilterGrammes = 25; %size, proportional to connectedness of spicules, around 35
StepFilterDbw = [1 2]; %vector of width of main wave, range widths to detect, around 1.5

lowLimit = -5.5;  %orientation, moving clockwise: horizontal (-6) --> vert (0) --> horizontal (6)
StepFilterStep = pi/24; %angle step
highLimit = 5;  %must be >lowLimit
StepFilterDTheta = [pi-highLimit*StepFilterStep:StepFilterStep:pi-lowLimit*StepFilterStep]; %angle vector, pi vertical, pi/2 horizontal,

StepFilterType = 0; %type of step filter in {0,1} : 1 step function, 0 Polynomial filter

BinThresh = 700; %decreasing allows more spicules

smallThresh = 80; %remove objects with less than n pixels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 4 TRACKING SPICULES

newCrop2 = 0;  %0 uses saved crop, 1 uses new crop

newTrackReg = 0;  %0 uses saved tracking region, 1 uses new tracking region
orientationLimit = 10;  %limit of chance of spicule orientation in degrees (used to track the correct spicule)

matchIterations = 5;  %higher iterations track spicules more accurately, lower iteration allow marginally faster processing
showIndices = 0;  %1 shows detected spicule indices on a figure per frame processed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 5 PLOTTING STATISTICS

plotStats = 1;

plotCentOsc = 1;
plotTopOsc = 1;
plotBottOsc = 1;

plotCentVelOsc = 1;
plotTopVelOsc = 1;
plotBottVelOsc = 1;

plotDimensions = 1;
plotOrientation = 1;

plotHistograms = 0;
binNumber = 15;

numMatches = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END OF PARAMETERS


samplingFreq = framesPerSecond;  %sampling frequency in hertz = frames per second

pixelsPerMeter = vid.width/realFrameWidth;
pixFac = 1/pixelsPerMeter;

frameFac = 1/framesPerSecond;



if newCrop == 1 || ~exist('rect') %prompts crop if no original crop or specified to use new crop
    fprintf('\n----------------\nCrop to desired area of interest')
    [~,rect] = imcrop(vidFrame);
end


clear spicDataCell
clear spicLengthData
clear spicWidthData
clear spicOrientData

ii = 1; %set index
spicInstances = 0;
totWidth = 0;
totLength = 0;
totOrient = 0;


while ii<=numFrames
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %converting frame to grayscale
    
    vidFrame2 = vid.readFrame();
    croppedFrame = imcrop(vidFrame2, rect); %apply original crop
    grayVidFrame = rgb2gray(croppedFrame); %covert to gray
    
    %testing
    
    if numFrames == 1
        figure;
        imshow(grayVidFrame);
        %         title('Cropped gray frame')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % contrast enhancement techniques
    
    if enhanceMethod ==0
        enhancedFrame = grayVidFrame; %original
    end
    
    if enhanceMethod ==1
        enhancedFrame = imadjust(grayVidFrame);
    end
    
    if enhanceMethod ==2
        enhancedFrame = histeq(grayVidFrame);
    end
    
    if enhanceMethod ==3
        enhancedFrame = adapthisteq(grayVidFrame);
    end
    
    temp1(ii,:,:,:) = enhancedFrame; %store frame in temp
    
    
    %testing
    
    if numFrames == 1
        figure;
        imshow(enhancedFrame)
        %         title('Enhanced frame')
    end
    
    RESULT = enhancedFrame;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % THRESHOLDING METHOD
    
    % BEGINNING OF METHOD 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optional manual thresholding
    
    if threshMethod ==1
        
        enhancedFrameBW = RESULT < grayThreshold;
        threshFrame = RESULT .* uint8(enhancedFrameBW);
        
        temp2(ii,:,:,:) = threshFrame; %store frame in temp
        
        %testing
        
        if numFrames == 1
            figure;
            imshow(threshFrame), title('Manually thresholded frame');
        end
        
        RESULT = threshFrame;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END OF METHOD 1
    
    
    
    % BEGINNING OF METHOD 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optional automatic thresholding
    
    if threshMethod ==2
        
        level = graythresh(RESULT);
        BW = imbinarize(RESULT,level);
        BW = ~BW;
        BW = 1-BW;
        BW = (BW == 0);
        threshFrame = RESULT .* uint8(BW);
        
        temp2(ii,:,:,:) = threshFrame; %store frame in temp
        
        %testing
        
        if numFrames == 1
            figure;
            imshow(threshFrame),title('Automatically thresholded frame');
        end
        
        RESULT = threshFrame;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END OF METHOD 2
    
    
    
    % Spicule detection
    
    % BEGINNING OF TECHNIQUE 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % edge detection
    
    if detectionMethod == 1
        
        %automatic edge detection
        
        if manualEdge == 0
            bw0 = edge(RESULT, edgeMethod);
        end
        
        %manual parameters edge detection
        
        if manualEdge == 1
            
            if edgeMethod == "sobel" || edgeMethod == "prewitt"
                bw0 = edge(RESULT, edgeMethod, edgeThresh, sobPrewDirect);
            end
            
            if edgeMethod == "roberts" || edgeMethod == "zerocross" || edgeMethod == "approxcanny"
                bw0 = edge(RESULT, edgeMethod, edgeThresh);
            end
            
            if edgeMethod == "log" || edgeMethod == "canny"
                bw0 = edge(RESULT, edgeMethod, edgeThresh, logCanSigma);
            end
            
        end
        
        linesDetected = imfill(RESULT .* uint8(bw0)); %fill enclosed edges
        
        temp3(ii,:,:,:) = linesDetected; %store frame in temp
        
        
        %testing
        
        if numFrames == 1
            rgb0 = imoverlay_old(RESULT, bw0, [0 1 0]); %green overlay
            figure;
            imshow(rgb0),title('Detected edges');
            
            figure;
            imshow(linesDetected),title('Filled and detected edges');
            
            [xx, yy]=find(linesDetected);
            figure;
            plot(xx,yy, '.');
            title('Coordinates of detected lines');
            xlabel('X Coordinate')
            ylabel('Y coordinate')
        end
        
        RESULT = linesDetected;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END OF TECHNIQUE 1
    
    
    
    % BEGINNING OF TECHNIQUE 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if detectionMethod == 2
        
        
        % Finding curves in frame
        
        grammes = StepFilterGrammes;
        Dbw = StepFilterDbw;
        DTheta = StepFilterDTheta;
        
        [linesDetected]= step_pol_Filtering(RESULT,grammes,DTheta,Dbw,StepFilterType);
        
        temp3Doub(ii,:,:,:) = linesDetected; %store double frame in tempDoub
        
        %testing
        
        if numFrames == 1
            figure;
            imagesc(squeeze(linesDetected))
            set(gca,'visible','off')
            set(gca,'xtick',[])
            %             title('Detected curves')
            colormap('gray');
        end
        
        RESULT = linesDetected;
        
        
        
        % binarization
        
        RESULT(RESULT<=BinThresh)=0; %set values less than or equal to thresh to 0
        RESULT(RESULT>BinThresh)=1;
        
        %testing
        
        if numFrames == 1
            figure;
            imshow(RESULT)
            title('Binarized detected curves')
        end
        
        
        %removing small areas
        
        RESULT = bwareaopen(RESULT,smallThresh);
        
        temp3(ii,:,:,:) = RESULT; %store binary frame in temp
        
        
        %testing small patch removal
        
        if numFrames == 1
            figure;
            imshow(RESULT)
            %             title('Small areas removed')
            
            %plotting
            
            %             F = sparse(RESULT);
            %             figure;
            %             spy(F);  %plots sparsity pattern
            %             xlabel(' ')
            %             ylabel(' ')
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END OF TECHNIQUE 2
    
    
    
    % Spicule analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %new crop
    if (newCrop2 == 1 && ii==1) || ~exist('rect2') %prompts crop if no original crop or user specifies use of new crop on first iteration
        fprintf('\nCrop to desired area of analysis')
        [~,rect2] = imcrop(RESULT);
    end
    
    secondCropArea = rect2(3)*rect(4);
    
    %apply crop to enhancedFrame
    
    newEnhancedFrame = imcrop(enhancedFrame,rect2);
    
    if numFrames == 1
        figure
        imshow(newEnhancedFrame)
    end
    
    
    
    %apply crop to binary frame
    
    RESULT = imcrop(RESULT, rect2);  %apply second crop
    
    if numFrames == 1
        figure
        imshow(RESULT)
        figure
        imshow(uint8(RESULT).*newEnhancedFrame)  %show detected spicules from enhanced frame
    end
    
    
    
    
    
    clear spiculeInfo
    spiculeInfo = regionprops(RESULT, 'Centroid', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');  %find object centroids
    
    
    
    if numFrames == 1
        %display spicule centroids
        figure;
        imshow(RESULT)
        %         title('Spicule centroids')
        hold on
        for k = 1:numel(spiculeInfo)
            cent = spiculeInfo(k).Centroid;
            %             scatter(cent(1),cent(2),'r','filled') % plot these coordinates as a scatter, on top of the image, in red
            text(cent(1), cent(2), sprintf('%d', k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', '[1,0,0]');
        end
        %         hold off
        
        %display box bounds
        %         figure;
        %         imshow(RESULT),title('Spicule bounding boxes');
        %         hold on;
        for k = 1:numel(spiculeInfo)
            box = spiculeInfo(k).BoundingBox;
            rectangle('Position', box, 'EdgeColor', 'r');
        end
        hold off;
    end
    
    
    
    
    %save centroid, box bounds coordinates, length, width
    for k = 1:numel(spiculeInfo)
        
        %centroid
        cent = spiculeInfo(k).Centroid;
        spiculeData(k,1) = cent(1); %x coordinate
        spiculeData(k,2) = cent(2); %y
        
        %bounding box
        box = spiculeInfo(k).BoundingBox;
        spiculeData(k,3) = box(1);  %x of upper left
        spiculeData(k,4) = box(2);  %y of upper left
        spiculeData(k,5) = box(3); %width of box around object
        spiculeData(k,6) = box(4); %height of box around object
        
        %length and width
        spiculeData(k,7) = spiculeInfo(k).MajorAxisLength;
        spiculeData(k,8) = spiculeInfo(k).MinorAxisLength;
        
        %orientation
        spiculeData(k,9) = spiculeInfo(k).Orientation;
    end
    
    
    spicDataCell{ii}=spiculeData;
    %each cell contains rows (each a spicule) with: x,y of centroids, x,y of upper left bounding box, width and height of box, length, width, orientation
    
    
    
    %set threshold region for tracking spicule COG
    
    if (newTrackReg == 1 && ii==1) || ~exist('trackReg') %prompts crop if no original crop or user specifies use of new crop on first iteration
        fprintf('\nSelect acceptable tracking region (applied to all spicules)')
        [~,trackReg] = imcrop(newEnhancedFrame);
        
        close all
    end
    
    
    
    
    
    %Tracking spicules
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ii>1
        
        
        %deleting duplicates
        
        [prevNumRows,~] = size(spicDataCell{ii-1}); %get number of rows of previous frame spicules
        [currentNumRows,~] = size(spicDataCell{ii});
        
        clear tempNaNRemoval
        index=1;
        for j = 1:currentNumRows
            %if all fields of a current iteration row match some row from previous iteration --> current row becomes NaN
            for k = 1:prevNumRows
                if spicDataCell{ii}(j,1) == spicDataCell{ii-1}(k,1)...
                        && spicDataCell{ii}(j,2) == spicDataCell{ii-1}(k,2)...
                        && spicDataCell{ii}(j,3) == spicDataCell{ii-1}(k,3)...
                        && spicDataCell{ii}(j,4) == spicDataCell{ii-1}(k,4)...
                        && spicDataCell{ii}(j,5) == spicDataCell{ii-1}(k,5)...
                        && spicDataCell{ii}(j,6) == spicDataCell{ii-1}(k,6)...
                        && spicDataCell{ii}(j,7) == spicDataCell{ii-1}(k,7)...
                        && spicDataCell{ii}(j,8) == spicDataCell{ii-1}(k,8)...
                        && spicDataCell{ii}(j,9) == spicDataCell{ii-1}(k,9)...
                        
                    spicDataCell{ii}(j,:) = NaN; %fields set to zero if all fields match previous frame
                end
            end
            if spicDataCell{ii}(j,1)>0
                tempNaNRemoval(index,:) = spicDataCell{ii}(j,:);
                index=index+1;
            end
        end
        
        spicDataCell{ii}= tempNaNRemoval;  %new data with removed duplicates (from prev iteration)
        
        
        %redefine sizes of current and previous iterations
        
        [prevNumRows, ~] = size(spicDataCell{ii-1}); %get number of rows of previous frame spicules
        [currentNumRows, ~] = size(spicDataCell{ii});
        
        
        %allow centroids within defined tracking region
        
        clear possMat
        index=1;
        for k = 1:prevNumRows
            for j = 1:currentNumRows
                if (spicDataCell{ii}(j,1)>(spicDataCell{ii-1}(k,1)-trackReg(3)/2)...         %current x > previous x - width/2 &
                        &&spicDataCell{ii}(j,1)<(spicDataCell{ii-1}(k,1)+trackReg(3)/2))...  %current x < previous x + width/2
                        &&(spicDataCell{ii}(j,2)>(spicDataCell{ii-1}(k,2)-trackReg(4)/2)...  %current y > previous y - height/2 &
                        &&spicDataCell{ii}(j,2)<(spicDataCell{ii-1}(k,2)+trackReg(4)/2))...  %current y < previous y + height/2
                        &&spicDataCell{ii}(j,9)<spicDataCell{ii-1}(k,9)+orientationLimit...  %current orientation below previous orientation + upper limit
                        &&spicDataCell{ii}(j,9)>spicDataCell{ii-1}(k,9)-orientationLimit...  %current orientation above previous orientation - lower limit
                        
                    possMat(index, 1)=k;
                    possMat(index, 2)=j;  %save previous array index and current array possible match
                    index=index+1;
                end
            end
        end
        
        
        
        %find the distance between the potential matches
        
        [possMatRows, ~] = size(possMat);
        for k=1:possMatRows
            x1=spicDataCell{ii-1}(possMat(k,1),1);
            y1=spicDataCell{ii-1}(possMat(k,1),2);
            
            x2=spicDataCell{ii}(possMat(k,2),1);
            y2=spicDataCell{ii}(possMat(k,2),2);
            
            possMat(k,3)=sqrt((x2-x1)^2+(y2-y1)^2);
        end
        
        
        
        
        
        finalMatchMat=possMat;
        
        for gg=1:matchIterations
            
            %write to matrix smallest matches for each previous row (current rows can match with more than one previous row, NaN exist as minDist)
            
            [finalMatchMatRows, ~] = size(finalMatchMat);
            index=1;
            clear match
            
            for num=1:finalMatchMatRows  %loop for minimum of match and maximum matrix
                minDist = 2*(trackReg(3)+trackReg(4));  %arbitrary comparison, =perimeter of tracking region box
                couple = NaN;
                for k=1:finalMatchMatRows
                    if finalMatchMat(k,1) == num && finalMatchMat(k,3)<minDist  %find least distance for given value
                        minDist = finalMatchMat(k,3);
                        couple = finalMatchMat(k,2);
                    end
                end
                match(index,1) = num;  %record data
                match(index,2) = couple;
                match(index,3) = minDist;
                index=index+1;
                
            end
            
            
            %label NaN
            
            clear match2
            [matchRows, ~] = size(match);
            index=1;
            for k=1:matchRows
                if match(k,3)< 2*(trackReg(3)+trackReg(4))
                    match2(index,:) = match(k,:);
                else
                    match2(index,1:3)=NaN;
                end
                index=index+1;
            end
            
            
            
            %match current with only one (previous row), (NaN exist as minDist)
            
            clear match3
            index=1;
            [match2Rows, ~] = size(match2);
            for num=1:match2Rows  %loop for minimum of match and maximum matrix
                minDist = 2*(trackReg(3)+trackReg(4));  %arbitrary comparison, =perimeter of tracking region box
                couple = NaN;
                for k=1:match2Rows
                    if match2(k,2) == num && match2(k,3)<minDist && match2(k,3)>0  %find least distance for given value
                        minDist = match2(k,3);
                        couple = match2(k,1);
                    end
                end
                match3(index,2) = num;  %record data
                match3(index,1) = couple;
                match3(index,3) = minDist;
                index=index+1;
            end
            
            
            %label NaN
            
            clear match4
            [match3Rows, ~] = size(match3);
            index=1;
            for k=1:match3Rows
                if match3(k,3)< 2*(trackReg(3)+trackReg(4))
                    match4(index,:) = match3(k,:);
                else
                    match4(index,1:3) = NaN;
                end
                index=index+1;
            end
            
            
            
            %reintroduce previously eliminated rows (whose smallest distance corresponding row had a smaller distance corresponding row)
            
            clear match5
            index=1;
            [match4Rows, ~] = size(match4);
            
            for k=1:possMatRows
                exists=0;
                for j=1:match4Rows
                    if possMat(k,1)== match4(j,1) || possMat(k,2)== match4(j,2)  %if either the previous or current index exists
                        exists=1;
                        j = match4Rows;  %exit
                    end
                end
                if exists == 0 && possMat(k,3)>0
                    match5(index,:)=possMat(k,:);
                elseif exists == 1
                    match5(index,1:3)=NaN;
                end
                index=index+1;
                
            end
            
            clear finalMatchMat
            finalMatchMat=[match4; match5];  %stack match4 on top of match5 to be thrown back for sizing
            
        end
        
        
        
        %match indices of current iteration with previous iteration to track spicules
        
        [prevNumRows, ~] = size(spicDataCell{ii-1});
        [currentNumRows, currentNumCols] = size(spicDataCell{ii});
        
        clear tempSpicData
        for k=1:currentNumRows
            tempSpicData(prevNumRows + k,:)=spicDataCell{ii}(k,:);  %writes current frame data after its row size
        end
        
        index=currentNumRows+1;
        
        for k=1:prevNumRows
            matchExists=0;
            for j=1:match4Rows
                if k == match4(j,1)  %checks if previous frame row has match with current frame row
                    matchExists=1;
                    matchCurrentIndex = match4(j,2);  %saves matching current row index
                    j=match4Rows;  %exits loop
                end
            end
            
            if matchExists==0
                tempSpicData(k,1:currentNumCols) = NaN;  %labels current frame row as not a number if corresponding previous frame row has no match (to indicate end of lifetime)
            elseif matchExists==1
                tempSpicData(k,:) = spicDataCell{ii}(matchCurrentIndex,:);  %sets current frame row (with match from corresponding previous frame row) as the matching value
                tempSpicData(prevNumRows + matchCurrentIndex,:) = 999.999;  %flag moved rows
            end
        end
        
        [tempRows, ~] = size(tempSpicData);
        index=1;
        clear spicDataCell{ii}
        
        for k=1:tempRows
            if tempSpicData(k,1)~=999.999  %only writes data that has not been flagged as moved
                spicDataCell{ii}(index,1:currentNumCols) = tempSpicData(k,1:currentNumCols);
                index=index+1;
            end
        end
        
        
        
    end  %close if ii>1
    
    
    
    % write length, width, and oritentation data to separate arrays
    
    [currentNumRows, ~] = size(spicDataCell{ii});
    index = 1;
    for k=1:currentNumRows
        if spicDataCell{ii}(k,1) > 0
            spicInstances = spicInstances+1;
            
            spicLengthData(index,1) = spicDataCell{ii}(k,7);
            spicWidthData(index,1) = spicDataCell{ii}(k,8);
            spicOrientData(index,1) = spicDataCell{ii}(k,9);
            
            index = index + 1;
        end
    end
    
    
    
    %display new spicule indices (based on previous frames)
    if showIndices ==1
        
        figure
        imshow(RESULT),title('Spicule indices');
        hold on  %allows all spicules to be numbered throughout iterations
        
        [currentNumRows, ~] = size(spicDataCell{ii});
        
        for k = 1:currentNumRows
            
            %write indices of spicules
            if spicDataCell{ii}(k,1)>0
                text(spicDataCell{ii}(k,1), spicDataCell{ii}(k,2), sprintf('%d', k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', 'r');
            end
        end
        
        hold off
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ii = ii + 1; %move to next frame
    
    
    %frames processed counter
    
    if ii==2
        fprintf('\n-----------\nFrames processed: ');
    end
    
    if ii>1
        for j=0:log10(ii-2)
            fprintf('\b'); % delete previous counter display
        end
        fprintf('%d', ii-1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ii-1==numFrames
        
        
        while plotStats ==1 && numFrames>1
            fprintf('\n-----------\nEnter spicule index to be plotted (1 to %d):',currentNumRows);
            specificSpicPlot = input(' ');  %get user input for spicule index to be presented
            
            
            
            %finding max and min x and y axes to uniformly plot chosen spicule
            
            appearedBefore = 0;
            
            
            for k = 1:numFrames
                [currentNumRows,~] = size(spicDataCell{k});
                if (specificSpicPlot <= currentNumRows) && (spicDataCell{k}(specificSpicPlot,1)>0)
                    
                    %set initial top and bottom limits of x and y axes
                    if appearedBefore == 0
                        minX = spicDataCell{k}(specificSpicPlot,3);  %left boundbox x
                        minY = spicDataCell{k}(specificSpicPlot,4);  %top bounbox y
                        
                        maxX = spicDataCell{k}(specificSpicPlot,3)+spicDataCell{k}(specificSpicPlot,5);  %left boundbox x + boundbox width
                        maxY = spicDataCell{k}(specificSpicPlot,4)+spicDataCell{k}(specificSpicPlot,6);  %top boundbox y + boundbox height
                        
                        appearedBefore = 1;
                        
                    elseif appearedBefore ==1
                        
                        %compare min X
                        if spicDataCell{k}(specificSpicPlot,3) < minX
                            minX = spicDataCell{k}(specificSpicPlot,3);
                        end
                        
                        %compare min Y
                        if spicDataCell{k}(specificSpicPlot,4) < minY
                            minY = spicDataCell{k}(specificSpicPlot,4);
                        end
                        
                        %compare max X
                        if (spicDataCell{k}(specificSpicPlot,3)+spicDataCell{k}(specificSpicPlot,5)) > maxX
                            maxX = spicDataCell{k}(specificSpicPlot,3)+spicDataCell{k}(specificSpicPlot,5);
                            
                        end
                        
                        %compare max Y
                        if (spicDataCell{k}(specificSpicPlot,4)+spicDataCell{k}(specificSpicPlot,6)) > maxY
                            maxY = spicDataCell{k}(specificSpicPlot,4)+spicDataCell{k}(specificSpicPlot,6);
                        end
                        
                    end
                end
            end
            
            %convert pixels to meters
            
            minX = minX * pixFac;
            minY = minY * pixFac;
            
            maxX = maxX * pixFac;
            maxY = maxY * pixFac;
            
            
            
            %plotting spicule tops, centroids, and bottoms with previously found overal max and min x and y axes limits
            
            
            for k = 1:numFrames
                [currentNumRows,~] = size(spicDataCell{k});
                if (specificSpicPlot <= currentNumRows) && (spicDataCell{k}(specificSpicPlot,1)>0)
                    
                    figure
                    hold on
                    
                    grid on
                    title(sprintf('Plot of top, centroid, and bottom of spicule of index %d in frame %d\n', specificSpicPlot,k))
                    set(gca, 'YDir', 'reverse')
                    
                    xlabel('x-coordinate in m')
                    ylabel('y-coordinate in m')
                    
                    %using obtained min and maxes to set plot limits
                    xlim([minX, maxX]);
                    ylim([minY, maxY]);
                    
                    
                    %plot centroid
                    
                    xCent = spicDataCell{k}(specificSpicPlot,1);
                    yCent = spicDataCell{k}(specificSpicPlot,2);
                    
                    %convert to meters and plot
                    scatter(xCent*pixFac,yCent*pixFac,'filled','r')  %plot the centroid as a point on a scatter plot
                    
                    
                    if spicDataCell{k}(specificSpicPlot,9)>0  %positive orientation
                        xTop = spicDataCell{k}(specificSpicPlot,3)+spicDataCell{k}(specificSpicPlot,5);  %top left x + width of boundbox
                        yTop = spicDataCell{k}(specificSpicPlot,4);  %top left y
                        
                        xBott = spicDataCell{k}(specificSpicPlot,3);  %top left x
                        yBott = spicDataCell{k}(specificSpicPlot,4)+spicDataCell{k}(specificSpicPlot,6);  %top left y + height of boundbox
                        
                        %convert to meters
                        scatter(xTop*pixFac,yTop*pixFac,'filled','g')  %plot the top as a point on a scatter plot
                        scatter(xBott*pixFac,yBott*pixFac,'filled','b')  %plot the bottom as a point on a scatter plot
                        
                    elseif spicDataCell{k}(specificSpicPlot,9)<0  %negative orientation
                        xTop = spicDataCell{k}(specificSpicPlot,3);  %top left x
                        yTop = spicDataCell{k}(specificSpicPlot,4);  %top left y
                        
                        xBott = spicDataCell{k}(specificSpicPlot,3)+spicDataCell{k}(specificSpicPlot,5);  %top left x + width of boundbox
                        yBott = spicDataCell{k}(specificSpicPlot,4)+spicDataCell{k}(specificSpicPlot,6);  %top left y + height of boundbox
                        
                        %convert to meters
                        scatter(xTop*pixFac,yTop*pixFac,'filled','g')  %plot the top as a point on a scatter plot
                        scatter(xBott*pixFac,yBott*pixFac,'filled','b')  %plot the bottom as a point on a scatter plot
                    end
                    
                    
                    legend({'Centroid','Top','Bottom'})
                    hold off
                    
                end
                
            end
            
            
            
            
            
            
            
            
            
            
            
            %variables that store oscillation frequencies of all spicules
            
            clear spicCentFreq
            clear spicTopFreq
            clear spicBottFreq
            
            clear spicCentVelFreq
            clear spicTopVelFreq
            clear spicBottVelFreq
            
            clear allSpicCentVel
            clear allSpicTopVel
            clear allSpicBottVel
            
            clear spicLengthFreq
            clear spicWidthFreq
            
            clear spicOrientFreq
            
            
            
            %go through all detected spicules
            
            [totNumSpic,~] = size(spicDataCell{ii-1});
            
            for currentSpicAnalysis = 1:totNumSpic
                
                index=1;
                
                %clearing variables that store displacement from average centroid, top, bottom, respectively
                
                clear spicCentOsc
                clear spicTopOsc
                clear spicBottOsc
                
                appearedBefore = 0;
                
                sumCentXPos = 0;
                sumCentYPos = 0;
                
                sumTopXPos = 0;
                sumTopYPos = 0;
                
                sumBottXPos = 0;
                sumBottYPos = 0;
                
                counter = 0;
                
                
                if (plotCentOsc ==1 || plotTopOsc ==1 || plotBottOsc ==1) && currentSpicAnalysis == specificSpicPlot
                    figure
                    hold on
                    grid on
                    grid minor
                end
                
                
                for k=1:numFrames
                    [currentNumRows,~] = size(spicDataCell{k});
                    
                    if currentSpicAnalysis <= currentNumRows && spicDataCell{k}(currentSpicAnalysis,1) > 0  %if the chosen spicule to be analyzed is greater than number of rows then skip cell since spicule hasn't appeared yet
                        
                        %get average centroid, top, bottom, coordinates to set as reference to measure oscillation
                        
                        if appearedBefore == 0
                            
                            %start loop from first frame spicule appears until last frame
                            for gg=k:numFrames
                                
                                if spicDataCell{gg}(currentSpicAnalysis,1) > 0
                                    
                                    sumCentXPos = sumCentXPos + spicDataCell{gg}(currentSpicAnalysis,1);
                                    sumCentYPos = sumCentYPos + spicDataCell{gg}(currentSpicAnalysis,2);
                                    
                                    if spicDataCell{gg}(currentSpicAnalysis,9)>0  %positive orientation
                                        sumTopXPos = sumTopXPos + spicDataCell{gg}(currentSpicAnalysis,3)+spicDataCell{gg}(currentSpicAnalysis,5);  %top left x + width of boundbox
                                        sumTopYPos = sumTopYPos + spicDataCell{gg}(currentSpicAnalysis,4);  %top left y
                                        
                                        sumBottXPos = sumBottXPos + spicDataCell{gg}(currentSpicAnalysis,3);  %top left x
                                        sumBottYPos = sumBottYPos + spicDataCell{gg}(currentSpicAnalysis,4)+spicDataCell{gg}(currentSpicAnalysis,6);  %top left y + height of boundbox
                                        
                                    elseif spicDataCell{gg}(currentSpicAnalysis,9)<0  %negative orientation
                                        sumTopXPos = sumTopXPos + spicDataCell{gg}(currentSpicAnalysis,3);
                                        sumTopYPos = sumTopYPos + spicDataCell{gg}(currentSpicAnalysis,4);
                                        
                                        sumBottXPos = sumBottXPos + spicDataCell{gg}(currentSpicAnalysis,3)+spicDataCell{gg}(currentSpicAnalysis,5);  %top left x + width of boundbox
                                        sumBottYPos = sumBottYPos + spicDataCell{gg}(currentSpicAnalysis,4)+spicDataCell{gg}(currentSpicAnalysis,6);  %top left y + height of boundbox
                                        
                                    end
                                    
                                    counter = counter + 1;
                                end
                            end
                            
                            %setting reference points (average centroid, top, bottom for spicule)
                            
                            
                            
                            x1 = sumCentXPos/counter;
                            y1 = sumCentYPos/counter;
                            
                            x1top = sumTopXPos/counter;
                            y1top = sumTopYPos/counter;
                            
                            x1bott = sumBottXPos/counter;
                            y1bott = sumBottYPos/counter;
                            
                            
                            appearedBefore = 1;  %mark spicule as appeared prior to avoid repeating process above
                            
                            
                        end
                        
                        
                    if appearedBefore == 1
                            
                            %get the x and y centroid coordinates of the current frame
                            
                            x2 = spicDataCell{k}(currentSpicAnalysis,1);
                            y2 = spicDataCell{k}(currentSpicAnalysis,2);
                            
                            spicDisp = sqrt((x2-x1)^2+(y2-y1)^2);  %find the distance between the current and average centroid of that spicule
                            
                            if x2>x1 || (x1==x2 && y2>y1)  %moved horizontally right or no horizontal movement but vertically upward
                                spicDeltaDisp = spicDisp;  %keep the distance change positive if the current x coordinate is to the right of the previous x coordinate
                            elseif x2<x1 || (x1==x2 && y2<y1)
                                spicDeltaDisp = -1*spicDisp;  %make the distance change negative if the current x coordinate is to the left of the previous x coordinate
                            end
                            
                            %convert to meters
                            spicCentOsc(index,1) = spicDeltaDisp*pixFac;  %store the current displacement from the mean centroid
                            
                            
                            
                            %calculating top and bottom spicule tips displacement
                            
                            if spicDataCell{k}(currentSpicAnalysis,9)>0  %if positive orientation (between neutral and clockwise up to horizontal)
                                %top coordinates
                                x2top = spicDataCell{k}(currentSpicAnalysis,3)+spicDataCell{k}(currentSpicAnalysis,5);
                                y2top = spicDataCell{k}(currentSpicAnalysis,4);
                                
                                %bottom coordinates
                                x2bott = spicDataCell{k}(currentSpicAnalysis,3);
                                y2bott = spicDataCell{k}(currentSpicAnalysis,4)+spicDataCell{k}(currentSpicAnalysis,6);
                                
                            elseif spicDataCell{k}(currentSpicAnalysis,9)<0
                                %top coordinates
                                x2top = spicDataCell{k}(currentSpicAnalysis,3);
                                y2top = spicDataCell{k}(currentSpicAnalysis,4);
                                
                                %bottom coordinates
                                x2bott = spicDataCell{k}(currentSpicAnalysis,3)+spicDataCell{k}(currentSpicAnalysis,5);
                                y2bott = spicDataCell{k}(currentSpicAnalysis,4)+spicDataCell{k}(currentSpicAnalysis,6);
                            end
                            
                            
                            %top coordinate
                            spicDisp = sqrt((x2top-x1top)^2+(y2top-y1top)^2);  %find the distance between the current and previous frame centroids of that spicule
                            
                            if x2top>x1top || (x1top==x2top && y2top>y1top)
                                spicDeltaDisp = spicDisp;  %keep the distance change positive if the current x coordinate is to the right of the previous x coordinate
                            elseif x2top<x1top || (x1top==x2top && y2top<y1top)
                                spicDeltaDisp = -1*spicDisp;  %make the distance change negative if the current x coordinate is to the left of the previous x coordinate
                            end
                            
                            %convert to meters
                            spicTopOsc(index,1) = spicDeltaDisp*pixFac;
                            
                            
                            
                            %bottom coordinate
                            
                            spicDisp = sqrt((x2bott-x1bott)^2+(y2bott-y1bott)^2);  %find the distance between the current and previous frame centroids of that spicule
                            
                            if x2bott>x1bott || (x1bott==x2bott && y2bott>y1bott)
                                spicDeltaDisp = spicDisp;  %keep the distance change positive if the current x coordinate is to the right of the previous x coordinate
                            elseif x2bott<x1bott || (x1bott==x2bott && y2bott<y1bott)
                                spicDeltaDisp = -1*spicDisp;  %make the distance change negative if the current x coordinate is to the left of the previous x coordinate
                            end
                            
                            %convert to meters
                            spicBottOsc(index,1) = spicDeltaDisp*pixFac;
                            
                            
                            
                            %plotting
                            
                            if plotCentOsc == 1 && currentSpicAnalysis == specificSpicPlot
                                scatter(k,spicCentOsc(index,1),'filled','r')  %plot the current displacement as a point on a scatter plot
                            end
                            
                            if plotTopOsc == 1 && currentSpicAnalysis == specificSpicPlot
                                scatter(k,spicTopOsc(index,1),'filled','g')
                            end
                            
                            if plotBottOsc == 1 && currentSpicAnalysis == specificSpicPlot
                                scatter(k,spicBottOsc(index,1),'filled','b')
                            end
                            
                        
                        end
                        
                        index = index+1;
                        
                        
                        
                    end
                    
                end
                
                
                if index-1>1
                    %centroid oscillation
                    
                    dataLength = length(spicCentOsc);
                    oscFourier = fft(spicCentOsc);
                    
                    oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                    [~,actualFreq] = max(abs(oscFourier));
                    freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                    
                    spicCentFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                    
                    
                    %top oscillation
                    
                    dataLength = length(spicTopOsc);
                    oscFourier = fft(spicTopOsc);
                    
                    oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                    [~,actualFreq] = max(abs(oscFourier));
                    freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                    
                    spicTopFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                    
                    
                    %bottom oscillation
                    
                    dataLength = length(spicBottOsc);
                    oscFourier = fft(spicBottOsc);
                    
                    oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                    [~,actualFreq] = max(abs(oscFourier));
                    freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                    
                    spicBottFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                    
                    
                    %plot specific spicule properties
                    
                    if (plotCentOsc ==1 || plotTopOsc ==1 || plotBottOsc ==1) && currentSpicAnalysis == specificSpicPlot
                        title(sprintf('Displacement vs frame number plot of spicule of index %d\nLifetime = %d s\nCentroid oscillation frequency = %d Hz\nTop oscillation frequency = %d Hz\nBottom oscillation frequency = %d Hz', specificSpicPlot, (index-1)*frameFac, spicCentFreq(currentSpicAnalysis,1), spicTopFreq(currentSpicAnalysis,1), spicBottFreq(currentSpicAnalysis,1)))  %finds most likely frequency
                        xlabel('Frame number')
                        ylabel('Displacement from mean location in m')
                        legend({'Centroid','Top','Bottom'}, 'location', 'northwest')
                        set(gca, 'YDir', 'reverse')
                        hold off
                    end
                    
                elseif (plotCentOsc ==1 || plotTopOsc ==1 || plotBottOsc ==1) && currentSpicAnalysis == specificSpicPlot
                    hold off
                    close
                end
                
                
                
                
                
                
                
                
                %storing all spicule velocity information and plotting selected spicule velocity
                
                index=1;
                
                %clearing variables that store displacement from previous frame centroid, top, bottom, respectively
                
                clear spicCentVelOsc
                clear spicTopVelOsc
                clear spicBottVelOsc
                
                appearedBefore = 0;
                
                
                if (plotCentVelOsc ==1 || plotTopVelOsc ==1 || plotBottVelOsc ==1) && currentSpicAnalysis == specificSpicPlot
                    figure
                    hold on
                    grid on
                    grid minor
                end
                
                
                
                for k=1:numFrames
                    [currentNumRows,~] = size(spicDataCell{k});
                    
                    if currentSpicAnalysis <= currentNumRows && spicDataCell{k}(currentSpicAnalysis,1) > 0  %if the chosen spicule to be analyzed is greater than number of rows then skip cell since spicule hasn't appeared yet
                        
                        %get average centroid, top, bottom, coordinates to set as reference to measure oscillation
                        
                        if appearedBefore == 0
                            
                            
                            %mark first position as no velocity as it has just appeared
                            
                            spicCentVelOsc(index,1) = 0;
                            spicTopVelOsc(index,1) = 0;
                            spicBottVelOsc(index,1) = 0;
                            
                            
                            appearedBefore = 1;  %mark spicule as appeared prior to avoid repeating process above
                            
                            
                        elseif appearedBefore == 1  %if existed before
                            
                            %get previous frame coordinates
                            
                            x1 = spicDataCell{k-1}(currentSpicAnalysis,1);
                            y1 = spicDataCell{k-1}(currentSpicAnalysis,2);
                            
                            
                            %get the x and y centroid coordinates of the current frame
                            
                            x2 = spicDataCell{k}(currentSpicAnalysis,1);
                            y2 = spicDataCell{k}(currentSpicAnalysis,2);
                            
                            distance = sqrt((x2-x1)^2+(y2-y1)^2);  %find the distance between the current and average centroid of that spicule
                            
                            %convert pixels to meters and frames to seconds
                            spicCentVelOsc(index,1) = distance*pixFac/frameFac;  %velocity = distance/time
                            
                            
                            
                            %getting x and y coordinates of top and bottom according to orientation of spicule in previous frame
                            
                            if spicDataCell{k-1}(currentSpicAnalysis,9)>0  %if positive orientation (between neutral and clockwise up to horizontal)
                                %top coordinates
                                x1top = spicDataCell{k-1}(currentSpicAnalysis,3)+spicDataCell{k-1}(currentSpicAnalysis,5);
                                y1top = spicDataCell{k-1}(currentSpicAnalysis,4);
                                
                                
                                %bottom coordinates
                                x1bott = spicDataCell{k-1}(currentSpicAnalysis,3);
                                y1bott = spicDataCell{k-1}(currentSpicAnalysis,4)+spicDataCell{k-1}(currentSpicAnalysis,6);
                                
                                
                            elseif spicDataCell{k-1}(currentSpicAnalysis,9)<0
                                %top coordinates
                                x1top = spicDataCell{k-1}(currentSpicAnalysis,3);
                                y1top = spicDataCell{k-1}(currentSpicAnalysis,4);
                                
                                %bottom coordinates
                                x1bott = spicDataCell{k-1}(currentSpicAnalysis,3)+spicDataCell{k-1}(currentSpicAnalysis,5);
                                y1bott = spicDataCell{k-1}(currentSpicAnalysis,4)+spicDataCell{k-1}(currentSpicAnalysis,6);
                                
                            end
                            
                            
                            
                            %getting x and y coordinates of top and bottom according to orientation of spicule in current frame
                            
                            if spicDataCell{k}(currentSpicAnalysis,9)>0  %if positive orientation (between neutral and clockwise up to horizontal)
                                %top coordinates
                                x2top = spicDataCell{k}(currentSpicAnalysis,3)+spicDataCell{k}(currentSpicAnalysis,5);
                                y2top = spicDataCell{k}(currentSpicAnalysis,4);
                                
                                %bottom coordinates
                                x2bott = spicDataCell{k}(currentSpicAnalysis,3);
                                y2bott = spicDataCell{k}(currentSpicAnalysis,4)+spicDataCell{k}(currentSpicAnalysis,6);
                                
                                
                            elseif spicDataCell{k}(currentSpicAnalysis,9)<0
                                %top coordinates
                                
                                x2top = spicDataCell{k}(currentSpicAnalysis,3);
                                y2top = spicDataCell{k}(currentSpicAnalysis,4);
                                
                                %bottom coordinates
                                
                                x2bott = spicDataCell{k}(currentSpicAnalysis,3)+spicDataCell{k}(currentSpicAnalysis,5);
                                y2bott = spicDataCell{k}(currentSpicAnalysis,4)+spicDataCell{k}(currentSpicAnalysis,6);
                            end
                            
                            
                            
                            %calculating top and bottom spicule tips displacement
                            
                            %top coordinate
                            distance = sqrt((x2top-x1top)^2+(y2top-y1top)^2);  %find the distance between the current and previous frame centroids of that spicule
                            spicTopVelOsc(index,1) = distance*pixFac/frameFac;
                            
                            
                            
                            %bottom coordinate
                            
                            distance = sqrt((x2bott-x1bott)^2+(y2bott-y1bott)^2);  %find the distance between the current and previous frame centroids of that spicule
                            spicBottVelOsc(index,1) = distance*pixFac/frameFac;
                            
                            
                            %plotting
                            
                            if plotCentVelOsc == 1 && currentSpicAnalysis == specificSpicPlot
                                scatter(k,spicCentVelOsc(index,1),'filled','r')  %plot the current displacement as a point on a scatter plot
                            end
                            
                            if plotTopVelOsc == 1 && currentSpicAnalysis == specificSpicPlot
                                scatter(k,spicTopVelOsc(index,1),'filled','g')
                            end
                            
                            if plotBottVelOsc == 1 && currentSpicAnalysis == specificSpicPlot
                                scatter(k,spicBottVelOsc(index,1),'filled','b')
                            end
                            
                            
                        end
                        
                        
                        index = index+1;
                        
                        
                    end
                    
                end
                
                %storing mean velocities of one spicule
                
                allSpicCentVel(currentSpicAnalysis,1) =  mean(spicCentVelOsc);
                allSpicTopVel(currentSpicAnalysis,1) =  mean(spicTopVelOsc);
                allSpicBottVel(currentSpicAnalysis,1) =  mean(spicBottVelOsc);
                
                
                if index-1>1
                    %centroid oscillation
                    
                    dataLength = length(spicCentVelOsc);
                    oscFourier = fft(spicCentVelOsc);
                    
                    oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                    [~,actualFreq] = max(abs(oscFourier));
                    freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                    
                    spicCentVelFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                    
                    
                    %top oscillation
                    
                    dataLength = length(spicTopVelOsc);
                    oscFourier = fft(spicTopVelOsc);
                    
                    oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                    [~,actualFreq] = max(abs(oscFourier));
                    freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                    
                    spicTopVelFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                    
                    
                    %bottom oscillation
                    
                    dataLength = length(spicBottVelOsc);
                    oscFourier = fft(spicBottVelOsc);
                    
                    oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                    [~,actualFreq] = max(abs(oscFourier));
                    freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                    
                    spicBottVelFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                    
                    
                    
                    %plot specific spicule properties
                    
                    if (plotCentVelOsc ==1 || plotTopVelOsc ==1 || plotBottVelOsc ==1) && currentSpicAnalysis == specificSpicPlot
                        title(sprintf('Velocity vs frame number plot of spicule of index %d\nMean centroid velocity = %d m/s\nMean top velocity = %d m/s\nMean bottom velocity = %d m/s\nCentroid velocity oscillation frequency = %d Hz\nTop velocity oscillation frequency = %d Hz\nBottom velocity oscillation frequency = %d Hz', specificSpicPlot, mean(spicCentVelOsc), mean(spicTopVelOsc), mean(spicBottVelOsc), spicCentVelFreq(currentSpicAnalysis,1), spicTopVelFreq(currentSpicAnalysis,1), spicBottVelFreq(currentSpicAnalysis,1)))  %finds most likely frequency
                        xlabel('Frame number')
                        ylabel('Velocity in m/s')
                        legend({'Centroid','Top','Bottom'}, 'location', 'north')
                        set(gca, 'YDir', 'reverse')
                        hold off
                    end
                    
                    
                elseif (plotCentVelOsc ==1 || plotTopVelOsc ==1 || plotBottVelOsc ==1) && currentSpicAnalysis == specificSpicPlot
                    
                    hold off
                    close
                    
                end
                
                
                
                
                
                %storing all spicule dimensions and plotting chosen spicule
                
                index = 1;
                clear spicLengthOsc
                clear spicWidthOsc
                
                if plotDimensions ==1 && currentSpicAnalysis == specificSpicPlot
                    figure;
                    hold on
                    grid on
                    grid minor
                    xlabel('Frame number')
                    ylabel('Dimension in m')
                end
                
                for k=1:numFrames
                    [currentNumRows,~] = size(spicDataCell{k});
                    
                    if currentSpicAnalysis <= currentNumRows && spicDataCell{k}(currentSpicAnalysis,1) > 0  %if the chosen spicule to be analyzed is greater than number of rows then skip cell since spicule hasn't appeared yet
                        
                        %convert to meters
                        spicLengthOsc(index,1)= spicDataCell{k}(currentSpicAnalysis,7)*pixFac;
                        spicWidthOsc(index,1)= spicDataCell{k}(currentSpicAnalysis,8)*pixFac;
                        
                        
                        if plotDimensions ==1 && currentSpicAnalysis == specificSpicPlot
                            scatter(k,spicLengthOsc(index,1),'filled','b')  %plot the current displacement as a point on a scatter plot
                            scatter(k,spicWidthOsc(index,1),'filled','r')
                        end
                        
                        index = index+1;
                    end
                end
                
                allSpicLength(currentSpicAnalysis,1) =  mean(spicLengthOsc);
                allSpicWidth(currentSpicAnalysis,1) =  mean(spicWidthOsc);
                
                
                
                %length oscillation
                
                dataLength = length(spicLengthOsc);
                oscFourier = fft(spicLengthOsc);
                
                oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                [~,actualFreq] = max(abs(oscFourier));
                freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                
                spicLengthFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                
                
                %width oscillation
                
                dataLength = length(spicWidthOsc);
                oscFourier = fft(spicWidthOsc);
                
                oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                [~,actualFreq] = max(abs(oscFourier));
                freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                
                spicWidthFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                
                
                %plot
                
                if plotDimensions ==1 && currentSpicAnalysis == specificSpicPlot
                    title(sprintf('Length and width vs frame plot of spicule of index %d\nMean length = %d m\nMean width = %d m\nSpicule length osccilation frequency = %d Hz\nSpicule width oscillation frequency = %d Hz', specificSpicPlot, mean(spicLengthOsc), mean(spicWidthOsc), spicLengthFreq(currentSpicAnalysis,1), spicWidthFreq(currentSpicAnalysis,1)));  %plot the current displacement as a point on a scatter plot
                    legend({'Length','Width'}, 'location', 'northwest')
                    hold off
                end
                
                
                
                
                
                
                
                
                %recording and plotting orienation data
                
                index = 1;
                clear spicOrientOsc
                
                if plotOrientation ==1 && currentSpicAnalysis == specificSpicPlot
                    figure;
                    hold on
                    grid on
                    grid minor
                    xlabel('Frame number')
                    ylabel('Orientation in degrees')
                end
                
                for k=1:numFrames
                    [currentNumRows, ~] = size(spicDataCell{k});
                    
                    if currentSpicAnalysis <= currentNumRows && spicDataCell{k}(currentSpicAnalysis,1) > 0  %if the chosen spicule to be analyzed is greater than number of rows then skip cell since spicule hasn't appeared yet
                        
                        spicOrientOsc(index,1)= spicDataCell{k}(currentSpicAnalysis,9);
                        
                        if plotOrientation ==1 && currentSpicAnalysis == specificSpicPlot
                            scatter(k,spicOrientOsc(index,1),'filled','r')
                        end
                        
                        index = index+1;
                    end
                end
                
                
                allSpicOrient(currentSpicAnalysis,1) =  mean(spicOrientOsc);
                
                
                
                dataLength = length(spicOrientOsc);
                oscFourier = fft(spicOrientOsc);
                
                oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                [~,actualFreq] = max(abs(oscFourier));
                freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                
                spicOrientFreq(currentSpicAnalysis,1) = freq(actualFreq);  %store oscillation frequency
                
                
                %plot
                if plotOrientation ==1 && currentSpicAnalysis == specificSpicPlot
                    title(sprintf('Orientation vs frame plot of spicule of index %d\nMean spicule orientation = %d degrees\nSpicule orientation oscillation frequency = %d Hz', specificSpicPlot, mean(spicOrientOsc), spicOrientFreq(currentSpicAnalysis,1)));  %plot the current displacement as a point on a scatter plot
                    hold off
                end
                
                
                
                
                
                
                
                
            end  %for currentSpicAnalysis = 1:totNumSpic
            
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %HISTOGRAMS
            
            
            if plotHistograms == 1
                
                
                %new spicules and new spicule density
                
                clear spicBirthData
                clear spicBirthDensity
                
                
                %set initial to 0 since no new births
                
                spicBirthData(1,1) = 0;
                spicBirthDensity(1,1) = 0;
                
                
                for k=2:numFrames
                    [prevNumRows, ~] = size(spicDataCell{k-1});
                    [currentNumRows, ~] = size(spicDataCell{k});
                    
                    spicBirthData(k,1)= (currentNumRows - prevNumRows)*framesPerSecond;  %new spicules per second
                    
                    spicBirthDensity(k,1) = spicBirthData(k,1)/(secondCropArea*pixFac^2);  %new spicules per second per meter square
                end
                
                
                %find oscillation
                
                dataLength = length(spicBirthDensity);
                oscFourier = fft(spicBirthDensity);
                
                oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                [~,actualFreq] = max(abs(oscFourier));
                freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                
                spicBirthDensityFreq = freq(actualFreq);  %store oscillation frequency
                
                
                
                
                
                
                %all active spicules and spicule density
                
                clear activeSpicules
                clear spicDensity
                
                for k=1:numFrames
                    spiculeCounter = 0;
                    [currentNumRows, ~] = size(spicDataCell{k});
                    
                    for j = 1:currentNumRows
                        if spicDataCell{k}(j,1)>0
                            spiculeCounter = spiculeCounter + 1;
                        end
                    end
                    
                    activeSpicules(k,1) = spiculeCounter;  %active spicules in a frsme
                    
                    spicDensity(k,1) = spiculeCounter/(secondCropArea*pixFac^2);  %spicules per meter squared
                end
                
                
                %find oscillation
                
                dataLength = length(spicDensity);
                oscFourier = fft(spicDensity);
                
                oscFourier(floor(1:dataLength/2+1));  %only need to search 1/2 of oscFourier for the max because x is real-valued
                [~,actualFreq] = max(abs(oscFourier));
                freq = 0:samplingFreq/dataLength:samplingFreq/2;  %gives range of possible frequencies
                
                spicDensityFreq = freq(actualFreq);  %store oscillation frequency
                
                
                
                
                
                
                %spicule lifetime
                
                clear spicLifetimeData
                prevNumRows = 0;
                index=1;
                
                for gg = 1:numFrames  %cycle through all frames
                    if gg>1  %get previous frame number of rows (starting after first iteration)
                        [prevNumRows, ~] = size(spicDataCell{gg-1});
                    end
                    
                    [currentNumRows, ~] = size(spicDataCell{gg});  %get current number of rows
                    
                    for k=prevNumRows+1:currentNumRows  %cycle starting from first index after previous fram indices up till current greatest index
                        
                        aliveFrameCount = 0;  %reset 'alive' frame counter
                        
                        for j=gg:numFrames  %cycle starting from current frame til the end of frames
                            if spicDataCell{j}(k,1)>0
                                aliveFrameCount=aliveFrameCount+1;  %count as alive if not NaN
                            end
                        end
                        
                        spicLifetimeData(index,1) = aliveFrameCount*frameFac;  %convert frames to seconds
                        index=index+1;
                        
                    end
                    
                end
                
                
                
                
                
                
                
                %mean calculations (displayed on histograms)
                
                
                [currentNumRows, ~] = size(spicDataCell{ii-1});
                totNumSpic = currentNumRows;
                
                avLength = mean(allSpicLength);
                avWidth = mean(allSpicWidth);
                avOrient = mean(allSpicOrient);
                
                
                
                avLifetime = mean(spicLifetimeData);
                
                
                
                avSpicPerFrame = mean(activeSpicules);
                avDensity = mean(spicDensity);
                
                avNewSpic = mean(spicBirthData);
                avNewSpicDensity = mean(spicBirthDensity);
                
                
                
                avSpicCentVel = mean(allSpicCentVel);
                avSpicTopVel = mean(allSpicTopVel);
                avSpicBottVel = mean(allSpicBottVel);
                
                
                
                avOscFreq = mean(spicTopFreq);
                avTopOscFreq = mean(spicTopFreq);
                avBottOscFreq = mean(spicBottFreq);
                
                avCentVelFreq = mean(spicCentVelFreq);
                avTopVelFreq = mean(spicTopVelFreq);
                avBottVelFreq = mean(spicBottVelFreq);
                
                avOrientOscFreq = mean(spicOrientFreq);
                
                
                
                
                
                
                
                %plot length
                realLength = spicLengthData*pixFac;  %convert to meters
                minimum = min(realLength(realLength>0));
                maximum = max(realLength);
                
                data = realLength(realLength>0);
                
                figure; tiledlayout(1,2)
                nexttile
                fitmethis(data,'pdist',numMatches,'output','off');
                title(sprintf('Distribution fitting output (top 4 matches)\n')) 
                
                nexttile
                distribution = 'loglogistic';
                pd = fitdist(data,distribution);
                
                histfit(data,16,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.1f');
                grid on
                title(sprintf('Spicule length histogram\nMean spicule length = %d m\nDistribtion: %s\nmu = %d\nsigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
%                 title(sprintf('Spicule length histogram\nMean spicule lwngth = %d meters', avLength))
                xlabel('Length in m')
                ylabel('Probability')
                
                %plot length osc freq
                maximum = max(spicLengthFreq);
                if maximum>0
                    minimum = min(spicLengthFreq(spicLengthFreq>0));
                    figure
                    title('Spicule length oscillation frequency histogram')
                    xlabel('Length Oscillation Frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule length')
                end
                
                
                
                %plot width
                realWidth = spicWidthData*pixFac;
                minimum = min(realWidth(realWidth>0));
                maximum = max(realWidth);
                
                data = realWidth(realWidth>0);
                
                figure; tiledlayout(1,2)
                nexttile
                fitmethis(data,'pdist',numMatches,'output','off');
                title(sprintf('Distribution fitting output (top 4 matches)\n')) 
                
                nexttile
                distribution = 'loglogistic';
                pd = fitdist(data,distribution);
                
                histfit(data,23,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.1f');
                grid on
                title(sprintf('Spicule width histogram\nMean spicule width = %d m\nDistribtion: %s\nmu = %d\nsigma = %d', avWidth, distribution, pd.mu, pd.sigma))
                
%                 title(sprintf('Spicule width histogram\nMean spicule width = %d meters', avWidth))
                xlabel('Width in m')
                ylabel('Probability')
                
                %plot width osc freq
                maximum = max(spicWidthFreq);
                if maximum>0
                    minimum = min(spicWidthFreq(spicWidthFreq>0));
                    figure
                    histfit(spicWidthFreq(spicWidthFreq>0),binNumber,'normal');
                    title('Spicule width oscillation frequency histogram')
                    xlabel('Width oscillation Frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule width')
                end
                
                
                
                %plot orientation
                minimum = min(spicOrientData);
                maximum = max(spicOrientData);
                figure
                histogram(spicOrientData,[minimum:(maximum-minimum)/binNumber:maximum],'normalization','probability');
                title(sprintf('Spicule orientation histogram\nMean spicule orientation = %d degrees', avOrient))
                xlabel('Orientation in Degrees')
                ylabel('Probability')
                
                %plot orientation osc freq
                maximum = max(spicOrientFreq);
                if maximum>0
                    minimum = min(spicOrientFreq(spicOrientFreq>0));
                    figure
                    histfit(spicOrientFreq(spicOrientFreq>0),binNumber,'normal');
                    title(sprintf('Spicule orentation oscillation frequency histogram\nMean spicule orientation oscillation frequency = %d Hz', avOrientOscFreq))
                    xlabel('Orientation Oscillation Frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule orientation')
                end
                
                
                
                %plot centroid velocity
                minimum = min(allSpicCentVel(allSpicCentVel>0));
                maximum = max(allSpicCentVel);
                
                data = allSpicCentVel(allSpicCentVel>0);
                
                figure; tiledlayout(1,2)
                nexttile
                fitmethis(data,'pdist',numMatches,'output','off');
                title(sprintf('Distribution fitting output (top 4 matches)\n')) 
                
                nexttile
                distribution = 'gp';
                pd = fitdist(data,distribution);
                
                histfit(data,11,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.1f');
                grid on
                %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
                title(sprintf('Spicule centroid velocity histogram\nMean spicule centroid velocity\n= %d m/s\nDistribtion: %s\nk = %d\nsigma = %d\ntheta = %d', avSpicCentVel, distribution, pd.k, pd.sigma, pd.theta))
                xlabel('Velocity in m/s')
                xtickangle(90)
                ylabel('Probability')
                
                %plot centroid velocity osc freq
                maximum = max(spicCentVelFreq);
                if maximum>0
                    minimum = min(spicCentVelFreq(spicCentVelFreq>0));
                    figure
                    histfit(spicCentVelFreq(spicCentVelFreq>0),binNumber,'normal');
                    title(sprintf('Spicule centroid velocity oscillation frequency histogram\nMean spicule centroid velocity oscillation frequency = %d Hz', avCentVelFreq))
                    xlabel('Centroid Oscillation Frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule centroid velocity')
                end
                
                
                
                
                %plot top velocity
                minimum = min(allSpicTopVel(allSpicTopVel>0));
                maximum = max(allSpicTopVel);
                
                data = allSpicTopVel(allSpicTopVel>0);
                
                figure; tiledlayout(1,2)
                nexttile
                fitmethis(data,'pdist',numMatches,'output','off');
                title(sprintf('Distribution fitting output (top 4 matches)\n')) 
                
                nexttile
                distribution = 'gp';
                pd = fitdist(data,distribution);
                
                histfit(data,11,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.1f');
                grid on
                %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
                title(sprintf('Spicule top velocity histogram\nMean spicule top velocity\n= %d m/s\nDistribtion: %s\nk = %d\nsigma = %d\ntheta = %d', avSpicTopVel, distribution, pd.k, pd.sigma, pd.theta))
                xlabel('Velocity in m/s')
                xtickangle(90)
                ylabel('Probability')
                
                %plot top velocity osc freq
                maximum = max(spicTopVelFreq);
                if maximum>0
                    minimum = min(spicTopVelFreq(spicTopVelFreq>0));
                    figure
                    histfit(spicTopVelFreq(spicTopVelFreq>0),binNumber,'normal');
                    title(sprintf('Spicule top velocity oscillation frequency histogram\nMean spicule top velocity oscillation frequency = %d Hz', avTopVelFreq))
                    xlabel('Top Oscillation Frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule top velocity')
                end
                
                
                
                
                %plot bottom velocity
                minimum = min(allSpicBottVel(allSpicBottVel>0));
                maximum = max(allSpicBottVel);
                
                data = allSpicBottVel(allSpicBottVel>0);
                
                figure; tiledlayout(1,2)
                nexttile
                fitmethis(data,'pdist',numMatches,'output','off');
                title(sprintf('Distribution fitting output (top 4 matches)\n')) 
                
                nexttile
                distribution = 'gp';
                pd = fitdist(data,distribution);
                
                histfit(data,11,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.1f');
                grid on
                %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
                title(sprintf('Spicule bottom velocity histogram\nMean spicule bottom velocity\n= %d m/s\nDistribtion: %s\nk = %d\nsigma = %d\ntheta = %d', avSpicBottVel, distribution, pd.k, pd.sigma, pd.theta))
                xlabel(sprintf('Velocity in m/s'))
                xtickangle(90)
                ylabel('Probability')
                
                plot bottom velocity osc freq
                maximum = max(spicBottVelFreq);
                if maximum>0
                    minimum = min(spicBottVelFreq(spicBottVelFreq>0));
                    figure
                    histfit(spicBottVelFreq(spicBottVelFreq>0),binNumber,'normal');
                    title(sprintf('Spicule bottom velocity oscillation frequency histogram\nMean spicule bottom velocity oscillation frequency = %d Hz', avBottVelFreq))
                    xlabel('Bottom Oscillation Frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule bottom velocity')
                end
                
                
                
                
                %plot centroid osc freq
                maximum = max(spicCentFreq);
                if maximum>0
                    minimum = min(spicCentFreq(spicCentFreq>0));
                    
                    data = spicCentFreq(spicCentFreq>0);
                    
                    figure; 
                    
%                     tiledlayout(1,2)
%                     nexttile
%                     fitmethis(data,'pdist',numMatches,'output','off');
%                     title(sprintf('Distribution fitting output (top 4 matches)\n')) 
%                     
%                     nexttile
                    distribution = 'ev';
                    pd = fitdist(data,distribution);
                    
                    histfit(data,4,distribution);
                    yt = get(gca, 'YTick');
                    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                    xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.3f');
                    grid on
                    %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                    
                    title(sprintf('Spicule centroid oscillation histogram\nMean centroid oscillation frequency = %d Hz\nMean centroid oscillation period = %d s\nDistribtion: %s\nmu = %d\nsigma = %d', avOscFreq, 1/avOscFreq, distribution, pd.mu, pd.sigma))
                    xlabel('Oscillation frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule centroid')
                end
                
                
                
                %plot top osc freq
                maximum = max(spicTopFreq);
                if maximum>0
                    minimum = min(spicTopFreq(spicTopFreq>0));
                    
                    data = spicTopFreq(spicTopFreq>0);
                    
                    figure; 
%                     tiledlayout(1,2)
%                     nexttile
%                     fitmethis(data,'pdist',numMatches,'output','off');
%                     title(sprintf('Distribution fitting output (top 4 matches)\n')) 
%                     
%                     nexttile
                    distribution = 'ev';
                    pd = fitdist(data,distribution);
                    
                    histfit(data,4,distribution);
                    yt = get(gca, 'YTick');
                    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                    xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.3f');
                    grid on
                    %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                    
                    title(sprintf('Spicule top oscillation histogram\nMean top oscillation frequency = %d Hz\nMean top oscillation period = %d s\nDistribtion: %s\nmu = %d\nsigma = %d', avTopOscFreq, 1/avTopOscFreq, distribution, pd.mu, pd.sigma))
                    xlabel('Top oscillation frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule top position')
                end
                
                
                
                %plot bottom osc freq
                maximum = max(spicBottFreq);
                if maximum>0
                    minimum = min(spicBottFreq(spicBottFreq>0));
                    
                    data = spicBottFreq(spicBottFreq>0);
                    
                    figure; 
                    
%                     tiledlayout(1,2)
%                     nexttile
%                     fitmethis(data,'pdist',numMatches,'output','off');
%                     title(sprintf('Distribution fitting output (top 4 matches)\n')) 
%                     
%                     nexttile
                    distribution = 'ev';
                    pd = fitdist(data,distribution);
                    
                    histfit(data,4,distribution);
                    yt = get(gca, 'YTick');
                    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                    xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.3f');
                    grid on
                    %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                    
                    title(sprintf('Spicule bottom oscillation histogram\nMean bottom oscillation frequency = %d Hz\nMean top oscillation period = %d s\nDistribtion: %s\nmu = %d\nsigma = %d', avBottOscFreq, 1/avBottOscFreq, distribution, pd.mu, pd.sigma))
                    xlabel('Bottom oscillation frequency in Hz')
                    ylabel('Probability')
                else
                    fprintf('\nNo oscillation found in any spicule bottom position')
                end
                
                
                
                %plot lifetime
                minimum = min(spicLifetimeData(spicLifetimeData>0));
                maximum = max(spicLifetimeData);
                
                data = spicLifetimeData(spicLifetimeData>0);
                
                figure; 
                
%                 tiledlayout(1,2)
%                 nexttile
%                 fitmethis(data,'pdist',numMatches,'output','off');
%                 title(sprintf('Distribution fitting output (top 4 matches)\n')) 
%                 
%                 nexttile
                distribution = 'loglogistic';
                pd = fitdist(data,distribution);
                
                histfit(data,47,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, 250]); xtickformat('%d');
                grid on
                %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
                title(sprintf('Spicule lifetime histogram\nMean spicule lifetime = %d s\nDistribtion: %s\nmu = %d\nsigma = %d', avLifetime, distribution, pd.mu, pd.sigma))
                xlabel('Lifetime in s')
                ylabel('Probability')
                
                
                
                
                %plot spic density
                minimum = min(spicDensity(spicDensity>0));
                maximum = max(spicDensity);
                
                data = spicDensity(spicDensity>0);
                
                figure; 
                
%                 tiledlayout(1,2)
%                 nexttile
%                 fitmethis(data,'pdist',numMatches,'output','off');
%                 title(sprintf('Distribution fitting output (top 4 matches)\n')) 
%                 
%                 nexttile
   
                distribution = 'gev';
                pd = fitdist(data,distribution);
                
                histfit(data,5,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.2f');
                grid on
                %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
                title(sprintf('Spicule density histogram\nMean spicule density = %d spic/m^2\nActive spicule oscillation frequency = %d Hz\nDistribtion: %s\nk = %d\nmu = %d\nsigma = %d', avDensity, spicDensityFreq, distribution, pd.k, pd.mu, pd.sigma))
                xlabel('Density in spic/m^2')
                ylabel('Probability')
                
                
                
                
                %plot birth rate density
                minimum = min(spicBirthDensity(spicBirthDensity>0));
                maximum = max(spicBirthDensity);
                
                data = spicBirthDensity(spicBirthDensity>0);
                
                figure; 
%                 tiledlayout(1,2)
%                 nexttile
%                 fitmethis(data,'pdist',numMatches,'output','off');
%                 title(sprintf('Distribution fitting output (top 4 matches)\n')) 
%                 
%                 nexttile
                distribution = 'gev';
                pd = fitdist(data,distribution);
                
                histfit(data,7,distribution);
                yt = get(gca, 'YTick');
                set(gca, 'YTick', yt, 'YTickLabel', yt/numel(data));
                xlim([0, maximum+(maximum-minimum)*0.1]); xtickformat('%.1f');
                grid on
                %                 title(sprintf('Spicule length histogram\nMean spicule length = %d meters\nDistribtion: %s\nMu = %d\nSigma = %d', avLength, distribution, pd.mu, pd.sigma))
                
                title(sprintf('Spicule birth rate density histogram\nMean birth rate density = %d new spic/s/m^2\nSpicule birth rate oscillation frequency = %d Hz\nDistribtion: %s\nk = %d\nmu = %d\nsigma = %d', avNewSpicDensity, spicBirthDensityFreq, distribution, pd.k, pd.mu, pd.sigma))
                xlabel('Birth rate density in new spic/s/m^2')
                ylabel('Probability')
                
                
                
                plotHistograms = 0;
                
                
                fprintf('\n-----------\nTotal number of unique spicules analyzed = %d',totNumSpic);
                
                
                
            end
            
            
            plotStats = 0;
            
            fprintf('\n-----------\nEnter 1 to plot another spicule, enter 0 to end script:');
            plotStats = input(' ');  %get user input for spicule index to be presented
            
        end
        
        
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end  %close while ii<=numFrames


%% play original frames

clear;
close all;


vid = VideoReader(videoName); %reads video and gives information
vid.Duration;
vid.FrameRate;
totalFrames = vid.FrameRate()*vid.Duration;
readVideo;



%INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vid.CurrentTime = 0;
playSpeed = vid.FrameRate; %vid.FrameRate for normal speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%maintains axis size
figure;
ax = axes;


while vid.hasFrame()
    temp = vid.readFrame();
    imshow(temp,'Parent', ax); %shows current temp frame and allows window adjustment
    pause(1.0/playSpeed); %play video at its framerate
end

%% extract single frame at given time

clear;
close all;

vid = VideoReader(videoName); %reads video and gives information
vid.Duration;
vid.FrameRate;
totalFrames = vid.FrameRate()*vid.Duration;
readVideo;

%INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vid.CurrentTime = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%show single frame at given time

vidFrame = vid.readFrame();

%shows read video frame
figure;
imshow(vidFrame);

%% play uint8 frames

%INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loopLimit = 3;
playSpeed = 3; %vid.FrameRate for normal speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii=1;
loop = 0;
array2play = temp3;

figure;
ax = axes;

while ii<=numFrames && loop <= loopLimit %loop video
    
    imshow(squeeze(array2play(ii,:,:,:)),'Parent', ax);
    pause(1.0/playSpeed);
    ii=ii+1;
    
    
    %return to first frame after reaching final frame
    
    if(ii==numFrames)
        ii=1;
        loop = loop+1;
    end
    
end

%% play double type frames

%INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loopLimit = 3;
playSpeed = vid.FrameRate; %vid.FrameRate for normal speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ii=1;
loop = 0;

figure;
ax = axes;

while ii<=numFrames && loop <= loopLimit %loop video
    
    colormap('gray');
    imagesc(squeeze(temp3Doub(ii,:,:,:)),'Parent', ax);
    
    pause(1.0/playSpeed);
    ii=ii+1;
    
    
    %return to first frame after reaching final frame
    
    if(ii==numFrames)
        ii=1;
        loop = loop+1;
    end
    
end

%% show nth uint8 frame

%INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displayFrame = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

array2play = temp3;

figure;
imshow(squeeze(array2play(displayFrame,:,:,:)));

%% show nth double type frame

%INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displayFrame = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
imagesc(squeeze(temp3Doub(displayFrame,:,:,:)));
colormap('gray');
