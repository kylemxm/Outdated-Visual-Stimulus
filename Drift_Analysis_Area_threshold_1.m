% program is designed to import ROI intesity data from multimeasure
% imageJ plugin and connect it with the timestamps from the XML file and
% the CSV file to get the data on the stimuli out. event identity stored as
% a csv file stimevents.m

%Modified 7-19-16 by KJ to take the value for each ROI as the area under
%the curve from 0-5 seconds of the average trace, (stimulus duration) and 
%calculate how many STDs the mean of average trace from 0-5 seconds is from
%F mean (0).
 
clear all
%% open and import .csv file
[filename, pathname] = uigetfile('*.csv', 'Select your .csv ROI data file');
%cancel if user clicks cancel
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
else
    %set current directory to pathname
    cd(pathname);
    %set file to path string
   file = [pathname filename]; 
   %exports numeric data from sheet 1 of the excel sheet to a matrix
   ROIdata_raw = csvread(file);
end

%% open and import XML file
  [filename, pathname] = uigetfile('*.xml', 'Select your .xml file from Prairie data');
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
else
   file = [pathname filename]; 
   data = xml2struct(file);
end
%determine length
framelength = length(data.PVScan.Sequence.Frame);
%create cell of proper size
timestamps=zeros(framelength,1);
%loop to extract timestamps, either relative or absolute
if (data.PVScan.Sequence.Frame{1, 2}.Attributes.relativeTime>0)
for i = 1:framelength
    timestamps(i,1)=str2double(data.PVScan.Sequence.Frame{1, i}.Attributes.relativeTime);
    end
else
  for i = 1:framelength
    timestamps(i,1)=str2double(data.PVScan.Sequence.Frame{1, i}.Attributes.absoluteTime);
  end  
end
clear i


%% Select Electrophys data
%We now have our ROI data (data) and our timestamps (timestamps). Now we
%need the event timestamps.
%open file
[filename, pathname] = uigetfile('*.csv', 'Select your .CSV file for electrophysiology');
if isequal(filename,0) || isequal(pathname,0)
    disp('action canceled')
else
   file = [pathname filename]; 
   electrophystime = csvread(file,2);
end
   
%by shifting files up and down, I can find the start and end of square
%waves
roundedtimes=round(electrophystime);
shiftdown=[0,0;roundedtimes];
shiftup=[roundedtimes;0,0];
%find differences
firstevent=(eq(shiftdown(:,2),shiftup(:,2)));
%1 now means there was an event
secondevent=1.-firstevent;
%event will now be same size as original roundedtimes matrix
thirdevent=secondevent(1:end-1);
%multiply together with rounded times column 1 to get timestamps where
%events happened. because it is relative time, this corresponds to
%handles.timestamps for image data (once divided by 1000 to be seconds).
event_timestamps=(roundedtimes(:,1).*thirdevent)/1000;
%remove zeros from event_timestamps 
Events_notzero=event_timestamps(event_timestamps~=0);
%take every twentieth value of a matrix, this is the start of each block of
%stimuli
startevent=Events_notzero(1:2:end);
%take every twentieth value from 20, this is the blank frames
blankframes=Events_notzero(2:2:end);

%add timestamps to replace first column of data
ROIdata_raw(:,1)=timestamps;

%% find F for each ROI and convert data to delta F/F
%calculate number of ROIs
[~,numROIs]=size(ROIdata_raw);
%account for column of timestamps
numROIs=numROIs-1;
%preallocate matrix for blank frame data (no visual stim period b/t stims)
F_data=zeros(1,(numROIs+1));
%loop to collect data (assumes FIVE seconds of inter stimulus interval)
for ii=1:length(blankframes-1)
    Chunk = (ROIdata_raw(ROIdata_raw(:,1) >=(blankframes(ii))&ROIdata_raw(:,1) <=(blankframes(ii)+5),:));
    F_data=[F_data;Chunk];
end
%calculate mean of each ROI (first column is time, so cut it off)
meanF=mean(F_data(2:end,:),1);
meanF=meanF(2:end);

%Find standard deviation for each ROI
standarddev=std(F_data(2:end,:),1);
standarddev=standarddev(2:end);
standarddev=standarddev./meanF.*100;  %still working on this. should it be -1 then *100?

%normalize every ROI by its F
ROIdata=ROIdata_raw;
clear ii jj
for ii=1:numROIs
ROIdata(:,ii+1)=((ROIdata_raw(:,ii+1)./meanF(ii))-1).*100;
end

% %% DeltaF/F ROI data
% %use first 20 frames. Alternatively I could use the mean of all the inter
% %stimulus periods, but that is a hell of a lot of coding. Need to do it
% %though. 
% [~,numROIs]=size(ROIdata_raw);
% for i=1:numROIs-1
%     ROIdata(:,i+1)=((ROIdata_raw(:,i+1)./mean(ROIdata_raw(1:20,i+1)))-1).*100;
% end
% 

%% process by event
%.mat file that has the stimulus order saved. 
load('stimevents');

%% updating to deal with stimevents as a matrix, assuming every row is a new repeat 6-29-16
%find # of events
stimnumber=numel(stimevents);
%rotate 90 degrees counterclockwise 3 times
stimevents=rot90(stimevents,3);
%flip it left to right. now basically the stimevents, rather than reading
%left to right, top to bottom like english reads top to bottom, left to
%right, which is the matlab indexing.
stimevents=fliplr(stimevents);
%reshape into vector.
stimevents=reshape(stimevents,1,stimnumber);

%organize into matrices
eventmatrix(:,1)=stimevents;
eventmatrix(:,2)=startevent(1:length(stimevents));

%number of unique stims
uniquestims=unique(eventmatrix(:,1));
%presentations per stim
counts = histc(eventmatrix(:,1), uniquestims);

% extract time series data. uneven lengths, so store data as a cell. I want 5
% seconds before to 5 seconds after each stim ASSUMES A FIVE SECOND STIM
for i=1:length(uniquestims)
    %timestamps for uniquestim i
    stim=uniquestims(i);
    uniquetimestamps=eventmatrix(eventmatrix(:,1)==stim,2);
    %cell array to store data
    Datacell(i,1)={stim};
    %loop to get all events for each unique stim
    for j=1:counts(i)
       %data range for each stim
        Datacell(i,j+1)={ROIdata(ROIdata(:,1) >=(uniquetimestamps(j)-5)&ROIdata(:,1) <=(uniquetimestamps(j)+10),:)};
        %normalize times of each event to be on a scale with 0 as the event
        %start
        Datacell{i,j+1}(:,1)=Datacell{i,j+1}(:,1)-uniquetimestamps(j);
    end
end
 clear i j l 
 
 
 
 %normalize times of each event to be on a -5 to 15 scale

 %% Seperate by ROI  
%plot by ROI with an interruption so I can exit the loop
%figure with panels for the unique stims

%blank matrix to save max fluorescence data and std distance
Area_Under_Curve_Data=zeros(numROIs,length(uniquestims));
STD_from_Mean=zeros(numROIs,length(uniquestims));

%loop through the ROIs. 
for l=1:numROIs
%create new figure with ROI name
figure('name',sprintf('Plot of ROI %d',l),'numbertitle','off','position',[250,500,1000,700])
    %loop through stims
for i=1:length(uniquestims)
    %loop through the number of stims
    for j=1:counts(i)
        %calculate size of graph (requires even number of unique stims.
 graphsize=length(uniquestims)/2;     
subplot(2,graphsize,i);
plot(Datacell{i,j+1}(:,1),Datacell{i,j+1}(:,l+1))

hold on
    end
%plot an average curve, 3rd party function. 
[avgH,avgData]=plotAverage(gcf);
averageline=avgData{1,1}(:,1:2);
%find data b/t 0 and FIVE seconds
databetween=averageline(averageline(:,1)>=0&averageline(:,1)<=5,:);
%take area under the curve
AUC=trapz(databetween(:,1),databetween(:,2));
%save area under the curve
Area_Under_Curve_Data(l,i)=AUC;

%% find mean for databetween
mean_databetween=mean(databetween(:,2));
%STDs away from mean F
STD_Distance=(mean_databetween/standarddev(l));
%save
STD_from_Mean(l,i)=STD_Distance;


%% Make it obvious which orientations are 2 std or more above mean 
legend=legend(sprintf('AUC = %0.3f',AUC),sprintf('STD = %0.3f',STD_Distance)); 
if STD_Distance>=0.95;
    legend.FontWeight='bold';
    set(gca,'Color',[0.4 0.8 0.8]);
end
clear averageline avgdata AUC STD_Distance legend mean_databetween
end
    w=waitforbuttonpress;
    if w==0
        close all
        continue
    else
      break
    end
end
   
   
    clear avgData avgH blankframes Chunk counts data databetween Datacell electrophystime event_timestamps eventmatrix Events_notzero F_data file 
    clear filename firstevent framelength graphsize i ii j l meanF numROIs pathname ROIdata_raw roundedtimes secondevent shiftdown 
 clear shiftup standarddev startevent stim stimevents stimnumber thirdevent timestamps uniquestims uniquetimestamps w ROIdata
 
 save('AUC and STD data');