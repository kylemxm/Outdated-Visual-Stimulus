%% Analysis of max std data from Drifting_Grating_analysis_STD function
%This program is intended to analyze data from the driting frating program,
%specifically using the maximum peak value for each orientation and the #
%of standard deviations from the mean each peak maximum is. 

%updated 8/22/16 to work with area under curve rather than peak

%inputs:
% Area_Under_Curve_Data, STD_from_Mean

%% Functions to add: ROI deletion, setting threshold, Half width at half max, better OSI calculation 

%% Select Contra data
[Filename, Pathname]=uigetfile('*.mat', 'Select your Contra data');
load([Pathname Filename]);
Max_Contra=Area_Under_Curve_Data;
STD_Contra=STD_from_Mean;
 %set current directory to pathname
    cd(Pathname);
%% Select Ipsi data
[Filename, Pathname]=uigetfile('*.mat', 'Select your Ipsi data');
load([Pathname Filename]);
Max_Ipsi=Area_Under_Curve_Data;
STD_Ipsi=STD_from_Mean;
clear Area_Under_Curve_Data STD_from_Mean Filename

%record total number of cells
Cell_total=length(STD_Contra);

%set all negative #s to 0 
Max_Contra(Max_Contra<0)=0;
Max_Ipsi(Max_Ipsi<0)=0;

%% Add feature to choose ROIs to exclude


%% FInd ROIs that are visually responsive, make smaller matrices
%set non-sig STD orientations to 0 threshold CURRENTLY 1 (have to change
%lines 40, 41, 78 and 79 to change threshold). 
STD_Contra(STD_Contra<.95)=0;
STD_Ipsi(STD_Ipsi<.95)=0;
% Collapse to a vector
SumC=sum(STD_Contra,2);
SumI=sum(STD_Ipsi,2);
%Add contra and ipsi vectors
SumCI=SumC+SumI;
SumCI(SumCI>0)=1;
%Make new matrices of Max and std data from only responsive ROIs
Rindices=find(SumCI==1);
MaxC=Max_Contra(Rindices,:);
MaxI=Max_Ipsi(Rindices,:);
STDC=STD_Contra(Rindices,:);
STDI=STD_Ipsi(Rindices,:);
clear SumC SumI SumCI

%% Separate ROIs into contra, ipsi or binoc
%find all the sig ROIs again
% Collapse to a vector
SumC=sum(STDC,2);
SumI=sum(STDI,2);
%set all >1 to 1
SumC(SumC>0)=1;
SumI(SumI>0)=1;
%Add contra and ipsi vectors
SumCI=SumC+SumI;
%set values less than two (monocular ROIs) to 0
SumCI(SumCI<2)=0;
%select binoc ROIs
SumCI(SumCI==2)=1;
%Find Binocular ROIs
Binoc=find(SumCI==1);
%Find Contra ROIs
Contra=find(SumCI==0 & SumC==1);
%Find Ipsi ROIs
Ipsi=find(SumCI==0 & SumI==1);

%% Find Peak for each ROI 
%set non-sig STD orientations to 0
STDC(STDC<.95)=0;
STDI(STDI<.95)=0;
%set sig STD orientations to 1
STDC(STDC>0)=1;
STDI(STDI>0)=1;

%max matrices with all non sig orientations set to 0
Cmax=MaxC.*STDC;
Imax=MaxI.*STDI;
%Find max indice of each eye WARNING FOR ROWS WITH ALL ZEROS THIS SETS THE
%MAX INDEX TO 1. IT SHOULD NOT BE A PROBLEM, AS I WILL NEVER REFERENCE IT
%BECAUSE I CALL INDICES FROM THESE MATRICES USING BINOC CONTRA AND IPSI,
%BUT IT IS SOMETHING TO BE AWARE OF. I COULD MAKE ALL 0'S NAN. 
[~,C_max]=max(Cmax,[],2);
[~,I_max]=max(Imax,[],2);
%find max indices for binoc, biggest index either ipsi or contra. Means ODI
%is weird for values with an offset as it settles on biggest orientation,
%not using biggest orientation for both ipsi and contra.
B_max=zeros(numel(Binoc),1);
for ii=1:numel(Binoc)
    [~,K]=max(Cmax(Binoc(ii),:),[],2);
    [~,L]=max(Imax(Binoc(ii),:),[],2);
    M=max(K,L);
    B_max(ii)=M;
end

%% Calculate and spit out ODI,OSI,DSI,and matching data

%% binoc ROIs 
%calculate ODI for each ROI, as ipsi/(Contra+ipsi). This could also be done
%with simple indexing, but a for statement should work as well and still
%relatively fast. small matrices. 
ODI_B=zeros(numel(Binoc),1);
for jj=1:numel(Binoc)
    ODI_B(jj)=MaxI(Binoc(jj),B_max(jj))/(MaxI(Binoc(jj),B_max(jj))+MaxC(Binoc(jj),B_max(jj)));
end

%% Find binocular matching offset
%find contra and ipsi preferred orientation, column 1 is contra, column 2
%is ipsi
Match_B=zeros(numel(Binoc),2);
for kk=1:numel(Binoc)
    Match_B(kk,1)=C_max(Binoc(kk))*30-30;
    Match_B(kk,2)=I_max(Binoc(kk))*30-30;
end
%Set on 0 to 180 scale (if # higher than 180, subtract 180).
for jj=1:numel(Binoc)*2;
    if Match_B(jj)>=180;
        Match_B(jj)=Match_B(jj)-180;
    end
end

%find absolute value of difference
Match_B=abs(Match_B(:,1)-Match_B(:,2));
%if difference is greater than 90 degrees, subtract from 180.
clear kk jj
for jj=1:numel(Match_B)
    if Match_B(jj)>90
        Match_B(jj)=180-Match_B(jj);
    end
end
        

%% Contra ROIs
%calculate ODI for each ROI, as ipsi/(Contra+ipsi).
ODI_C=zeros(numel(Contra),1);
for jj=1:numel(Contra)
    ODI_C(jj)=MaxI(Contra(jj),C_max(Contra(jj)))/(MaxI(Contra(jj),C_max(Contra(jj)))+MaxC(Contra(jj),C_max(Contra(jj))));
end
%calculate OSI
OSI_C=zeros(numel(Contra),1);
%shift contra data so largest value 1st
[rows,columns]=size(MaxC);
aligned_Contra=zeros(rows,columns);
for ii=1:length(C_max)
    aligned_Contra(ii,:)=circshift(MaxC(ii,:),[0,(columns-C_max(ii)+1)]);
end
%Actual OSI calculation (max-orthogonal)/(max+orthogonal)
for jj=1:numel(Contra)
    OSI_C(jj)=(aligned_Contra(Contra(jj),1)-mean([aligned_Contra(Contra(jj),4),aligned_Contra(Contra(jj),10)]))/(aligned_Contra(Contra(jj),1)+mean([aligned_Contra(Contra(jj),4),aligned_Contra(Contra(jj),10)]));
end
%Direction selectivity (max-opposite)/(max+opposite)
Dir_C=zeros(numel(Contra),1);
for jj=1:numel(Contra)
    Dir_C(jj)=(aligned_Contra(Contra(jj),1)-aligned_Contra(Contra(jj),7))/(aligned_Contra(Contra(jj),1)+aligned_Contra(Contra(jj),7));
end

%% Ipsi ROIs
%calculate ODI for each ROI, as ipsi/(Contra+ipsi).
ODI_I=zeros(numel(Ipsi),1);
for jj=1:numel(Ipsi)
    ODI_I(jj)=MaxI(Ipsi(jj),I_max(Ipsi(jj)))/(MaxI(Ipsi(jj),I_max(Ipsi(jj)))+MaxC(Ipsi(jj),I_max(Ipsi(jj))));
end
%calculate OSI
OSI_I=zeros(numel(Ipsi),1);
%shift ipsi data so largest value 1st
[rows,columns]=size(MaxI);
aligned_Ipsi=zeros(rows,columns);
for ii=1:length(C_max)
    aligned_Ipsi(ii,:)=circshift(MaxI(ii,:),[0,(columns-I_max(ii)+1)]);
end
%Actual OSI calculation (max-orthogonal)/(max+orthogonal)
for jj=1:numel(Ipsi)
    OSI_I(jj)=(aligned_Ipsi(Ipsi(jj),1)-mean([aligned_Ipsi(Ipsi(jj),4),aligned_Ipsi(Ipsi(jj),10)]))/(aligned_Ipsi(Ipsi(jj),1)+mean([aligned_Ipsi(Ipsi(jj),4),aligned_Ipsi(Ipsi(jj),10)]));
end
%Direction selectivity (max-opposite)/(max+opposite)
Dir_I=zeros(numel(Ipsi),1);
for jj=1:numel(Ipsi)
    Dir_I(jj)=(aligned_Ipsi(Ipsi(jj),1)-aligned_Ipsi(Ipsi(jj),7))/(aligned_Ipsi(Ipsi(jj),1)+aligned_Ipsi(Ipsi(jj),7));
end

%% Binoc OSI and direction (for each eye).
%OSI, column 1 is contra column 2 is Ipsi
OSI_B=zeros(numel(Binoc),2);
for jj=1:numel(Binoc)
    OSI_B(jj,1)=(aligned_Contra(Binoc(jj),1)-mean([aligned_Contra(Binoc(jj),4),aligned_Contra(Binoc(jj),10)]))/(aligned_Contra(Binoc(jj),1)+mean([aligned_Contra(Binoc(jj),4),aligned_Contra(Binoc(jj),10)]));
    OSI_B(jj,2)=(aligned_Ipsi(Binoc(jj),1)-mean([aligned_Ipsi(Binoc(jj),4),aligned_Ipsi(Binoc(jj),10)]))/(aligned_Ipsi(Binoc(jj),1)+mean([aligned_Ipsi(Binoc(jj),4),aligned_Ipsi(Binoc(jj),10)]));
end
% direction selectivity. Column 1 is Contra Column 2 is Ipsi
Dir_B=zeros(numel(Binoc),2);
for jj=1:numel(Binoc)
    Dir_B(jj,1)=(aligned_Contra(Binoc(jj),1)-aligned_Contra(Binoc(jj),7))/(aligned_Contra(Binoc(jj),1)+aligned_Contra(Binoc(jj),7));
    Dir_B(jj,2)=(aligned_Ipsi(Binoc(jj),1)-aligned_Ipsi(Binoc(jj),7))/(aligned_Ipsi(Binoc(jj),1)+aligned_Ipsi(Binoc(jj),7));
end

%% Calculate average deltaF/F for each category (binoc/Contra/Ipsi

%binoc
Binoc_F=zeros(numel(Binoc),2);
for jj=1:numel(Binoc)
   Binoc_F(jj,1)=MaxC(Binoc(jj),C_max(Binoc(jj)));
   Binoc_F(jj,2)=MaxI(Binoc(jj),I_max(Binoc(jj)));
end
Binoc_F=mean(Binoc_F,1);

%Contra
Contra_F=zeros(numel(Contra),1);
for jj=1:numel(Contra)
   Contra_F(jj)=MaxC(Contra(jj),C_max(Contra(jj)));
end
Contra_F=mean(Contra_F);

%Ipsi
Ipsi_F=zeros(numel(Ipsi),1);
for jj=1:numel(Ipsi)
   Ipsi_F(jj)=MaxI(Ipsi(jj),I_max(Ipsi(jj)));
end
Ipsi_F=mean(Ipsi_F);

clear I_max ii Imax jj k kk L M Max_Contra Max_Ipsi MaxC MaxI rows STD_Contra STD_Ipsi STDC STDI SumC SumCI SumI aligned_Contra aligned_Ipsi C_max columns 
clear B_max Cmax K Pathname Rindices