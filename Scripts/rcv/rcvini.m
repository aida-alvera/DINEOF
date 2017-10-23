%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter definition file %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add simplification by streaming info from preprecol through files to
% read, then moove these param towards rcv

clear all
close all

%LIST OF INPUT FILES possibly required by RCV
%Recolour Visualisation

% Outputs of PREPRECOL :
filenameinvectdates=['vectseldates.gher']; 
filenameinNodata1=['vistpabsdatall.gher'];
filenameinNodata2=['vistpabsdatspsel.gher'];
filenameinMask=['belcolour_subregion_mask.gher'];
filenameinMasktempred=['belcolour_region_tempreducedmask.gher'];
filenameinBackgd=['2Dbelcolour_region_period_background.gher'];
filenameinAnomaly=['2Dbelcolour_region_period_anomaly.gher'];
filenameinOridat=['2Dbelcolour_region_period_dat.gher'];

% Outputs of DINEOF analysis:
filenameinmeandata=['meandata.val'];
filenameinneofretained=['neofretained.val'];
filenameinFilleddat=['F2Dbelcolour_region_period_datfilled.gher']; 

% positions of stations to be extracted:
filenameinstations=['coordstations.txt'];

% DEFINING THE DATASET PROCESSED BY DINEOF  % stream from preprecol and
% list.DAT
datareal2D=0     % 
%%%%%SENSOR%%%%%%%
sensor='M'; %for MERIS
%sensor='A'; %for MODIS AQUA
%%%%%PARAMETER%%%%%%%
%param='CHL';
param='TSM';
%param='SST';
%%%%%UNITS%%%%%%%
%unit='(°C)';
unit='mg/l';
%unit='µg/l';
%%%%%PERIOD%%%%%%%
period='all'; % produce numb images from preprecol.ini
%%%%%REGION%%%%%%%
%stream lat long ini from preprecol.ini according to cropped zone!!!!
%for all belcolour focus
 %latini=60.002; 
 %longini=-4.000;
%for recolour focus
 latini=52.500; 
 longini=-4.000;
% grid dx dy
dlat=0.009;
dlong=0.014;

% DEFINING THE PREPROCESSING APPLIED BY PREPRECOL  % stream from preprecol
workedinlog=1   %!!!!!!!!comments to adapt   =1 to rebuild and record data in real units if it was processed in log10; =0 to output in log10 if processed in log10% =0 if processed in raw units and need to be outputed in raw units.
workedinanomaly=1    % option for rebuilding filled and original fields as anomaly only (workedinanomaly==0) or as real values ,with bckgrd addedd (workedinanomaly==1)
AL=1  %worked on the anomaly of the log10 of the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST OF OPTIONS ACCORDING TO USER REQUESTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VISUALISATION OF STATISTICS OF ABSENT/ELIMINATED/BACKGROUND DATA 
absentdata=0
vistempelim=0
background=0  % rebuild, then display or save background
extractserie=1 % =1 if extract time series at particular stations
readplotexts=0
tempaverage=1 % =1 if calculate multitemporal averages
perioday=7  % parameter in days for the multitemporal averages and station extraction, adapt to a variable number of days/month if monthly mean requested


% COMPARING ORIGINAL VERSUS FILLED DATA
oridata=1    % not used for now: adapt to read 2D_data and 2D_anomaly slice/sclice, and display according to request
filleddata=0 % not in used for now: adapt to read F2D slice/slice, rebuilding done for now internally through eofs (next param)
displayby2=0 % if display 2/2 after producing the bin
nbtimeperfig=5;
norimin=1
norimax=0   % if ==0 ...all

errorfields=0


% REBUILDING THE EOFS %
reverteof1=1  % =1 to revert the sign of the first eofs (both in space and time)
revvect=[1 2] % vector mentionning the eofs to revert if previous parameter is set to 1
%%%%%%%%%
spatialmodes=0% =1 to rebuild + save or display the spatial eofs
showneofsp=0  % number max of spatial modes do display, if ==0, then takes the number of eof retained
eofmvlsng=1
%%%%%%%%%
tempmodes=0  % =1 to save or display the temporal eofs
datetpeofs=1 % option to add : =1 if display only at dates of selected images; =2 if display only at interpolated dates; =3 if display both

if sensor=='M'
    tpyi=2003
    tpyf=2007.1
    climtpeofmin=-0.2
    climtpeofmax=0.2
elseif sensor=='A'
    tpyi=2002.4
    tpyf=2007.1
    climtpeofmin=-0.1
    climtpeofmax=0.1
end
showneof2=3 % nb of tp eofs to display replace in code showneof by useneof for reconstructions, and (showneof by showneof2)  everywhere to 'useneof'%%%

%REBUILDING THE FILLED DATA %(add extraction at selected dates)
%for diagnosis : checking eofs contributions at selected dates
% recomposing images at selected dates (internally from the eofs) 
rebuildnoi=0  %without interpolation
checkcontrib=0 % displaying the contribution of each eof to the image reconstruction
nnoimin=170 %image 170
nnoimax=170

% to recompose images at regular time steps (from temporal eofs interpolation) 
% and to produced multitemporal averages 
rebuild=0  % with interpolation
ntintmin=0  %  % (if=0 then it takes automatically the date of the first non nan interpolated eof)
ntintmax=0 %351  % (if=0 then it takes automatically the date of the last non nan interpolated eof, if different then 0,it sets the number of interpolated images to produce)
djint=1 %delta days between rebuilt images from interpolated temporal eofs

rebuildmissim=0 % if requested to rebuild from interpolated eofs the images eliminated for quality check

%DEFINING THE STATIONS FOR WHICH TIME SERIES SHOULD BE EXTRACTED !! adapt
%for multistations loaded in a file  : station name, lat, long     % now in
%rcv
%station  names and order:         
%1: Scheelde turbidity maximum
%2: station 230                  
%3: station 330                   
%4: NorthDogger           %%%%%not used in recolour focus :lat 55.683, long
%2.280
%5: OysterGrounds         %%%%%not used in recolour focus: lat 54.414, long
%4.042
%6: Warp                  
%7: WestG                 
% lats     longs   in decimal degree: 
%latcheck=[51.47,51.308,51.433,51.526,51.990]; %north
%longcheck=[3.23,2.849,2.809,1.030,2.090]; %east
%nstmax=size(latcheck,2);
%for nst=1:nstmax   % adapt latini and longini according to cropped zone  and moove this calcul to rcv
%    linecheck(nst)=round(-(latcheck(nst)-latini)/dlat)+1; %+1; %+1 because in matlab lines and col start by the index "1"??
%    colcheck(nst)=round((longcheck(nst)-longini)/dlong)+1; %+1 because in matlab lines and col start by the index "1"??        
%end
  
%%SAVING OPTIONS        ?!!!!!!!!!!! exclusive for now (1 of the 3), because of transpositions required before saving, to adapt
visu=0
writegher=0
writemumm=0

%SCALE AND DISPLAY OPTIONS
% video
mkmovies=0
fps=2;
% graphs
showinlog10=1       % =1 if... %adapt so as conditionnal 
climdat=1           % add option
visunodat=0
formatimg='png'     % add in the code and in the subroutine v2D
% a-priori scale limits according to parameter and data  !!!! update + check uniform
% in rcv
if param=='CHL'
    if workedinlog==1
        climindata=0
        climaxdata=30
    else
        climindata=0.1  %to erase ..
        climaxdata=2  %to erase and also for TSM...
    end
    if showinlog10==1
        climindata=-1
        climaxdata=2
    end
    climinbckd=0                    %adapt?
    climaxbckd=30                   % adapt?
    climineof1=-15
    climaxeof1=15
    climineof=-15
    climaxeof=15
end
if param=='TSM'
    if workedinlog==1
        climindata=0
        climaxdata=50
    else      
        climindata=-2  % to adapt
        climaxdata=2  % to adapt
    end    
    if showinlog10==1
        climindata=-1
        climaxdata=2
    end
    if sensor=='A'
        climinbckd=0                    % adapt?
        climaxbckd=30                   % adapt?
        climineof1=-20
        climaxeof1=20
        climineof=-15
        climaxeof=15
    elseif sensor=='M'
        climinbckd=0                    % adapt?
        climaxbckd=30                   % adapt?
        climineof1=-10
        climaxeof1=10
        climineof=-8
        climaxeof=8
    end
end

readRebuildini=1