%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter definition file %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add simplification by streaming info from preprecol through files to read, then moove these param towards oneD22Dv25 itself

clear all
close all

vsubset=[1:1:3];

%filenames
filenameinNodata1=['vistpabsdatall.gher'];
filenameinNodata2=['vistpabsdatspsel.gher'];
filenameinMask=['belcolour_subregion_mask.gher'];
filenameinBackgd=['2Dbelcolour_region_period_background.gher'];
filenameinAnomaly=['2Dbelcolour_region_period_anomaly.gher'];
filenameinMasktempred=['belcolour_region_tempreducedmask.gher'];
%filenameinFilleddat=['F2Dbelcolour_region_period_datfilledsub1TSM.gher']; %after dineof
filenameinFilleddat=['F2Dbelcolour_region_period_datfilled.gher']; %after dineof

%filenameinOridat1=['2Dbelcolour_region_period_dat.gher']; %before dineof
filenameinOridat1=['2Dbelcolour_region_period_datsub1CHL.gher']; %before dineof
filenameinOridat2=['2Dbelcolour_region_period_datsub2CHL.gher']; %before dineof
filenameinOridat3=['2Dbelcolour_region_period_datsub3CHL.gher']; %before dineof

filenameinvectdates=['vectseldates.gher']; 
filenameinmeandata=['meandata.val'];
filenameinneofretained=['neofretained.val'];

% options, add comments

datareal2D=0

%sensor, param, period, region
sensor='M'; %for MERIS
%sensor='A'; %for MODIS AQUA

%param='CHL';
param='TSM';
%param='SST';

%unit='(°C)';
unit='mg/l';
%unit='µg/l';

period='all';
subset=1;   %move /adapt /cancel when filled data adpated for run on 3subset at once...like oridat.
%region='all belcolour';

%stream lat long ini from preprecol.ini according to cropped zone!!!!
%for all belcolour focus
 %latini=60.002; 
 %longini=-4.000;
%for recolour focus
 latini=52.500; 
 longini=-4.000;

dlat=0.009;
dlong=0.014;

extractserie=0
latcheck=51.47; %north
longcheck=3.23; %east
linecheck=round(-(latcheck-latini)/dlat)+1; %+1; %+1 because in matlab lines and col start by the index "1"??
colcheck=round((longcheck-longini)/dlong)+1; %+1 because in matlab lines and col start by the index "1"??

%%%
%link both tp and sp eofs as modif realized to bring sign from 1st eof as 'positive habit'.
spatialmodes=1
tempmodes=1
showneof2=15 % replace in code showneof by useneof for reconstructions, and (showneof by showneof2)  everywhere to 'useneof'%%%
reverteof1=1
filleddata=0

oridata=0
climdat=1 %add option
%%%
vistempelim=0
visunodat=0
mkmovies=0
fps=2;

workedinlog=0   %to adapt with comments    =1 to rebuild and record data in real units if it was processed in log10; =0 to output in log10 if processed in log10% =0 if processed in raw units and need to be outputed in raw units.
showinlog10=0 % =1 if...
%%%

rebuildnoi=0  %without interpolation
checkcontrib=0 % displaying the contribution of each eof to the image reconstruction
nnoimin=1
nnoimax=0

rebuild=0  % with interpolation
ntintmin=8  %  % (if=0 then max number of images interp is taken)
ntintmax=8 %351  % (if=0 then max number of images interp is taken)

background=1  % rebuild, then display or save background

workedinanomaly=1    % option for rebuilding filled and original fields as anomaly only (workedinanomaly==0) or as real values ,with bckgrd addedd (workedinanomaly==1)
AL=1  %worked on the anomaly of the log10 of the data

djint=1; %delta days between rebuilt images from interpolated temporal eofs

%%output options, exclusive
visu=1
writegher=0
writemumm=0

%if rebuild==1
%    tempmodes=1
%    spatialmodes=1
%end

temprebuild=1  % initialisation of counter for rebuilding

moydata=load(filenameinmeandata);
useneof=load(filenameinneofretained);

%filenameinmeandata
%fidD = fopen(filenameinmeandata,'r','ieee-le'); 
%[moydata,countD] = fread(fidD,'float32');
%[moydata] = fread(fidD);

%fidD = fopen(filenameinneofretained,'r','ieee-le'); 
%[useneof,countD] = fread(fidD,'uint32');
%[useneof,countD] = fread(fidD);

showneof=useneof; % change show by useneof in the rest of the code, or make distinction

moydata
showneof

if param=='CHL'
    if workedinlog==1
        climindata=0
        climaxdata=30
    else
        climindata=0.1  %to erase ..
        climaxdata=2  %to errase and also for TSM...
        climineof1=-2
        climaxeof1=2
        climineof=-0.007
        climaxeof=0.007
    end
    if showinlog10==1
        climindata=-1
        climaxdata=2
    end
end
if param=='TSM'
    if workedinlog==1
        climindata=0
        climaxdata=50
    else      
        climindata=-2
        climaxdata=2
        climineof1=-0.004
        climaxeof1=0.004
        climineof=-0.003
        climaxeof=0.003
    end    
    if showinlog10==1
        climindata=-1
        climaxdata=2
    end
end

readRebuildini=1
