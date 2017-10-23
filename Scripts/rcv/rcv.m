oneD22D_started=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% below is the program itself       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

temprebuild=1  % initialisation of counter for rebuilding
if workedinanomaly==1
    moydata=load(filenameinmeandata);
end
useneof=load(filenameinneofretained);
showneof=useneof; % change show by useneof in the rest of the code, or make distinction

%useneof=6; showneof=6; % here special

if mkmovies==1
    if oridata==1
    fig=figure;
    set(fig,'DoubleBuffer','on');
    set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off')
    filenameout=[sensor,'_',param,'_',period,'_ori.avi'];
    mov1 = avifile(filenameout,'fps',fps,'quality',100,'compression','None')
    end
    if filleddata==0    
    fig=figure;
    set(fig,'DoubleBuffer','on');
    set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off')
    filenameout=[sensor,'_',param,'_',period,'_f.avi']    
    mov2 = avifile(filenameout,'fps',fps,'quality',100,'compression','None')
    end
    if rebuild==1
    fig=figure;
    set(fig,'DoubleBuffer','on');
    set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off')
    filenameoutmov3=[sensor,'_',param,'_',period,'_fi.avi'];
    mov3 = avifile(filenameout,'fps',fps,'quality',100,'compression','None')
    end
end

for i=0:9
    datemonth(i+1,:)=[int2str(0) int2str(i)];
    dateday(i+1,:)=[int2str(0) int2str(i)];
    datehh(i+1,:)=[int2str(0) int2str(i)];
    datemin(i+1,:)=[int2str(0) int2str(i)];
end
for i=10:12
    datemonth(i+1,:)=[int2str(i)];
end
for i=10:31
    dateday(i+1,:)=[int2str(i)];
end
for i=10:24
    datehh(i+1,:)=[int2str(i)];
end
for i=10:60
    datemin(i+1,:)=[int2str(i)];
end

%OPENING AND CREATING MASKS (preparing visu of temporal elimination)
mask=greadsingle(filenameinMasktempred);         %%%%%%%%%modified to single
mask=round(mask);
openedmaskred=1
maskori=greadsingle(filenameinMask);      %%%%%%%%%modified to single
maskori=round(maskori);
openedmaskori=1
%maskorir=imrotate((maskori'),90);
maskorir=rot90(maskori');  %octave modif    !!!!!!!!!!+ modif to load gis coast mask instead 

%here only gis coast no mask from it:
    coast=load('17876.dat'); % belcolour region
   % coast=[1 1]; %en attendant le bon fichier
    %coast=load('22622.dat'); % Corsica region
%if plotcoastusgs==1
%end

%[m,n,k]=size(mask)
[m,n]=size(mask) %octave modif
notmattpred=not(mask);

%STATISTICS ON INPUT DATA
if vistempelim==1
    visutpelim=maskori+mask;
   if writegher==1
        filenameout=['maskpixelim.gher'];
        gwrite(filenameout,visutpelim);
   end
   if writemumm==1
            %visutpelimb=visutpelim';
            filenameout=['maskpixelim.bin'];
            fidD = fopen(filenameout,'w','ieee-le'); 
            count = fwrite(fidD,visutpelim,'float32');
            status = fclose(fidD);            
   end
   if visu==1
     figure
     imagesc(visutpelim); 
     %xlabel('longitude, degree east','FontWeight','bold','FontSize',14) %here special
     %ylabel('latitude, degree north','FontWeight','bold','FontSize',14) %here special
     title(['Mask of pixels considered for the analysis of ',sensor,' ',param,' ',period],'FontWeight','bold','FontSize',16)
     filenameout=['visutempelim'];
     saveas(gcf,filenameout,formatimg)
   end
end

if absentdata==1
    nodat1=greadsingle(filenameinNodata1);  
    nodat2=greadsingle(filenameinNodata2);  
    minnodat1=min(min(nodat1))
    maxnodat1=max(max(nodat1))
    meannodat1=mean(mean(nodat1))    
    minnodat2=min(min(nodat2))
    maxnodat2=max(max(nodat2))
    meannodat2=mean(mean(nodat2))
    if visu==1
     figure
     imagesc(nodat1);
     set(gca,'CLim',[60,100]); % add as option
     title(['map of temporal % of absent data for all images of ',sensor,' ',param,' ',period],'FontWeight','bold','FontSize',16)
     shading flat, colorbar; 
     filenameout=['vistpabsdatall'];
     saveas(gcf,filenameout,formatimg)
     figure
     imagesc(nodat2); 
     set(gca,'CLim',[60,100]); % add as option
     title(['map of temporal % of absent data for selected images of ',sensor,' ',param,' ',period],'FontWeight','bold','FontSize',16)
     shading flat, colorbar; 
     filenameout=['vistpabsdatspsel'];
     saveas(gcf,filenameout,formatimg)
   end
end

if background==1
    data=greadsingle(filenameinBackgd);       %%%%%%%%%modified to single
    openedbackground=1
    minbckgd=min(min(data))
    maxbckgd=max(max(data))
    %if workedinlog==1
    %    data=10.^(data);
    %end
    %data=permute(data,[1 3 2]);
    [pmax,qmax,k]=size(data)
    if AL==0
        if showinlog10==1
            for q=1:qmax
                data(:,q,1)=log10(data(:,q,1));
            end
        end
    end
    [indicei,indicej,sbis] = find(mask);  % move out of background calc as uniform for all the run outputs ?
    %size(indicei)  %
    %size(indicej)  %
    vectdat=(data(:,1)); 
    %size(vectdat)  %      
    % below is the rebuild to 2D / display / save module
    filenameout1=[sensor,'_',param,'_',period,'_backgrd'];
    title1=[sensor,' ',param,' ',unit,', background'];
    climin=climinbckd;
    climax=climinbckd;
    [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    

if extractserie==1
    % loading ccord from a file
        coordst=load(filenameinstations);
        [nstmax,ncoord]=size(coordst) ;
        latcheck=coordst(:,1);
        longcheck=coordst(:,2);
    % defining coord internally   
        %latcheck=[51.47,51.308,51.433,51.526,51.990]; %north for 
        %longcheck=[3.23,2.849,2.809,1.030,2.090]; %east
        %1: Scheelde turbidity maximum
        %2: station 230                  
        %3: station 330                   
        %4: NorthDogger           %%%%%not used in recolour focus :lat 55.683, long
        %2.280
        %5: OysterGrounds         %%%%%not used in recolour focus: lat 54.414, long
        %4.042
        %6: Warp                  
        %7: WestG   
        nstmax=size(latcheck,1);
        for nst=1:nstmax   % adapt latini and longini according to cropped zone  and moove this calcul to rcv
            linecheck(nst)=round(-(latcheck(nst)-latini)/dlat)+1; %+1; %+1 because in matlab lines and col start by the index "1"??
            colcheck(nst)=round((longcheck(nst)-longini)/dlong)+1; %+1 because in matlab lines and col start by the index "1"??        
        end
end

%LOADING TEMPORAL INFORMATION and preparing data for temporal interpolations
% opening vectdates and streaming to corresponding variable
vectseldates=greadsingle(filenameinvectdates);                     %%%%%%%%%modified to single
nbimax=size(vectseldates,1);
vsd1=sortrows(vectseldates,[2 3 4 5 6 7]); % en fait enough check on col 2 as it is nimused itself ; check
first=0
for ncountimjet=1:nbimax
      if vsd1(ncountimjet,2)==1
         if first==0
            nbimaxunused=ncountimjet-1
            first=1
        end
      end
end
vsdini=vsd1(nbimaxunused+1:nbimax,:);
%vsd2=vsd1(nbimaxunused+1:nbimax,:); % option subset virée, ici vsd2=vsdini
%check ok on rebuildnoi
nbimused=nbimax-nbimaxunused  %check

%preparing the x time labels for display of time series turn this as a
%subroutine for tempeof and for extracted series
nm=0;
nbt=0;
dm=12;
for yeartab=1997:2006
    for monthtab=1:dm:12
        nbt=nbt+1;
        vtick(nbt,1)=yeartab;
        vtick(nbt,2)=monthtab;
        vtick(nbt,3)=1;
    end
end
vtimetick=datenum([vtick(:,1) vtick(:,2) vtick(:,3)]);
zeros(1:nbimused)=0;
%datevector=[vsdini(:,3) vsdini(:,4) vsdini(:,5) vsdini(:,6) vsdini(:,7) zeros.']; 
serialdate=datenum([vsdini(:,3) vsdini(:,4) vsdini(:,5) vsdini(:,6) vsdini(:,7) zeros.']);
sdtpy=datenum([ [tpyi tpyf].' [1 1].' [1 1].' [0 0].' [0 0].' [0.0 0.0].' ]);

%bla=0

if oridata==1   
   
if workedinanomaly==1
    databckgrd=greadsingle(filenameinBackgd);
    [pb,kb,qb]=size(databckgrd);      
    [flag,imax,jmax,kmax,valex,nbmots,nl,ir,pastread,fidslice] = greadslice(filenameinAnomaly)
    openedOridat=1
    %[pmax,k,qmax]=size(data)  
else  
     data=greadsingle(filenameinDat); 
     [flag,imax,jmax,kmax,valex,nbmots,nl,ir,pastread,fidslice] = greadslice(filenameinDat)
     openedOridat=1
     %[pmax,k,qmax]=size(data)
end     
       
[indicei,indicej,sbis] = find(mask);
size(indicei)
size(indicej)
if norimax==0
    norimax=kmax;
end



    for q=norimin:norimax %qmax % qmax instead qmax for short test
        tab = ureadslice(fidslice,imax,jmax,kmax,valex,nbmots,nl,ir,pastread);
        if workedinanomaly==1       
                tab(1:imax,1)=tab(1:imax,1)+databckgrd(1:imax,1,1);
        end       
        if (workedinlog==1) & (showinlog10==0) 
            tab=10.^(tab);
        end
        if (workedinlog==0) & (showinlog10==1)  
            for q=1:imax
                tab(q,1)=log10(tab(q,1));
            end
        end
%        if showinlog10==1   
            %for p=1:pmax
            %    if imagep2d(p)<=0
           %        imagep2d(p)=0.1;
                    %imagep2d(p)=NaN;
                    %    end
                    %end
 %           data(1:pmax)=log10(data(1:pmax));
 %      end              
        
        pmax=imax;
        vectdat=tab; 
        sv=size(vectdat)
        
        if (q==norimin)
            vectdatmax(1:sv)=-100000;
            vectdatmin(1:sv)=100000;
            vectdatmax=vectdatmax.';
            vectdatmin=vectdatmin.';
        end
        %vectdatmax(vectdat > vectdatmax)=vectdat;
        %vectdatmin(vectdat < vectdatmin)=vectdat;
        
        mask=find(vectdat > vectdatmax);
        vectdatmax(mask)=vectdat(mask);
        mask=find(vectdat < vectdatmin); 
        vectdatmin(mask)=vectdat(mask);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
        dateq2=[num2str(vsdini(q,5)),'/',num2str(vsdini(q,4)),'/',num2str(vsdini(q,3)),' at ',num2str(vsdini(q,6)),'h',num2str(vsdini(q,7)),'min']
        dateq3=[int2str(vsdini(q,3)),'_',datemonth(vsdini(q,4)+1,:),'_',dateday(vsdini(q,5)+1,:),'_',datehh(vsdini(q,6)+1,:),'_',datemin(vsdini(q,7)+1,:)]
        qq=num2str(q);
        %if visu==1   
            filenameout1=['unfilled/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];
            %else
            filenameout2=['unfilled/',sensor,'_',dateq3,'_',param,'_v0.1_ori00'];
            %end       
        title1=[sensor,' ',param,', date : ',dateq2,'             ',unit];
        climin=climindata
        climax=climaxdata
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        close all    
        if extractserie==1
            for nst=1:nstmax
               extractst(q,nst)=MOut(linecheck(nst),colcheck(nst)); % adapt to stations with boucle sur nst if not ok.
            end 
        end
    end
    fclose(fidslice); 
    %if mkmovies==1  %%adapt later
    %    mov1 = close(mov1);
    %end
    if extractserie==1   % adapt for multistations display according to the file of coordinates
       % figure
        %ext=10.^extractst(:,1:nstmax);
        plot(serialdate,10.^extractst(:,1:nstmax),'*','MarkerSize',3)
        set(gca,'XTick',vtimetick) %,'fontweight','bold')   %to adapt to 
        set(gca,'XTickLabel',datestr(vtimetick,20)) 
        set(gca,'Xlim',[sdtpy(1),sdtpy(2)])         
        %adapt to load stations name in same order
        %h = legend('Schelde turb. max. ','st. 230','st.
        %330','Warp','WestG');
%        h = legend('Schelde turb. max. ','st. 230','st. 330','NorthDogger','OysterGrounds','Warp','WestG');
        h = legend('G6','N2','N10','N20','W2','W20','W70','WAR');
%        set(h,'FontWeight','bold','FontSize',18);
        title(['Original data extracted at various stations; ',sensor,' ',param,' ',period],'FontWeight','bold','FontSize',20)
        filenameout=['unfilled/',sensor,'_',param,'_extracted_ori'];
        xlabel('time ','FontWeight','bold','FontSize',18) %here special
%        ylabel([param,' log10(',unit,')'], 'FontWeight','bold','FontSize',18) %here special
        ylabel([param,' ',unit], 'FontWeight','bold','FontSize',18) %here special
        set(gca,'FontWeight','bold','FontSize',18)
        saveas(gcf,filenameout,'jpg')
        %filenameout1=[filenameout,'.gher'];
        %gwrite(filenameout,extractst);
%        filenameout1=[filenameout];
        save(filenameout,'extractst','-ascii');        
    end
    
    %saving min max files and visu
    visutemp=visu; visu=1;
    writemummtemp=writemumm; writemumm=1;  
        mask=find(vectdatmax == -100000);   % oulai lo tsarikh etze arbaa
        vectdatmax(mask)=NaN;
        mask=find(vectdatmin== 100000); 
        vectdatmin(mask)=NaN;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
            filenameout1=['unfilled/',sensor,'_',param,'_v0.1_max'];
            %else
            filenameout2=['unfilled/',sensor,'_',param,'_v0.1_max'];
            %end       
        title1=['Map of maximum values',sensor,' ',param,' ',unit];
        climin=climindata
        climax=climaxdata
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdatmax,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
            filenameout1=['unfilled/',sensor,'_',param,'_v0.1_min'];
            %else
            filenameout2=['unfilled/',sensor,'_',param,'_v0.1_min'];
            %end       
        title1=['Map of minimum values, ',sensor,' ',param,' ',unit];
        climin=climindata
        climax=climaxdata
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdatmin,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        visu=visutemp;
        writemumm=writemummtemp;
end

if filleddata==1   % adapt all from oridat to double here
   
if workedinanomaly==1
    databckgrd=greadsingle(filenameinBackgd);
    [pb,kb,qb]=size(databckgrd);      
    [flag,imax,jmax,kmax,valex,nbmots,nl,ir,pastread,fidslice] = greadslice(filenameinFilleddat)
    openedFilled=1
    %[pmax,k,qmax]=size(data)  
else  
     data=greadsingle(filenameinDat); 
     [flag,imax,jmax,kmax,valex,nbmots,nl,ir,pastread,fidslice] = greadslice(filenameinFilleddat)
     openedFilled=1
     %[pmax,k,qmax]=size(data)
end     
       
[indicei,indicej,sbis] = find(mask);
size(indicei)
size(indicej)
if norimax==0
    norimax=kmax;
end

    for q=norimin:norimax %qmax % qmax instead qmax for short test
        tab = ureadslicedouble(fidslice,imax,jmax,kmax,valex,nbmots,nl,ir,pastread);
        if workedinanomaly==1       
                tab(1:imax,1)=tab(1:imax,1)+databckgrd(1:imax,1,1)+moydata;
        end       
        if (workedinlog==1) & (showinlog10==0) 
            tab=10.^(tab);
        end
        if (workedinlog==0) & (showinlog10==1)  
            for q=1:imax
                tab(q,1)=log10(tab(q,1));
            end
        end
%        if showinlog10==1   
            %for p=1:pmax
            %    if imagep2d(p)<=0
            %        imagep2d(p)=0.1;
                    %imagep2d(p)=NaN;
                    %    end
                    %end
 %           data(1:pmax)=log10(data(1:pmax));
 %      end              
        
        pmax=imax;
        vectdat=tab; 
        sv=size(vectdat)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
        dateq2=[num2str(vsdini(q,5)),'/',num2str(vsdini(q,4)),'/',num2str(vsdini(q,3)),' at ',num2str(vsdini(q,6)),'h',num2str(vsdini(q,7)),'min']
        dateq3=[int2str(vsdini(q,3)),'_',datemonth(vsdini(q,4)+1,:),'_',dateday(vsdini(q,5)+1,:),'_',datehh(vsdini(q,6)+1,:),'_',datemin(vsdini(q,7)+1,:)]
        qq=num2str(q);
%        if visu==1   
            filenameout1=['filledDin/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];
            %       else
            filenameout2=['filledDin/',sensor,'_',dateq3,'_',param,'_v0.1_fil00'];
            %end       
        title1=[sensor,' ',param,', date : ',dateq2,'             ',unit];
        climin=climindata
        climax=climaxdata
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        close all        
    end
    fclose(fidslice); 
    %if mkmovies==1  %%adapt later
    %    mov1 = close(mov1);
    %end
end

%
if (errorfields==1)
    mask(:)=1 % or load transpmaskindineof
    [error_field,datarecOI]=dineof_errorfield(dataini,datarec,mask,crossval,mu2_eff,vlsng,lftvec,N)
    
end

% bla=1

if tempmodes==1      %add line zero or grid hztal + value and date clims from rcvini
    vlsng=load('outputEof.vlsng');
    temp=load('outputEof.rghvec'); 
    figure
    sizetemp1=size(temp,1)
    sizetemp2=size(temp,2)    
    if reverteof1==1 
        temp(:,revvect)=-temp(:,revvect); %to reverse first eof for user habit visualisation sense
    end
    plot(serialdate,temp(:,1:showneof2),'LineWidth',1)
    set(gca,'XTick',vtimetick) %,'fontweight','bold')   %to adapt to 
    set(gca,'XTickLabel',datestr(vtimetick,20)) 
    set(gca,'Xlim',[sdtpy(1),sdtpy(2)]) 
    %datetick('x',20)  %the data for the x axis must be serial date numbers
    %(as produced by datenum).
    %grid on
    h = legend('eof 1','eof 2','eof 3');%,'FontWeight','bold','FontSize',18);
    title(['Dominants temporal EOFS from ',sensorname,' ',param],'FontWeight','bold','FontSize',20)
%    title(['Temporal EOFS retained by the analysis of ',sensor,' ',param,' ',period],'FontWeight','bold','FontSize',20)
    filenameout=[sensor,'_',param,'_',period,'_tempeofs'];
    xlabel('time','FontWeight','bold','FontSize',16) %here special
    ylabel('temporal variation coefficient','FontWeight','bold','FontSize',16) %here special
    set(gca,'FontWeight','bold','FontSize',12)
    set(gca,'Ylim',[climtpeofmin,climtpeofmax])
    %datetick('x',20)  %the data for the x axis must be serial date numbers (as produced by datenum).
    %YGrid
    set(gca,'YGrid','on')
    saveas(gcf,filenameout,'jpg')
    %grid off
end

bla=3

if spatialmodes==1
    filenamein=['outputEof.lftvec']; %spatial modes 
    data=load(filenamein); 
    [pmax,qmax,k]=size(data)
    
    vlsng=load('outputEof.vlsng');   %%% no need here except for visu eof*vlsng in " sqrt(data unit)" .. to adapt
    for neof=1:useneof
        mvlsng(neof,neof)=vlsng(neof); 
    end
    
    %if showinlog10==1    %%% no need here except for visu eof*vlsng in log of " sqrt(data unit)"
    %    for q=1:qmax
    %        data(:,q,1)=log10(data(:,q,1));
    %    end
    % end
    
    if datareal2D==1  % set this out from this conditional loop as uniform for all rcv
        mask=mask.';
        [m,n,k]=size(mask);
        notmattpred=notmattpred.';
    end
    
    [indicei,indicej,sbis] = find(mask); 
    size(indicei)
    size(indicej)
    
    showinlog10bf=showinlog10;
    showinlog10=0; %forced here unless integration of vlsng
    
    if showneofsp==0
        showneofsp=useneof
    end
        
    for q=1:showneofsp              

        %first eof with singular values aready multiplied   
            if reverteof1==1 
                  speof=data; % option for *mvlsng(showneof,showneof)+moydata;??
                  speof(1:pmax,revvect)=-data(1:pmax,revvect);  
            else
                  speof=data; % option for *mvlsng(showneof,showneof)+moydata;??
            end  
        
        if eofmvlsng==1
            vectdat=(speof(:,q))*mvlsng(q,q);
            uniteof=' speof*vlsng';
        else
            vectdat=(speof(:,q));  
            uniteof=' coef.';
        end
        
        size(vectdat)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
        dateq=num2str(q);
        filenameout1=[sensor,'_',param,'_',period,'_eof',dateq];
        filenameout2=[sensor,'_',param,'_',period,'_eof',dateq];        
        %title1=['Spatial EOFs for ',sensor,' ',param,', ',period,', EOF n° ',dateq,uniteof];
        title1=['Spatial EOFs for ',sensorname,' ',param,', ',period,', EOF n° ',dateq,uniteof];
        if q==1
            climin=climineof1
            climax=climaxeof1
        else
            climin=climineof
            climax=climaxeof
        end
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        close all
    end
    showinlog10=showinlog10bf;
end

bla=44

if rebuild==1  
    data=load('outputEof.lftvec');
    [pmax,qmax,k]=size(data)
    vlsng=load('outputEof.vlsng');
    temp=load('outputEof.rghvec'); 
    done=1
    temp(1:2,1:2)
    [indicei,indicej,sbis] = find(mask); % move out up once and for all
    size(indicei)
    size(indicej) 
    
    for neof=1:useneof
            mvlsng(neof,neof)=vlsng(neof);  
    end
    
    %building of the interpolated temporal eofs with 3 subroutines  
    [nbjyear,nbjsince97,nbjpastmonth]=dayssince(referenceyear) ; % 
    [tini]=buildtini(vsdini,nbjyear,nbjsince97);
    [tint,tempintneof,ntintmax,ntintmin,vsd2]=interpeof(tini,temp,djint,nbjyear,nbjsince97,useneof,ntintmax,ntintmin,rebuildmissim);
               
     if workedinanomaly==1
            databckgrd=gread(filenameinBackgd); 
     else
            databckgrd=0;
     end
    
     nbperioday=1; 
     MOut2(1:m,1:n)=0;  
     count=0;

     % rebuilding the images from interpolated temporal eofs 
 
     for q=ntintmin:ntintmax       
        imagep2d(1:pmax)=databckgrd+moydata+data(1:pmax,1:useneof)*mvlsng(1:useneof,1:useneof)*tempintneof(1:useneof,q);
        
        if (workedinlog==1) & (showinlog10==0)    % check idem in rebuildnoi
           %for pix=1:pmax
                imagep2d(1:pmax)=10.^imagep2d(1:pmax);
           %end
        end   
        
        %%%%%%%%%%%%%%
        vectdat=imagep2d;
        %size(vectdat)   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
        qq=num2str(q);
        %if q<10    % to add one "0" for ordered display in the
        %first 10 images
        %    zerosi='0';
        %else
        %     zerosi='';
        % end 
        dateq2=[num2str(vsd2(q,5)),'/',num2str(vsd2(q,4)),'/',num2str(vsd2(q,3)),' at ',num2str(vsd2(q,6)),'h',num2str(vsd2(q,7)),'min']
        dateq3=[int2str(vsd2(q,3)),'_',datemonth(vsd2(q,4)+1,:),'_',dateday(vsd2(q,5)+1,:),'_',datehh(vsd2(q,6)+1,:),'_',datemin(vsd2(q,7)+1,:)]
        %if visu==1   
        %filenameout1=['filled_regtpstep/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];
        filenameout1=['filled_regtpstep/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];     % produce separately so as reb2d can produce both visu and file       
        %else
        filenameout2=['filled_regtpstep/',sensor,'_',dateq3,'_',param,'_v0.1_fil03'];
        %end       
        title1=[sensor,' ',param,', date : ',dateq2,'             ',unit];
        climin=climindata
        climax=climaxdata
        visu=0; %forced here
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,filenameout2,title1,showinlog10,climdat,climin,climax,maskorir);
        visu=1; %reset here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
        if mkmovies==1   % check if still ok after modifications
                    F = getframe(gcf);
                    mov3 = addframe(mov3,F);
                    %mov3 = getframe(gcf);
        end
        close all
        
        % extract value at cheking stations
        % adapt perioday : parameter in days ,  to a variable number of days/month if monthly mean requested
        
        if tempaverage==1
            MOut2=MOut2+MOut;
            count=count+1;
           % if (tint(q)+perioday+0.5)>(nbperioday*perioday+tint(ntintmin)) 
            if (count+0.5)>(perioday/djint) 
                MOut2=MOut2/perioday;
                if extractserie==1
                   for nst=1:nstmax
                       extractst(nbperioday,nst)=MOut2(linecheck(nst),colcheck(nst)); % adapt to stations with boucle sur nst if not ok.
                   end 
                end
                %save multitemporal mean fields
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                matrixOut=MOut2;
                sperioday=num2str(perioday);
                snbperioday=num2str(nbperioday);    
                filenameout2=['filled_tpaverage/',sensor,'_',sperioday,'_',snbperioday,'_',param,'_v0.1_fil03'];
                filenameout1=['filled_tpaverage/',snbperioday,'f_',sensor,'_',param,'_',period,'_avperiodn_',snbperioday];  
                title1=[sensor,' ',param,' ',sperioday,'days averages, period n°: ',snbperioday];
                climin=climindata; % adapt for different values then for instantaneous fields
                climax=climaxdata; % adapt for different values then for instantaneous fields
                if writegher==1
                    filenameout=[filenameout2,'.gher'];
                    gwrite(filenameout,matrixOut);
                end
                if writemumm==1
                    matrixOutm=matrixOut';
                    filenameout=[filenameout2,'.bin'];
                    fidD = fopen(filenameout,'w','ieee-le'); 
                    count = fwrite(fidD,matrixOutm,'float32');
        	        status = fclose(fidD);  	              
                end
                if visu==1  
                    %matrixOut=imrotate((matrixOut'),90);
                    matrixOut=imrotate((matrixOutm),90);
                        b=longini+dlong*[1:size(matrixOut,2)]; % moove out of subroutine as uniform for all the run reconstructions?
                        a=(latini-dlat*size(matrixOut,1))+dlat*[1:size(matrixOut,1)]; %oove out as uniform for all the run reconstructions?
                        figure
                        pcolor(b,a,matrixOut);  
                        title(title1,'FontWeight','bold','FontSize',14)
                        xlabel('longitude, degree east','FontWeight','bold','FontSize',14) %here special
                        ylabel('latitude, degree north','FontWeight','bold','FontSize',14) %here special
                        shading flat, colorbar; 
                        if showinlog10==1
                            set(gca,'CLim',[-1,2]) 
                            set(colorbar,'YTick',[-1:1:2],'YTickLabel',[0.1;1;10;100],'FontWeight','bold','FontSize',14)
                        else
                            if climdat==1
                                set(gca,'CLim',[climin,climax]) %here special
                            end
                            set(colorbar,'FontWeight','bold','FontSize',14)
                        end
                        hold on
                        %contour(b,a,maskorir,[0.5,0.5],'-k')
                        plot(coast(:,1),coast(:,2),'k','linewidth',0.5); %here special
                        %display buoys
                        hold on
                        longstations=[3.23,2.090]; %  adapt from rcvini longcheck
                        latstations=[51.47,51.990];
                        scatter(longstations,latstations,12,'k','filled')                       
                        %colorbar;
                        saveas(gcf,filenameout1,formatimg)
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                MOut2(:,:)=0;
                nbperioday=nbperioday+1;
                count=0;
            end      
        end
        close all
    end
    nbperioday=nbperioday-1; 
    if mkmovies==1
        %movie2avi(mov3,filenameoutmov3,'fps',2,'quality',100,'compression','None')
        mov3 = close(mov3);
    end
    if extractserie==1   % adapt for multistations display according to the file of coordinates
       % figure
        plot(extractst(1:nbperioday,1:nstmax),'LineWidth',3)
        h = legend('Schelde turb. max. ','st. 230','st. 330','Warp','WestG');
%        h = legend('Schelde turb. max. ','st. 230','st. 330','NorthDogger','OysterGrounds','Warp','WestG');
        set(h,'FontWeight','bold','FontSize',18);
        title([sperioday,' days mean value at various stations; ',sensor,' ',param,' ',period],'FontWeight','bold','FontSize',20)
        filenameout=[sensor,'_',param,'_',sperioday,'_days_meanst'];
        xlabel('time (week of the year)','FontWeight','bold','FontSize',18) %here special
        ylabel([param,' log10(',unit,')'], 'FontWeight','bold','FontSize',18) %here special
        set(gca,'FontWeight','bold','FontSize',18)
        saveas(gcf,filenameout,'jpg')
        filenameout1=[filenameout,'.gher'];
        gwrite(filenameout,extractst);
    end
end


if rebuildmissim==1  
    data=load('outputEof.lftvec');
    [pmax,qmax,k]=size(data)
    vlsng=load('outputEof.vlsng');
    temp=load('outputEof.rghvec'); 
    done=1
    temp(1:2,1:2)
    [indicei,indicej,sbis] = find(mask); 
    size(indicei)
    size(indicej) 
    
    for neof=1:useneof
            mvlsng(neof,neof)=vlsng(neof);  
    end
    
    %building of the interpolated temporal eofs with 3 subroutines  
    [nbjyear,nbjsince97,nbjpastmonth]=dayssince(referenceyear) ; % 
    [tini]=buildtini(vsdini,nbjyear,nbjsince97);
    
     %here special, adapt later
     ntintmin=1;
     ntintmax=2;
    
    [tint,tempintneof,ntintmax,ntintmin,vsd2]=interpeof(tini,temp,djint,nbjyear,nbjsince97,useneof,ntintmax,ntintmin,rebuildmissim);
               
     if workedinanomaly==1
            databckgrd=gread(filenameinBackgd); 
     else
            databckgrd=0;
     end
     

     
     % rebuilding the images from interpolated temporal eofs  
     for q=ntintmin:ntintmax       
        imagep2d(1:pmax)=databckgrd+moydata+data(1:pmax,1:useneof)*mvlsng(1:useneof,1:useneof)*tempintneof(1:useneof,q);
        
        if (workedinlog==1) & (showinlog10==0)    % check idem in rebuildnoi
           %for pix=1:pmax
                imagep2d(1:pmax)=10.^imagep2d(1:pmax);
           %end
        end   
        
        %%%%%%%%%%%%%%
%        vectdat=imagep2d.';    %%?test transposed?
        vectdat=imagep2d; 
        %size(vectdat)   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
        qq=num2str(q);
        %if q<10    % to add one "0" for ordered display in the
        %first 10 images
        %    zerosi='0';
        %else
        %     zerosi='';
        % end 
        dateq2=[num2str(vsd2(q,5)),'/',num2str(vsd2(q,4)),'/',num2str(vsd2(q,3)),' at ',num2str(vsd2(q,6)),'h',num2str(vsd2(q,7)),'min']
        dateq3=[int2str(vsd2(q,3)),'_',datemonth(vsd2(q,4)+1,:),'_',dateday(vsd2(q,5)+1,:),'_',datehh(vsd2(q,6)+1,:),'_',datemin(vsd2(q,7)+1,:)]
        if visu==1   
            %filenameout1=['filled_regtpstep/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];
            filenameout1=['missim/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];     % produce separately so as reb2d can produce both visu and file       
        else
            filenameout1=['missim/',sensor,'_',dateq3,'_',param,'_v0.1_fil03'];
        end       
        title1=[sensor,' ',param,', date : ',dateq2,'             ',unit];
        climin=climindata
        climax=climaxdata
        %display('rebelote')
        %min(indicei)
        %max(indicei)
        %min(indicej)
        %max(indicej)
        %size(vectdat)
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all
    end
end


bla=5

if rebuildnoi==1 
    rebuildnoi 
    data=load('outputEof.lftvec');
    [pmax,qmax,k]=size(data)
    vlsng=load('outputEof.vlsng');
    temp=load('outputEof.rghvec'); 
    [indicei,indicej,sbis] = find(mask); 
    size(indicei)
    size(indicej) 
    for neof=1:showneof
        mvlsng(neof,neof)=vlsng(neof); 
    end
    tempneof=temp.';
    if nnoimax==0
         nnoimax=size(tempneof,2)
    end
    if workedinanomaly==0
        databckgrd=0;
    else
        databckgrd=greadsingle(filenameinBackgd); 
        %sizebckgrd=size(databckgrd)
        %pmax
    end    
    for q=nnoimin:nnoimax %1 instead qmax for short test   
        
        qq=num2str(q);
        %if q<10    % to add one "0" for ordered display in the
        %first 10 images
        %    zerosi='0';
        %else
        %     zerosi='';
        % end 
        
        if checkcontrib==1
                climdatbf=climdat;
                showinlog10bf=showinlog10;
                climdat=0; showinlog10=0; % forced here as no sense to consider contribution as value yet as it has to be added together (with neg values) before logdisplay
            for neof=1:useneof
                neofstr=num2str(neof)
                vectdat=data(1:pmax,neof)*mvlsng(neof,neof)*tempneof(neof,q); 
                %size(vectdat)
                %size(indicei)
                %size(indicej)
                %mincontrib=min(vectdat)
                %maxcontrib=max(vectdat)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                % below is the module to rebuild to 2D / display / save         % adapt for beeing robust to other data preprocessing options (direct/logdirect/anomalydirect/...) 
                if (workedinlog==1) & (showinlog10==0)     
                    %for pix=1:pmax
               %         imagep2d(1:pmax)=10.^imagep2d(1:pmax);
                    %end
                end 
                filenameout1=['contrib_eof_',neofstr,'toimage_',qq];
                title1=['contrib of eof n° ',neofstr,' to image n° ',qq,' anomaly(log10(data))' ];  % unit?? as it is not yet direct since it should be added to background+meananomaly before taking the exp10 to be in original unit
                climin=1; % whatever, not used here
                climax=1;
                %if showinlog10==0
                %    if neof==1
                %        climin=climineof1;
                %        climax=climaxeof1;
                %    else
               %         climin=climineof;
               %         climax=climaxeof;
               %     end
               %  else
               %     climin=-1;
               %     climax=2;
               % end
                [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,title1,showinlog10,climdat,climin,climax,maskorir);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                close all 
            end
            climdat=climdatbf;
            showinlog10=showinlog10bf; 
        end
        
        %minim=min(databckgrd)
        %maxim=max(databckgrd)
        %moydata
        rebuildat=data(1:pmax,1:showneof)*mvlsng(1:showneof,1:showneof)*tempneof(1:showneof,q);
        %minim=min(rebuildat)
        %maxim=max(rebuildat)         
        imagep2d(1:pmax)=databckgrd+moydata+rebuildat;
        
        if (workedinlog==1)  & (showinlog10==0)     %to test + adapt with AL etc..
           %for pix=1:pmax
                imagep2d(1:pmax)=10.^imagep2d(1:pmax);
           %end
        end   
        
        %%%%%%%%%%%%%%%%%%%%%
        %% plot histog data
        %  [n, xout] = hist(imagep2d,20) ;
        %  figure
        %  bar(xout,n)
        %  minim=min(imagep2d)
        %  maxim=max(imagep2d)
        %%%%%%%%%%%%%%%%%%%%%
        
        if AL==0
            if showinlog10==1
                for p=1:pmax
                    if imagep2d(p)<=0
                        %imagep2d(p)=0.1;
                        imagep2d(p)=NaN;
                    end
                end
                imagep2d(1:pmax)=log10(imagep2d(1:pmax));
            end
        end
           
        vectdat=imagep2d; 
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % below is the module to rebuild to 2D / display / save 
        dateq2=[num2str(vsdini(q,5)),'/',num2str(vsdini(q,4)),'/',num2str(vsdini(q,3)),' at ',num2str(vsdini(q,6)),'h',num2str(vsdini(q,7)),'min']
        dateq3=[int2str(vsdini(q,3)),'_',datemonth(vsdini(q,4)+1,:),'_',dateday(vsdini(q,5)+1,:),'_',datehh(vsdini(q,6)+1,:),'_',datemin(vsdini(q,7)+1,:)]
        if visu==1   
            filenameout1=['filled/',qq,'f_',sensor,'_',param,'_',period,'_date_',dateq3];
        else
            filenameout1=['filled/',sensor,'_',dateq3,'_',param,'_v0.1_fil03'];
        end       
        title1=[sensor,' ',param,', date : ',dateq2,'             ',unit];
        climin=climindata
        climax=climaxdata
        [MOut,flag]=reb2d(coast,longini,dlong,latini,dlat,indicei,indicej,vectdat,m,n,pmax,notmattpred,writegher,writemumm,visu,filenameout1,title1,showinlog10,climdat,climin,climax,maskorir);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        if visu==1  %  test to see if ok after modifications
             if mkmovies==1
                    F = getframe(gcf);
                    mov3 = addframe(mov3,F);
                    %mov3 = getframe(gcf);
             end
             close all
        end 
    end
    if mkmovies==1
        %movie2avi(mov3,filenameoutmov3,'fps',2,'quality',100,'compression','None')
        mov3 = close(mov3);
    end
end


if displayby2==1
    
 visu=1
 lines=445 %to adapt
 col=644   
 nbtimecount=0;
 nbpage=0;
 %nbpage=50; % here forced test to finish CHL
 if norimax==0
     [norimax,blabla]=size(vsdini);
 end
 
    for q=norimin:norimax %qmax % qmax instead qmax for short test               
        
        if nbtimecount==0
            figure
            get(gcf,'Position')
            set(gcf,'Position',[232   010   560   600])
            niminplot=1;
            nbtimecount=1;
        elseif nbtimecount==nbtimeperfig
            nbpage=nbpage+1;
            qq=num2str(nbpage);
%            filenameout1=['compar/',qq,'c_',sensor,'_',param,'_',period,'_date_',dateq33];
            filenameout1=[qq,'c_',sensor,'_',param,'_',period,'_date_',dateq33];
            saveas(gcf,filenameout1,formatimg)
            close all
            figure
            set(gcf,'Position',[232   010   560   600])
            niminplot=1;
            nbtimecount=1;
        else
            nbtimecount=nbtimecount+1;
        end
        
        %dateq2=[num2str(vsdini(q,5)),'/',num2str(vsdini(q,4)),'/',num2str(vsdini(q,3)),' at ',num2str(vsdini(q,6)),'h',num2str(vsdini(q,7)),'min']
        dateq22=[num2str(vsdini(q,5)),'/',num2str(vsdini(q,4)),'/',num2str(vsdini(q,3))];
        dateq3=[int2str(vsdini(q,3)),'_',datemonth(vsdini(q,4)+1,:),'_',dateday(vsdini(q,5)+1,:),'_',datehh(vsdini(q,6)+1,:),'_',datemin(vsdini(q,7)+1,:)]
        dateq33=[int2str(vsdini(q,3)),'_',datemonth(vsdini(q,4)+1,:),'_',dateday(vsdini(q,5)+1,:)];
        qq=num2str(q);

        %plotting unfilled data
        filenamein1=['unfilled/',sensor,'_',dateq3,'_',param,'_v0.1_ori00.bin'];
        fidD = fopen(filenamein1,'r','ieee-le'); 
        [matrixOut,countD] = fread(fidD,[col,lines],'float32');
        fclose(fidD);
        %matrixOut=matrixOut.' % to build a subroutine visu 
        title1=[sensor,' ',param,', date : ',dateq22,'   original      ',unit];
        
        if visu==1
               % matrixOut=imrotate((matrixOut'),90);
                matrixOut=imrotate((matrixOut),90);
                b=longini+dlong*[1:size(matrixOut,2)]; % moove out of subroutine as uniform for all the run reconstructions?
                a=(latini-dlat*size(matrixOut,1))+dlat*[1:size(matrixOut,1)]; %oove out as uniform for all the run reconstructions?
                %figure
                subplot(nbtimeperfig,2,niminplot)
                pcolor(b,a,matrixOut); 
                %set(gca,'FontSize',8)
%                title(title1,'FontWeight','bold','FontSize',14)
                title(title1,'FontSize',8)
                %axes('FontSize',8)
                %xlabel('longitude, degree east','FontWeight','bold','FontSize',14) %here special
                %ylabel('latitude, degree north','FontWeight','bold','FontSize',14) %here special
                shading flat, colorbar; 
                if showinlog10==1
                    set(gca,'CLim',[-1,2]) 
                    set(colorbar,'YTick',[-1:1:2],'YTickLabel',[0.1;1;10;100],'FontSize',8)
                else
                    if climdat==1
                        set(gca,'CLim',[climin,climax]) %here special
                    end
                    set(colorbar,'FontSize',8)
                end
                hold on
                contour(b,a,maskorir,[0.5,0.5],'-k')
                niminplot=niminplot+1;
                %formatimg='png';
                %saveas(gcf,filenameout1,formatimg)
        end
        %plotting filled data
        filenamein2=['filledDin/',sensor,'_',dateq3,'_',param,'_v0.1_fil00.bin'];
        fidD = fopen(filenamein2,'r','ieee-le'); 
        [matrixOut,countD] = fread(fidD,[col,lines],'float32');
        fclose(fidD);
        title1=[sensor,' ',param,', date : ',dateq22,'     filled     ',unit];
        
        if visu==1
                matrixOut=imrotate((matrixOut),90);
                %matrixOut=imrotate((matrixOut).',90);
                b=longini+dlong*[1:size(matrixOut,2)]; % moove out of subroutine as uniform for all the run reconstructions?
                a=(latini-dlat*size(matrixOut,1))+dlat*[1:size(matrixOut,1)]; %oove out as uniform for all the run reconstructions?
                %figure
                subplot(nbtimeperfig,2,niminplot)
                pcolor(b,a,matrixOut);  
                %set(gca,'FontSize',8)
                %title(title1,'FontWeight','bold','FontSize',14)
                title(title1,'FontSize',8)

                %xlabel('longitude, degree east','FontWeight','bold','FontSize',14) %here special
                %ylabel('latitude, degree north','FontWeight','bold','FontSize',14) %here special
                shading flat, colorbar; 
                if showinlog10==1
                    set(gca,'CLim',[-1,2]) 
                    set(colorbar,'YTick',[-1:1:2],'YTickLabel',[0.1;1;10;100],'FontSize',8)
                else
                    if climdat==1
                        set(gca,'CLim',[climin,climax]) %here special
                    end
                    set(colorbar,'FontSize',8)
                end
                hold on
                contour(b,a,maskorir,[0.5,0.5],'-k')
                niminplot=niminplot+1;
                %formatimg='png';
                %saveas(gcf,filenameout1,formatimg)
        end
        
        %adapt below to display last needed??
        if (nbtimecount==nbtimeperfig) & (nbpage+1==floor((norimax-norimin+1)/nbtimeperfig))  
            nbpage=nbpage+1;
            qq=num2str(nbpage);
            formatimg='png';
            filenameout1=[qq,'c_',sensor,'_',param,'_',period,'_date_',dateq33];
            saveas(gcf,filenameout1,formatimg)
            close all
            figure
            niminplot=1;
            nbtimecount=1;
        end
        
          %adapt below to display last  needed ??
 
%            nbpage=nbpage+1;
%            qq=num2str(nbpage);
%            formatimg='png';
%            filenameout1=['compar/',qq,'c_',sensor,'_',param,'_',period,'_date_',dateq33];
%            saveas(gcf,filenameout1,formatimg)
            %close all
   end
    
end

if readplotexts==1
%plot extracted series


%preparing the x time labels for display of time series turn this as a
%subroutine for tempeof and for extracted series
nm=0;
nbt=0;
dm=6;
for y=1997:2006
    for m=1:dm:12
        nbt=nbt+1;
        vtick(nbt,1)=y;
        vtick(nbt,2)=m;
        vtick(nbt,3)=1;
    end
end
vtimetick=datenum([vtick(:,1) vtick(:,2) vtick(:,3)]);


    %building of the interpolated temporal eofs with 3 subroutines  later
    %load from saved tint, and adapt timing from dec days to datenum units
    %and reference
    temp=load('outputEof.rghvec'); 
    [nbjyear,nbjsince97,nbjpastmonth]=dayssince(referenceyear) ; % 
    [tini]=buildtini(vsdini,nbjyear,nbjsince97);
    [tint,tempintneof,ntintmax,ntintmin,vsd2]=interpeof(tini,temp,djint,nbjyear,nbjsince97,useneof,ntintmax,ntintmin,rebuildmissim);
    tmean=tint(ntintmin+4:7:ntintmax); % here specific to 7 days means
    tmeanadd=datenum([1997 1 1 0 0 0]);
    serialdate=tmean+tmeanadd;
    
%zeros(1:nbimused)=0;
%serialdate=datenum([vsdini(:,3) vsdini(:,4) vsdini(:,5) vsdini(:,6) vsdini(:,7) zeros.']);

sdtpy=datenum([ [tpyi tpyf].' [1 1].' [1 1].' [0 0].' [0 0].' [0.0 0.0].' ]);


        sensor2='MERIS';
        sensor2='MODIS';
        sperioday=num2str(7);

        bothparamin1graph=0;
        if bothparamin1graph==1
        
        param='TSM';
        filenamein=[sensor,'_',param,'_',sperioday,'_days_meanst'];  % stream from rcvini
        extractst=gread(filenamein);
        [nbperioday,nstmax]=size(extractst);      
        extractst=10.^(extractst);
        
        figure 
        plot(serialdate,extractst(1:nbperioday,1),'-b','LineWidth',1)
        hold on
        plot(serialdate,extractst(1:nbperioday,5),'-g','LineWidth',1) 
        
    %    also=0
    %if also==1
        param='CHL';        
        filenamein=[sensor,'_',param,'_',sperioday,'_days_meanst'];  % stream from rcvini
        extractst=gread(filenamein);
        [nbperioday,nstmax]=size(extractst);      
        extractst=10.^(extractst);
        hold on
        plot(serialdate,extractst(1:nbperioday,1),'-r','LineWidth',1)
        hold on
        plot(serialdate,extractst(1:nbperioday,5),'-m','LineWidth',1) 
    %end
        
        set(gca,'XTick',vtimetick) %,'fontweight','bold')   %to adapt to 
        set(gca,'XTickLabel',datestr(vtimetick,20)) 
        set(gca,'Xlim',[sdtpy(1),sdtpy(2)]) 
        
        
        param='TSM and CHL'; 
        h = legend('TSM at Schelde Turbidity Max. ','TSM at CEFAS buoy WestG','CHL at Schelde Turbidity Max. ','CHL at CEFAS buoy WestG');
%        h = legend('Schelde turb. max. ','st. 230','st. 330','NorthDogger','OysterGrounds','Warp','WestG');
        set(h,'FontWeight','bold','FontSize',12);
        title(['Weekly averaged TSM and CHL time series'],'FontWeight','bold','FontSize',20)
        filenameout=[sensor,'_',param,'_',sperioday,'_days_meanst'];
        xlabel('time','FontWeight','bold','FontSize',18) %here special
    %    ylabel([param,' log10(',unit,')'], 'FontWeight','bold','FontSize',18) %here special
        ylabel('TSM (mg/l) and CHL (µg/l)', 'FontWeight','bold','FontSize',18) %here special
        set(gca,'FontWeight','bold','FontSize',18)
        saveas(gcf,filenameout,formatimg)
        %filenameout1=[filenameout,'.gher'];
        %gwrite(filenameout,extractst);   

            
            
        end
        
        
        
        filenamein=[sensor,'_',param,'_',sperioday,'_days_meanst'];  % stream from rcvini
        extractst=gread(filenamein);
        [nbperioday,nstmax]=size(extractst);      
        extractst=10.^(extractst);
        
        figure 
        plot(serialdate,extractst(1:nbperioday,1:nstmax),'LineWidth',1)     
        set(gca,'XTick',vtimetick) %,'fontweight','bold')   %to adapt to 
        set(gca,'XTickLabel',datestr(vtimetick,20)) 
        set(gca,'Xlim',[sdtpy(1),sdtpy(2)]) 
        
        
     
        h = legend('Schelde Turbidity Max. ','st. 230','st. 330','Warp','WestG');
%        h = legend('Schelde turb. max. ','st. 230','st. 330','NorthDogger','OysterGrounds','Warp','WestG');
        set(h,'FontWeight','bold','FontSize',14);
        title([sperioday,' days mean value at various stations; ',sensor2,' ',param],'FontWeight','bold','FontSize',20)
        filenameout=[sensor,'_',param,'_',sperioday,'_days_meanst'];
        xlabel('time','FontWeight','bold','FontSize',18) %here special
    %    ylabel([param,' log10(',unit,')'], 'FontWeight','bold','FontSize',18) %here special
        ylabel([param,' ',unit], 'FontWeight','bold','FontSize',18) %here special
        set(gca,'FontWeight','bold','FontSize',18)
        saveas(gcf,filenameout,formatimg)
        %filenameout1=[filenameout,'.gher'];
        %gwrite(filenameout,extractst);   
        
end
        
        
finished_rcv=1
