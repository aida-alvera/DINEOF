function [tint,tempintneof,ntintmax,ntintmin,vsd2]=interpeof(tini,temp,djint,nbjyear,nbjsince97,useneof,ntintmax,ntintmin,rebuildmissim)

% building the vector (tint) of number of decimal days passed since the 01/01/1997 for each interpolation time requested
% and building the daily interpolated eof

        jint1=floor(djint/2)+0.5  ;  % convention : reinterpolate at mid-day
        nbjtot9707=nbjsince97(12) ;
        tint=[jint1:djint:nbjtot9707] ;
        
        if rebuildmissim==1
            %atemissim=load('listelim.txt'); % improove later to load
            %vsdmiss from listelim.txt
            %[A,count] = fscanf(fid,format,size) 
            vsdmiss(1,:)=[0 1 2003 03 14 11 02];
            vsdmiss(2,:)=[0 1 2004 09 02 10 55];           
            [tint]=buildtini(vsdmiss,nbjyear,nbjsince97);
        end
            

        tini=tini';
        tint=tint';
        
        % recreate date info from tint
        ntint=size(tint)
        year=0; %ref for 1996=year 0  ;  if need to be modified for another period than 1997-2007, adapt in respect to subroutine dayssince
        month=0;
        for t=1:ntint
            vsd2(t,1)=t;                %counter
            vsd2(t,2)=1;                % adapt flag 1=used, 2=unused, 3=interpolated
            while tint(t)>nbjsince97(year+1)
                  year=year+1;
                  month=0;
            end
            vsd2(t,3)=year+1996;         
            while tint(t)>(nbjyear(month+1,year)+nbjsince97(year))
                  month=month+1;
            end
            vsd2(t,4)=month;     %month            
            vsd2(t,5)=(tint(t)+0.5)-nbjsince97(year)-nbjyear(month,year);    %day of the month
        end      
        vsd2(1:ntint,6)=12; % in respect to mid-day interpolation convention
        vsd2(1:ntint,7)=00; % in respect to mid-day interpolation convention
        
        
        if rebuildmissim==1
            vsd2=vsdmiss;           
        end
        
        for neof=1:useneof
            tempun=temp(:,neof);
            %size(tini)
            %size(tempun)   
            %size(tint)
            tempint=interp1(tini,tempun,tint);
            ntint=size(tempint,1);
            %tint=tint';   
            tempintneof(neof,1:ntint)=tempint(1:ntint).';
        end

        % detect first and last non NAN interpolated eofs
        firstnonnan=1;
        tempintused=isnan(tempint);
        for nonnan=1:ntint     % replace by while loop easier
            if tempintused(nonnan)==0
                if firstnonnan==1
                      ntintautomin=nonnan
                      firstnonnan=2;
                end
                ntintautomax=nonnan;
            end
        end
        
        % calculate ntintmin and ntintmax 
        if rebuildmissim==0
            
            if ntintmin==0
                ntintmin=ntintautomin;  %  take first non null interpolated image 
            end
            if ntintmax==0
                ntintmax=ntintautomax;   % take last non null interpolated image 
            else
                ntintmax=ntintautomin+ntintmax; %  take first non null interpolated image + number of images requested
            end
        
        end
        
 %%end subroutine 