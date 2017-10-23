function [nbjyear,nbjsince97,nbjpastmonth]=dayssince(referenceyear)


% replace by function ' eomday(Y,M) '
%number of days/month of each year of the period of the study

% modif to allow variable period and reference year  +/or  load number of
% days/month/year from a file made once 

%  referenceyear;  % (unused for now, parameter for later)

    %here for all regular years
        nbjpastmonth(1)=0;
        nbjpastmonth(2)=31;
        nbjpastmonth(3)=28;
        nbjpastmonth(4)=31;
        nbjpastmonth(5)=30;
        nbjpastmonth(6)=31;
        nbjpastmonth(7)=30;
        nbjpastmonth(8)=31;
        nbjpastmonth(9)=31;
        nbjpastmonth(10)=30;
        nbjpastmonth(11)=31;
        nbjpastmonth(12)=30;
        nbjpastmonth(13)=31;
        nbjyearn(1)=nbjpastmonth(1);
        for i=2:13
           nbjyearn(i)=nbjyearn(i-1)+nbjpastmonth(i);
        end
   
    %here for all bissextil years
        nbjpastmonth(1)=0;
        nbjpastmonth(2)=31;
        nbjpastmonth(3)=29;
        nbjpastmonth(4)=31;
        nbjpastmonth(5)=30;
        nbjpastmonth(6)=31;
        nbjpastmonth(7)=30;
        nbjpastmonth(8)=31;
        nbjpastmonth(9)=31;
        nbjpastmonth(10)=30;
        nbjpastmonth(11)=31;
        nbjpastmonth(12)=30;
        nbjpastmonth(13)=31;
        nbjyearb(1)=nbjpastmonth(1);
        for i=2:13
           nbjyearb(i)=nbjyearb(i-1)+nbjpastmonth(i);
        end
        
        %compiling
        nbjyear(:,1)=nbjyearn(:);  %1997
        nbjyear(:,2)=nbjyearn(:);  %1998  
        nbjyear(:,3)=nbjyearn(:);  %1999   
        nbjyear(:,4)=nbjyearb(:);  %2000
        nbjyear(:,5)=nbjyearn(:);  %2001
        nbjyear(:,6)=nbjyearn(:);  %2002
        nbjyear(:,7)=nbjyearn(:);  %2003
        nbjyear(:,8)=nbjyearb(:);  %2004
        nbjyear(:,9)=nbjyearn(:);  %2005
        nbjyear(:,10)=nbjyearn(:);  %2006
        nbjyear(:,11)=nbjyearn(:);  %2007
        %size(nbjyear)

        nbjsince97(1)=0;
        for year=2:12
           %sum(nbjyear(1:13,year),1)
           nbjsince97(year)=nbjsince97(year-1)+nbjyear(13,year-1);
        end  
        
  %end subroutine