function [tini]=buildtini(vsdini,nbjyear,nbjsince97)

 %building the vector (tini) of number of decimal days passed since the 01/01/1997 at the moment
 %of each input image listed in vsdini
        njinimax=size(vsdini,1);
        for j=1:njinimax
           yeararg=vsdini(j,3)-1996; 
           montharg=vsdini(j,4);
           nbjdec=vsdini(j,5)+vsdini(j,6)/24+vsdini(j,7)/(24*60)+nbjyear(montharg,yeararg)+nbjsince97(yeararg);
           tini(j)=nbjdec;
        end