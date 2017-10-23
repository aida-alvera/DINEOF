function outliers_field=dineof_outliers(dataini,datarec,mask,crossval,mu2_eff,vlsng,lftvec,N)

%function outliers_field=dineof_outliers(dataini,datarec,mask,crossval,mu2_eff,vlsng,lftvec,N)
%
% Calculates an outlier field based on three tests:
%                  -EOF test
%                  -Proximity test
%                  -Median test
%
% Routines needed to run the detection of outliers:
%
%    mu2.m to calculate mu2_eff
%    dineof_outliers_config.m to determine parameters of outlier detection
%    dineof_errorfield.m for the EOF test
%    proximity_test.m for the proximity test
%    mediant_test.m for the median test
% 
% Output: outliers_field: a 2D field with the size of the input data set
%                         with: a value of 1 at pixels classiffied as outliers
%                               a value of 0 everywhere else
% 



dineof_outliers_config;

if C1 ~=0;
  disp('Difference respect to EOF analysis');
  [error_field,datarecOI,indi_field]=dineof_errorfield(dataini,datarec,mask,crossval,mu2_eff,vlsng,lftvec,N);
else
  indi_field=0;
end

if C2~=0;
  disp('Cloud proximity test');
  indi_field_cl = proximity_test(dataini,mask);
  indi_field_cl = indi_field_cl*3;
else
  indi_field_cl=0;
end

if C3~=0;
  disp('Median test');
  indi_field_med = median_test(dataini,mask,mbox);
else
  indi_field_med=0;
end

indi_combi = C1 * indi_field + C2 * indi_field_cl + C3 * indi_field_med;

outliers_field=double(indi_combi >= level);
