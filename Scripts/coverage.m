function cov=coverage(data,mask,ty)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% function cov = coverage(matrix,mask,type)                    %
% Calculates coverage in space or time of a given 3D matrix    %
%                                                              %
% Input: matrix: 3D matrix with NaN in land and clouds         %
%        mask with 1 on sea, 0 on land                         %
%        ty: String, type of average                           %
%                 ('sp' for a spatial distribution of coverage %
%               or 'tm' for temporal distribution of coverage) %
% Output: % of coverage                                        %
%                                                              %
% Aida Alvera-Azcarate, January 2005                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

size1=size(data,1);
size2=size(data,2);
size3=size(data,3);
 
if ty=='sp'
    count=zeros(size1,size2);
    for i=1:size1;
        for j=1:size2;
            for k=1:size3;
                if mask(i,j)==1 & isnan(data(i,j,k));
                    count(i,j)=count(i,j)+1;
                end
            end
        end
    end
    
    cov=count*100/size3;
    
elseif ty=='tm'
    count=zeros(size3,1);
    for i=1:size1;
        for j=1:size2;
            for k=1:size3;
                if mask(i,j)==1 & isnan(data(i,j,k));
                    count(k)=count(k)+1;
                end
            end
        end
    end
    
    points=sum(mask(:));
    
    cov=count*100/points;
   
    
    
else 'incorrect type of average, please indicate spatial average, sp, or temporal average, tm'
    
end 
