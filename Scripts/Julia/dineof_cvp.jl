"""

function dineof_cvp(fname,maskfname,outdir,nbclean)
 Extracts clouds from file 'fname' and adds them to the 'nbclean'
   cleanest images from the data set

 Input variables:
       fname (string): file name in the disk (gher format)
       maskfname (string): mask file name in teh disk (gher format)
       outdir (string): directory where new file will be written
       nbclean (integer): number of cleanest images to be covered with clouds
  
 File with added clouds will be written in a file with the same name
    as the initial file (and at the same location in the disk)
    but with the extension '.clouds' added. An additional file, 
    clouds.index is also written, which contains the indexes of 
    the added clouds (can be used as cross-validation points in 
    DINEOF)

 Reference:
 A. Alvera-Azcarate, A. Barth, M. Rixen, and J. M. Beckers. 
    Reconstruction of incomplete oceanographic data sets using 
    Empirical Orthogonal Functions. Application to the Adriatic 
    Sea Surface Temperature. Ocean Modelling, 9:325-346, 2005.

"""
function dineof_cvp(fname,maskfname,outdir,nbclean)

file,varname = split(fname,"#");
@show file 
ds = Dataset(String(file));
tmp = ds[String(varname)][:];
SST  = nomissing(tmp,NaN);
close(ds)

mfile,mvarname=split(maskfname,"#");
ds = Dataset(String(mfile));
mask = ds[String(mvarname)][:];
close(ds)

#SST
#mask = gread(maskfname);

for k=1:size(SST,3)
  tmp = SST[:,:,k];
  tmp[mask .== 0] .= NaN;
  SST[:,:,k] = tmp;
end

nbland = sum( mask[:] .== 0 );
mmax = sum( mask[:] .== 1 );

cloudcov = (sum(sum(isnan.(SST),dims=2),dims=1) .- nbland)/mmax;
cloudcov = cloudcov[:];


clean = sortperm(cloudcov);
clean = clean[1:nbclean];

N = length(cloudcov);

index = floor.(Int,N * rand(nbclean,1)).+1;

while (any(index .== clean))
  index = floor.(Int,N * rand(nbclean,1)).+1;
end

#to be checked /~isnan
SST2 = copy(SST);
SST2[:,:,clean] =  SST[:,:,clean] ./ .!isnan.(SST[:,:,index]);
SST2[isinf.(SST2)] .= NaN;

imax = size(SST,1);
jmax = size(SST,2);


mindex = zeros(Int,imax,jmax);
iindex = zeros(Int,mmax);
jindex = zeros(Int,mmax);

m=0;
for i=1:imax    
  for j=1:jmax
    if (mask[i,j] == 1)
      m = m+1;
      mindex[i,j] = m;
      iindex[m] = i;
      jindex[m] = j;
    end 
  end 
end 



indexex = findall(isnan.(SST2) .& .!isnan.(SST));
#iex,jex,kex = ind2sub(size(SST2),indexex);

nbpoints = length(indexex)
clouds_indexes = zeros(nbpoints,2);

for l=1:nbpoints
  iex = CartesianIndices(size(SST2))[indexex[l]] 
  clouds_indexes[l,1] = mindex[iex[1],iex[2]];
  clouds_indexes[l,2] = iex[3];
end  

cloudcov2 = (sum(sum(isnan.(SST2),dims=2),dims=1) .- nbland)/m;
cloudcov2 = cloudcov2[:];

nbgood = sum(.!isnan.(SST[:]));
nbgood2 = sum(.!isnan.(SST2[:]));

println("$(100*(nbgood-nbgood2)/nbgood) % of cloud cover added")


output = Dataset(joinpath(outdir,"clouds_index.nc"),"c");
defDim(output,"nbpoints",size(clouds_indexes,1))
defDim(output,"index",size(clouds_indexes,2))

ncCloud = defVar(output,"clouds_index",Int64,("nbpoints","index"));
ncCloud[:] = clouds_indexes;

close(output)

end
