function dataout = median_test(dataini,mask,box)

M = size(dataini,1);
N = size(dataini,2);
Z = size(dataini,3);

dataout=NaN*zeros(M,N,Z);
disp(['median test: calculating... Please be patient'])

%bottom border has a smaller median box
j = 1:box-1;
for i = 1:M-box
  temp = dataini(i:i+box,j,:);
  mid_size=size(temp,1)*size(temp,2);
  temp = reshape(temp,[mid_size Z]);
  med=nanmedian(temp,1);
  med_mat = repmat(med,[mid_size,1]);
  mad = 1.4826 * nanmedian(abs(temp-med_mat),1);
  dataout(i,j,:)=(abs(squeeze(dataini(i,j,:))'-med))./mad;
end

%right border has a smaller median box
i = M-box+1:M;
  for j = 1:N-box
    temp = dataini(i,j:j+box,:);
    mid_size=size(temp,1)*size(temp,2);
    temp = reshape(temp,[mid_size Z]);
    med=nanmedian(temp,1);
    med_mat = repmat(med,[mid_size,1]);
    mad = 1.4826 * nanmedian(abs(temp-med_mat),1);
    dataout(i,j,:)=(abs(squeeze(dataini(i,j,:))'-med))./mad;
  end


%top border has a smaller median box
j = N-box+1:N;
for i = box+1:M-box
  temp = dataini(i:i+box,j,:);
  mid_size=size(temp,1)*size(temp,2);
  temp = reshape(temp,[mid_size Z]);
  med=nanmedian(temp,1);
  med_mat = repmat(med,[mid_size,1]);
  mad = 1.4826 * nanmedian(abs(temp-med_mat),1);
  dataout(i,j,:)=(abs(squeeze(dataini(i,j,:))'-med))./mad;
end

%left border has a smaller median box
i = 1:box;
for j = box:N-box
  temp = dataini(i,j:j+box,:);
  mid_size=size(temp,1)*size(temp,2);
  temp = reshape(temp,[mid_size Z]);
  med=nanmedian(temp,1);
  med_mat = repmat(med,[mid_size,1]);
  mad = 1.4826 * nanmedian(abs(temp-med_mat),1);
  dataout(i,j,:)=(abs(squeeze(dataini(i,j,:))'-med))./mad;
end






%central domain
for i=box:M-box
  for j=box:N-box
     if mask(i,j)==1;
  %      for k=1:Z          
          temp = dataini(i-box+1:i+box,j-box+1:j+box,:);
          mid_size=size(temp,1)*size(temp,2);
          temp = reshape(temp,[mid_size Z]);
          med=nanmedian(temp,1);
          med_mat = repmat(med,[mid_size,1]);
          mad = 1.4826 * nanmedian(abs(temp-med_mat),1);
          dataout(i,j,:)=(abs(squeeze(dataini(i,j,:))'-med))./mad;
  %      end
     end
  end
end


          