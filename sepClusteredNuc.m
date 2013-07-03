function bw_nuc=sepClusteredNuc(bw_nuc,mean_nuc_size)
% this function takes a binary mask and tries to separate clustered nuclei

% seperate clustered nuclei based on shape
% distance transformation
tmp=bwdist(~bw_nuc);
% set minimum distance between two peaks
minDist=[mean_nuc_size mean_nuc_size];
% call function to find peak positions in the distance transformation
[xlmax,ylmax]=localMaximum(tmp,minDist,true);
% get rid of all the zeros in the background
i=1;
while (i <= max(size(xlmax)))
    if (tmp(xlmax(i),ylmax(i)) == 0)
        xlmax(i)=[];
        ylmax(i)=[];
    else
        i=i+1;
    end
end
% watershed the whole image
NUM=max(size(xlmax));
tmp=-tmp;
nuc_core=zeros(size(tmp));
for i=1:NUM
    nuc_core(xlmax(i),ylmax(i))=1;
end
tmp=imimposemin(tmp,nuc_core>0);
l_nuc=watershed(tmp);
l_nuc=bw_nuc.*double(l_nuc);
% binary image after seperating clustered objects
bw_nuc=l_nuc>0;

end
