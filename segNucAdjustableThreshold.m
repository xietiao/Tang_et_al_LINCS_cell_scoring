function bw_nuc=segNucAdjustableThreshold(im1,max_nuc_size,intial_nuc_int_cutoff,nuc_int_cutoff_scaler,min_nuc_int_cutoff,max_nuc_int_cutoff)
% A rough segmentation is performed based on local background estimation,
% followed by an adjusted thresholding determined from the intensity
% profile of individual nuclear patch. This approach gives a more accurate
% nuclear segmentation when dimmer objects are sitting next to really
% bright objects.

% initial segmentation of the DAPI channel
tmp=imopen(im1,strel('disk',max_nuc_size));
nuc_bg_scale=max_nuc_size*1.5;
if mod(nuc_bg_scale,2) == 0
    nuc_bg_scale=nuc_bg_scale+1;
end
[im1_bg,~]=localAvgStd2D(tmp,nuc_bg_scale);
bw_nuc=im1>(im1_bg*1.5);
l_nuc=bwlabel(bw_nuc);
nuc_cnt=max(l_nuc(:));

% refine nuclear segmentation by adjusting the threshold of individual
% nucleus based on its intensity profile
[ylim,xlim]=size(im1);
bw_nuc=zeros(size(im1));
% measure the bounding box and centroid for each nuclei
props = regionprops(l_nuc,im1,'BoundingBox');
box=cat(1,props.BoundingBox);
% calculate the intensity profile for each nucleus
int_profile=zeros(nuc_cnt,1);
for i=1:nuc_cnt
    int_profile(i)=quantile(im1(l_nuc==i),intial_nuc_int_cutoff);
end
min_int_profile=min(int_profile);
% cut patches to work with individual nucleus
for i=1:nuc_cnt
    % cut patch
    ymin=max(floor(box(i,2)),1);
    ymax=min(ceil(ymin+box(i,4)+1),ylim);
    xmin=max(floor(box(i,1)),1);
    xmax=min(ceil(xmin+box(i,3)+1),xlim);
    patch1=l_nuc(ymin:ymax,xmin:xmax);
    patch1=patch1==i;
    patch2=im1(ymin:ymax,xmin:xmax);
    patch2(~patch1)=0;
    % retrieve all the pixel values
    pv=patch2(patch1>0);
    % determine the cutoff percentile
    percent_cutoff=(int_profile(i)-min_int_profile)/max(min_int_profile,1/65535)*nuc_int_cutoff_scaler+min_nuc_int_cutoff;
    percent_cutoff=min(percent_cutoff,max_nuc_int_cutoff);
    % discard the lower 10% of all the pixels
    th_nuc=quantile(pv,percent_cutoff);
    % re-threshold the nucleus
    tmp=patch2>th_nuc;
    % update the segmentation of the whole image
    bw_nuc(ymin:ymax,xmin:xmax)=max(bw_nuc(ymin:ymax,xmin:xmax),tmp);
end

end
