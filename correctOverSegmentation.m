function bw_nuc=correctOverSegmentation(bw_nuc,mean_nuc_size)
% this function identifies small objects as a results of over-segmentation
% and tries to merge them back to their parent objects

min_nuc_area=round(mean_nuc_size^2/4);
% try to correct for over-segmentation by merging small objects back
% together
l_nuc=bwlabel(bw_nuc);
nuc_cnt=max(l_nuc(:));
props=regionprops(l_nuc,'Area','BoundingBox');
area=cat(1,props.Area);
box=cat(1,props.BoundingBox);
% cut patches to work with individual nucleus
for i=1:nuc_cnt
    % check if the area is below the cut-off value
    if (area(i) < min_nuc_area)
        % cut patch
        ymin=max(floor(box(i,2))-3,1);
        ymax=min(ceil(ymin+box(i,4)+1)+3,ylim);
        xmin=max(floor(box(i,1))-3,1);
        xmax=min(ceil(xmin+box(i,3)+1)+3,xlim);
        patch1=l_nuc(ymin:ymax,xmin:xmax);
        % check if the small object overlaps with its neighbors within
        % given pixel radius
        tmp=imdilate(patch1==i,strel('disk',4));
        tmp=unique(patch1(tmp));
        tmp(tmp==0)=[];
        if (length(tmp) > 2) % more than one neighbors
            tmp(tmp==i)=[];
            index=tmp(1);
            for j=2:length(tmp)
                % try to join the current object with its smaller neighbor
                if (area(tmp(j)) < area(index))
                    index=tmp(j);
                end
            end
        elseif (length(tmp) == 2) % only one neighbor
            tmp(tmp==i)=[];
            index=tmp;
        else % no neighbor
            index=-1;
        end
        
        % try to join the small object to its neighbor of choice
        if (index > 0)
            tmp=imdilate(patch1==i,strel('disk',2))+imdilate(patch1==index,strel('disk',2));
            tmp=tmp>1;
            tmp=tmp*index;
            patch1=max(patch1,tmp);
            patch1(patch1==i)=index;
            l_nuc(ymin:ymax,xmin:xmax)=max(l_nuc(ymin:ymax,xmin:xmax),patch1);
        end
    end
end
bw_nuc=l_nuc>0;
bw_nuc=bwareaopen(bw_nuc,min_nuc_area);

end
