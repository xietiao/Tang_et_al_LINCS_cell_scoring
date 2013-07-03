%score_cells_LTR_Nucview_assay(data_path) score cell phenotypes from Hoechst-LTR-NucView assay
%
% Inputs     data_path        : path where the data is stored
%
% Ouputs
%            
% Tiao Xie, May 2011 (last modified 07/02/2013)

function score_cells_LTR_Nucview_assay(data_path)

%% input section
% define I/O paths
data_path=[data_path filesep];
input_path=[data_path 'input_images' filesep];
output1_path=[data_path 'output1_images' filesep];
output2_path=[data_path 'output2_images' filesep];
output3_path=[data_path 'output3_images' filesep];

% some parameters for analysis
code_version='score_cells_LTR_Nucview_assay';
nsite=4;
max_nuc_size=50;
mean_nuc_size=20;
intial_nuc_int_cutoff=0.75;
min_nuc_int_cutoff=0.1;
nuc_int_cutoff_scaler=0.15;
max_nuc_int_cutoff=0.4;
LTR_th_scaler=0.9;
min_mito_ffactor=0.75;
min_mito_mal=20;
max_mito_mal=50;
min_mito_area=250;
LTR_blk_size=100;
mito_LTR_int_above_bg=300;
min_mito_LTR_coverage=0.75;
fitc_smooth_filter_size=2;
fitc_blk_size=50;
fitc_above_bg=25;
fitc_join_radius=5;
min_fitc_area=9;
nuc_dilate_radius=50;
min_LTR_int_above_bg=50;
min_LTR_area=50;
LTR_join_radius=5;
max_dead_nuc_LTR_coverage=0.5;
min_dead_nuc_area=1000;
% print parameters to file
fid1=fopen([data_path 'parameters.txt'],'wt');
fprintf(fid1,'code_version = %s\n',code_version);
fprintf(fid1,'nsite = %s\n',num2str(nsite));
fprintf(fid1,'max_nuc_size = %s\n',num2str(max_nuc_size));
fprintf(fid1,'mean_nuc_size = %s\n',num2str(mean_nuc_size));
fprintf(fid1,'intial_nuc_int_cutoff = %s\n',num2str(intial_nuc_int_cutoff));
fprintf(fid1,'min_nuc_int_cutoff = %s\n',num2str(min_nuc_int_cutoff));
fprintf(fid1,'nuc_int_cutoff_scaler = %s\n',num2str(nuc_int_cutoff_scaler));
fprintf(fid1,'max_nuc_int_cutoff = %s\n',num2str(max_nuc_int_cutoff));
fprintf(fid1,'LTR_th_scaler = %s\n',num2str(LTR_th_scaler));
fprintf(fid1,'min_mito_ffactor = %s\n',num2str(min_mito_ffactor));
fprintf(fid1,'min_mito_mal = %s\n',num2str(min_mito_mal));
fprintf(fid1,'max_mito_mal = %s\n',num2str(max_mito_mal));
fprintf(fid1,'min_mito_area = %s\n',num2str(min_mito_area));
fprintf(fid1,'LTR_blk_size = %s\n',num2str(LTR_blk_size));
fprintf(fid1,'mito_LTR_int_above_bg = %s\n',num2str(mito_LTR_int_above_bg));
fprintf(fid1,'min_LTR_percent_bright_area = %s\n',num2str(min_mito_LTR_coverage));
fprintf(fid1,'fitc_smooth_filter_size = %s\n',num2str(fitc_smooth_filter_size));
fprintf(fid1,'fitc_blk_size = %s\n',num2str(fitc_blk_size));
fprintf(fid1,'fitc_above_bg = %s\n',num2str(fitc_above_bg));
fprintf(fid1,'fitc_join_radius = %s\n',num2str(fitc_join_radius));
fprintf(fid1,'min_fitc_area = %s\n',num2str(min_fitc_area));
fprintf(fid1,'nuc_dilate_radius = %s\n',num2str(nuc_dilate_radius));
fprintf(fid1,'min_LTR_int_above_bg = %s\n',num2str(min_LTR_int_above_bg));
fprintf(fid1,'min_LTR_area = %s\n',num2str(min_LTR_area));
fprintf(fid1,'LTR_join_radius = %s\n',num2str(LTR_join_radius));
fprintf(fid1,'max_dead_nuc_LTR_coverage = %s\n',num2str(max_dead_nuc_LTR_coverage));
fprintf(fid1,'min_dead_nuc_area = %s\n',num2str(min_dead_nuc_area));
fclose(fid1);

%% load in filenames of input images
% text files that contains all the image file names
ch1_list=fopen([input_path 'blue_files.txt'],'r');
ch2_list=fopen([input_path 'green_files.txt'],'r');
ch3_list=fopen([input_path 'red_files.txt'],'r');
% read in all the data file names
i=0;
tline1=fgetl(ch1_list);
tline2=fgetl(ch2_list);
tline3=fgetl(ch3_list);
while ischar(tline1)
    i=i+1;
    ch1_filenames{i,1}=tline1;
    ch2_filenames{i,1}=tline2;
    ch3_filenames{i,1}=tline3;
    tline1=fgetl(ch1_list);
    tline2=fgetl(ch2_list);
    tline3=fgetl(ch3_list);
end
fclose(ch1_list);
fclose(ch2_list);
fclose(ch3_list);
% find out the number of images from the number of data file names
nimages=length(ch1_filenames);
nwell=nimages/nsite;

%% initialize output files
% write per site data to file
fid2=fopen([data_path 'Summary_site_data.csv'],'wt');
fprintf(fid2,'Well,site,Cell_cnt,Interphase_cnt,Apoptotic_cnt,Dead_cnt,Mitotic_cnt,');
fprintf(fid2,'Percent_interphase,Percent_apoptotic,Percent_mitotic,');
fprintf(fid2,'Median_interphase_area,Median_interphase_formfactor\n');

% write per well data to file
fid3=fopen([data_path 'Summary_well_data.csv'],'wt');
fprintf(fid3,'Well,Cell_cnt,Interphase_cnt,Apoptotic_cnt,Dead_cnt,Mitotic_cnt,');
fprintf(fid3,'Percent_interphase,Percent_apoptotic,Percent_mitotic,');
fprintf(fid3,'Median_interphase_area,Median_interphase_formfactor\n');

%% start analysis
% loop over all the images to perform segmentation and collect intensity
% information from control wells
flag_1st_loop=1;
im_out1=cell(4,1);
im_out2=cell(4,1);
im_out3=cell(4,1);
for w=1:nwell
    % initiate all the per well data array
    wn_nuc=0;
    wn_interphase=0;
    wn_apoptotic=0;
    wn_dead=0;
    wn_mitotic=0;
    w_area=[];
    w_formfactor=[];
    for s=1:nsite
        % calculate the image #
        n=(w-1)*nsite+s;
        
        % pick up well and site information from the filename
        filename=ch1_filenames{n};
        ntmp=strfind(filename,'_');
        ntmp_length=length(ntmp);
        well=filename(ntmp(ntmp_length-2)+1:ntmp(ntmp_length-1)-1);
        site=filename(ntmp(ntmp_length-1)+1:ntmp(ntmp_length)-1);
        
        % read in RGB images into three arrays
        im1=imread([input_path ch1_filenames{n}]);
        im2=imread([input_path ch2_filenames{n}]);
        im3=imread([input_path ch3_filenames{n}]);
        
        % convert images to double
        im1=im2double(im1);
        im2=im2double(im2);
        im3=im2double(im3);
        
        % define montage arrays
        if (flag_1st_loop == 1)
            [ylim,xlim]=size(im1);
            % preparation for montage
            ny=ceil(sqrt(nsite));
            nx=ceil(nsite/ny);
            xi=zeros(nx*ny,1);
            xf=zeros(nx*ny,1);
            yi=zeros(nx*ny,1);
            yf=zeros(nx*ny,1);
            stmp=0;
            for y=1:ny
                for x=1:nx
                    stmp=stmp+1;
                    yi(stmp)=(y-1)*ylim+1;
                    yf(stmp)=y*ylim;
                    xi(stmp)=(x-1)*xlim+1;
                    xf(stmp)=x*xlim;
                end
            end
            
            % initialize the montage image arrays and count array
            montage=zeros([ny*ylim nx*xlim 3]);
        end
        
        %% Hoechst channel
        % initial segmentation using an adjustable threshold approach
        bw_nuc=segNucAdjustableThreshold(im1,max_nuc_size,intial_nuc_int_cutoff,nuc_int_cutoff_scaler,min_nuc_int_cutoff,max_nuc_int_cutoff);
        bw_nuc=sepClusteredNuc(bw_nuc,mean_nuc_size);
        bw_nuc=correctOverSegmentation(bw_nuc,mean_nuc_size);
        
        %% LTR channel
        % rough segment of the LTR channel based on intensity watershed
        % combine DAPI and LTR, and flip the sign
        tmp=-(im1+im3);
        % impose all the nuclei as global minima
        tmp=imimposemin(tmp,bw_nuc);
        % watershed segmentation
        l_cell=watershed(tmp);
        bw_cell=l_cell>0;
        % nuclear count
        nuc_cnt=max(l_cell(:));
        % update nuclear index
        l_nuc=double(bw_nuc).*double(l_cell);
        
        % threshold LTR image
        cell_size=max_nuc_size;
        if mod(cell_size,2) == 0
            cell_size=cell_size+1;
        end
        [im3_bg,~]=localAvgStd2D(im3,cell_size);
        bw_cell2=im3>im3_bg+mito_LTR_int_above_bg/65535;
        bw_cell2=imopen(bw_cell2,strel('disk',2));
        bw_cell2=imfill(bw_cell2,'holes');
        bw_cell2(~bw_cell)=0;
        l_cell2=bwlabel(bw_cell2);
        tmp=l_cell2;
        
        % measure the form factor of each object and filter out mitotic cells
        % based on this measurement
        props=regionprops(l_cell2,'Area','Perimeter','MajorAxisLength');
        area=cat(1,props.Area);
        perim=cat(1,props.Perimeter);
        formfactor=4*pi*area./max(perim.^2,1);
        mal=cat(1,props.MajorAxisLength);
        tmp2=zeros(size(im1));
        for i=1:max(tmp(:))
            % filter out spherical LTR blobs based on shape parameters
            if (formfactor(i) >= min_mito_ffactor && mal(i) > min_mito_mal && mal(i) < max_mito_mal && area(i) > min_mito_area)
                tmp2(tmp==i)=1;
            end
        end
        mito_cell=tmp2>0;
        l_mito_cell=bwlabel(mito_cell);
        
        % score cells based on its LTR coverage within nuclei
        mito_cell=zeros(size(im1));
        mito_nuc=zeros(size(im1));
        props=regionprops(l_nuc,'Area');
        area=cat(1,props.Area);
        for i=1:max(l_mito_cell(:))
            tmp=double(l_mito_cell==i).*double(l_nuc);
            tmp2=unique(tmp);
            tmp2(tmp2==0)=[];
            if (~isempty(tmp2))
                for j=1:length(tmp2)
                    ratio=sum(sum(tmp==tmp2(j)))/area(tmp2(j));
                    if (ratio > min_mito_LTR_coverage)
                        mito_cell(l_mito_cell==i)=1;
                        mito_nuc(l_nuc==tmp2(j))=1;
                    end
                end
            end
        end
        other_nuc=bw_nuc-mito_nuc;
        mito_cell=mito_cell>0;
        mito_nuc=mito_nuc>0;
        other_nuc=other_nuc>0;
        
        %% NucView channel
        % detect NucView spots
        im2_smooth=filterGauss2D(im2,1);
        im2_bg=filterGauss2D(im2,cell_size);
        bw_fitc=im2_smooth>(im2_bg+fitc_above_bg/65535);
        bw_fitc=imclose(bw_fitc,strel('disk',fitc_join_radius));
        bw_fitc=bwareaopen(bw_fitc,min_fitc_area);
        
        % classify apoptotic cells
        l_other_nuc=bwlabel(other_nuc);
        apop_nuc=zeros(size(im1));
        inter_nuc=zeros(size(im1));
        for i=1:max(l_other_nuc(:))
            if (max(max(bw_fitc(l_other_nuc==i))) > 0)
                apop_nuc(l_other_nuc==i)=1;
            else
                inter_nuc(l_other_nuc==i)=1;
            end
        end
        apop_nuc=apop_nuc>0;
        inter_nuc=inter_nuc>0;
        
        %% detect late stage dead cell population
        % remove small fragments in the interphase population
        min_nuc_area=round(mean_nuc_size^2/4);
        inter_nuc=bwareaopen(inter_nuc,min_nuc_area);
        
        % score late stage dead cells in the interphase cell population by
        % looking for the absence of LTR signal
        dead_nuc=zeros(size(im1));
        
        % segment individual LTR block carved out based on previous watershed
        tmp2=ones(size(im1));
        for i=1:nuc_cnt
            tmp=im3(l_cell==i);
            % determine the threshold based on intensity profile within the
            % block
            th_cell=quantile(tmp,0.01);
            % slight lower the threshold by 10%
            tmp2(l_cell==i)=th_cell;
        end
        % a regionally segmented LTR image is achievement
        tmp=im3>tmp2+min_LTR_int_above_bg/65535;
        tmp=imopen(tmp,strel('disk',2));
        
        % clean up each LTR block and fill holes 
        props=regionprops(l_cell,'BoundingBox');
        box=cat(1,props.BoundingBox);
        for i=1:max(l_cell(:))
            ymin=max(floor(box(i,2))-1,1);
            ymax=min(ceil(ymin+box(i,4)+1)+1,ylim);
            xmin=max(floor(box(i,1))-1,1);
            xmax=min(ceil(xmin+box(i,3)+1)+1,xlim);
            patch1=tmp(ymin:ymax,xmin:xmax);
            patch2=l_cell(ymin:ymax,xmin:xmax);
            patch1=imclose(patch1,strel('disk',5));
            patch1=imfill(patch1,'holes');
            patch1(patch2~=i)=false;
            tmp(ymin:ymax,xmin:xmax)=max(tmp(ymin:ymax,xmin:xmax),patch1);
        end
        
        % score cells based on its LTR coverage within nuclei
        l_inter_nuc=bwlabel(inter_nuc);
        props=regionprops(l_inter_nuc,'Area');
        area=cat(1,props.Area);
        tmp=double(tmp).*double(l_inter_nuc);
        for i=1:max(l_inter_nuc(:))
            ratio=sum(sum(tmp==i))/area(i);
            if (ratio < max_dead_nuc_LTR_coverage  && area(i) > min_dead_nuc_area)
                dead_nuc(l_inter_nuc==i)=1;
                inter_nuc(l_inter_nuc==i)=0;
            end
        end
        dead_nuc=dead_nuc>0;
        
        %% measurements of individual cell for output
        % measure area and formfactor of all the interphase cells
        l_inter_nuc=bwlabel(inter_nuc);
        props=regionprops(l_inter_nuc,'Area','Perimeter');
        area=cat(1,props.Area);
        perim=cat(1,props.Perimeter);
        formfactor=4*pi*area./max(perim.^2,1);
        
        % count cells
        n_interphase=max(max(bwlabel(inter_nuc)));
        n_mitotic=max(max(bwlabel(mito_nuc)));
        n_apoptotic=max(max(bwlabel(apop_nuc)));
        n_dead=max(max(bwlabel(dead_nuc)));
        n_nuc=n_interphase+n_mitotic+n_apoptotic+n_dead;
        % calculate percentages
        percent_interphase=n_interphase/max(n_nuc,1);
        percent_mitotic=n_mitotic/max(n_nuc,1);
        percent_dead=(n_apoptotic+n_dead)/max(n_nuc,1);
        
        % write per site data to file
        fprintf(fid2,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            well,site,num2str(n_nuc),num2str(n_interphase),num2str(n_apoptotic),...
            num2str(n_dead),num2str(n_mitotic),num2str(percent_interphase),...
            num2str(percent_dead),num2str(percent_mitotic),...
            num2str(median(area)),num2str(median(formfactor)));
        
        % accumulate over the whole well
        wn_nuc=wn_nuc+n_nuc;
        wn_interphase=wn_interphase+n_interphase;
        wn_apoptotic=wn_apoptotic+n_apoptotic;
        wn_dead=wn_dead+n_dead;
        wn_mitotic=wn_mitotic+n_mitotic;
        w_area=cat(1,w_area,area);
        w_formfactor=cat(1,w_formfactor,formfactor);
        
        %% create output images
        % create output images
        im1_adj=imadjust(im1,[quantile(im1(:),0.002),quantile(im1(:),0.998)],[0 1]);
        im2_adj=imadjust(im2,[quantile(im2(:),0.002),quantile(im2(:),0.998)],[0 1]);
        im3_adj=imadjust(im3,[quantile(im3(:),0.002),quantile(im3(:),0.998)],[0 1]);
        % outlintes
        p_inter=bwperim(inter_nuc);
        p_mito=bwperim(mito_nuc);
        p_apop=bwperim(apop_nuc);
        p_dead=bwperim(dead_nuc);
        % output images
        im_out1{s}=cat(3,max(max(max(im1_adj,p_inter),p_mito),p_dead)-p_apop,max(max(max(im1_adj,p_inter),p_apop),p_dead)-p_mito,max(im1_adj,p_inter)-p_mito-p_apop-p_dead);
        im_out2{s}=cat(3,im2_adj-bwperim(bw_fitc),max(im2_adj,bwperim(bw_fitc)),im2_adj-bwperim(bw_fitc));
        im_out3{s}=cat(3,max(im3_adj,bwperim(mito_cell)),im3_adj-bwperim(mito_cell),im3_adj-bwperim(mito_cell));
        
        % flag the 1st loop is over
        flag_1st_loop=0;
    end
    
    %% 
    % calculate percentages
    percent_interphase=wn_interphase/max(wn_nuc,1);
    percent_dead=(wn_apoptotic+wn_dead)/max(wn_nuc,1);
    percent_mitotic=wn_mitotic/max(wn_nuc,1);
    
    % write per site data to file
    fprintf(fid3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
        well,num2str(wn_nuc),num2str(wn_interphase),num2str(wn_apoptotic),...
        num2str(wn_dead),num2str(wn_mitotic),num2str(percent_interphase),...
        num2str(percent_dead),num2str(percent_mitotic),...
        num2str(median(w_area)),num2str(median(w_formfactor)));
    
    % montage output images
    for s=1:nsite
        montage(yi(s):yf(s),xi(s):xf(s),:)=im_out1{s};
    end
    imwrite(montage,[output1_path well '.jpg'],'jpg');
    for s=1:nsite
        montage(yi(s):yf(s),xi(s):xf(s),:)=im_out2{s};
    end
    imwrite(montage,[output2_path well '.jpg'],'jpg');
    for s=1:nsite
        montage(yi(s):yf(s),xi(s):xf(s),:)=im_out3{s};
    end
    imwrite(montage,[output3_path well '.jpg'],'jpg');
end

fclose(fid2);
fclose(fid3);

end

