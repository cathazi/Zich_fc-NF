%% transform ROIs and get offline time series
clear all; close all; clc;


SPMPATH ='F:\FCNF\spm12\';addpath(SPMPATH);

MAINPATH = 'F:\FCNF\';
PATHIN1 = [MAINPATH 'rawdata\'];
PATHIN2 = [MAINPATH 'data\'];
PATHOUT = [MAINPATH 'data\ROI_MNI2\'];mkdir(PATHOUT)
SCRIPTPATH = [MAINPATH 'scripts_FCNF\'];

DICOM_ID={'20170315.F3T_2015_32_053.F3T_2015_32_053'};
SUBJ_ID={'s111'};

X_max=96;
Y_max=96;
Z_max=72;
%% transform ROI's from TBV to MNI

%spm fmri
for s=1:length(SUBJ_ID)
    %%
    PATHIN_SUBJ=[PATHIN1 DICOM_ID{s} '\TBVFiles\target_localizer\'];
    cd(PATHIN_SUBJ);
    vals=[];
    temp11=[];
    temp11=dir('*.roi'); % subject space
    
    %% find N of voxels per ROI
    fid = fopen(temp11(2).name);
    tline = fgets(fid);
    
    counter = 0;
    while ischar(tline)
        cd(SCRIPTPATH);
        if(strfind(tline, 'NrOfVoxels'))
            counter = counter+1;
            roisplit = strsplit(tline,':');
            roi_N(counter) = str2double(roisplit{2});
        end%if
        tline = fgets(fid);
    end%while
    
    fclose(fid);
    
    %% get ROI coordinates
    roi=[];
    cd(PATHIN_SUBJ);
    fid = fopen(temp11(2).name);
    tline = fgets(fid);
    
    while ischar(tline)
        if sum(tline)==0
            tline=fgetl(fid);
        elseif strfind(tline, 'FileVersion')
            tline=fgetl(fid);
        elseif strfind(tline, 'SaveVoxelsInROIs')
            tline=fgetl(fid);
        elseif strfind(tline, 'SaveSortedVoxelList')
            tline=fgetl(fid);
        elseif strfind(tline, 'SaveROIMaxPSC')
            tline=fgetl(fid);
        elseif strfind(tline, 'NrOfTCPlots')
            tline=fgetl(fid);
        elseif strfind(tline, 'NrOfTCs')
            tline=fgetl(fid);
        elseif strfind(tline, 'FromSlice')
            tline=fgetl(fid);
        elseif strfind(tline, 'Left')
            tline=fgetl(fid);
        elseif strfind(tline, 'Right')
            tline=fgetl(fid);
        elseif strfind(tline, 'Top')
            tline=fgetl(fid);
        elseif strfind(tline, 'Bottom')
            tline=fgetl(fid);
        elseif strfind(tline, 'NrOfVoxels')
            tline=fgetl(fid);
        else
            roi(end+1,:)=str2num(tline);
            tline=fgetl(fid);
        end%if
    end%while
    
    %% bring dimensions from subject TBV DICOM space to subject SPM NIFTI space
    roi_nii=[];
    roi_nii(:,1)=roi(:,1);
    roi_nii(:,2)=Y_max-roi(:,2);
    roi_nii(:,3)=roi(:,3);
    
    %% bring scale from TBV voxel scale to SPM mm scale
    PATHIN_SUBJ=[PATHIN2 SUBJ_ID{s} '/NF/'];
    cd(PATHIN_SUBJ);
    
    temp=dir('*RT*');
    SESS=[];
    for n=1:length(temp)
        SESS{end+1}=[temp(n).name,'/'];
    end
    PATHIN_SESS = [PATHIN_SUBJ SESS{1}];
    
    Vol = [spm_select('List',PATHIN_SESS,'^meanf*','.nii')]; % as realignment is performed to mean image
    Vols = strcat(PATHIN_SESS,Vol,',1');
    
    cd(SCRIPTPATH);
    roi_nii_mm= map_coords(roi_nii', Vols);
    
    %% seperate ROIs
    ROI1=roi_nii_mm(:,1:roi_N(1));
    ROI2=roi_nii_mm(:,roi_N(1)+1:roi_N(1)+roi_N(2));
    ROI3=roi_nii_mm(:,roi_N(1)+roi_N(2)+1:size(roi,1));
    
    %% save ROIs as .txt
    cd(PATHOUT);
    dlmwrite([SUBJ_ID{s},'_roi_1.txt'],ROI1,'delimiter','\t','newline','pc');
    dlmwrite([SUBJ_ID{s},'_roi_2.txt'],ROI2,'delimiter','\t','newline','pc');
    dlmwrite([SUBJ_ID{s},'_roi_3.txt'],ROI3,'delimiter','\t','newline','pc');
    
    for r=1:3
        %% transform to marsbar objects
        vals=[];
        vals = spm_load([SUBJ_ID{s},'_roi_',num2str(r),'.txt']);
        
        % get subj spec. dir and sessions names
        PATHIN_SUBJ=[PATHIN2 SUBJ_ID{s} '/NF/'];
        cd(PATHIN_SUBJ);
        temp=dir('*RT*');
        SESS=[];
        for n=1:length(temp)
            SESS{end+1}=[temp(n).name,'/'];
        end
        PATHIN_SESS = [PATHIN_SUBJ SESS{1}];
        
        Vol = [spm_select('List',PATHIN_SESS,'^meanf*','.nii')];
        Vols = strcat(PATHIN_SESS,Vol,',1');
        affine = spm_get_space(Vols);
        
        params = struct('XYZ',vals,'mat',affine);
        roi_final=maroi_pointlist(params,'mm');
        cd(PATHOUT);
        saveroi(roi_final,[SUBJ_ID{s},'_roi_',num2str(r),'_SPM_roi.mat']); %abbrevation _roi.mat must be included
        
        %% transform .mat to .nii file
        fname=[SUBJ_ID{s},'_roi_',num2str(r),'_SPM.nii'];
        sp=maroi('classdata','space');
        save_as_image(roi_final,fname,sp);
    end %roi
    
    %% coregister ROIs to MNI space using coregistration matrix from preprocessing
    cd(SCRIPTPATH);
    load norm_roi.mat
    
    PATHIN_struct=cellstr([PATHIN_SUBJ 'structure/']);
    norm_params=spm_select('List', PATHIN_struct,'y_','.nii');
    matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = strcat(PATHIN_struct,norm_params);
    
    % images to write = roi .nii images
    cd(PATHOUT);
    roi_nii_files =[];
    roi_nii_files =spm_select('List',PATHOUT,SUBJ_ID{s},'.nii');
    
    for file = 1:size(roi_nii_files,1)
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample{file,1}=strcat(PATHOUT,roi_nii_files(file,:),',1');
    end
    
    % seem to need one rf image as well
    cd(PATHIN_SESS);
    Vol = [spm_select('List',PATHIN_SESS,'^rf*','.nii')];
    Vols = strcat(PATHIN_SESS,Vol);
    copyfile(Vols(10,:),PATHOUT,'f');
    
    cd(PATHOUT);
    rf_file =[];
    rf_file =spm_select('List',PATHOUT,'^rf*','.nii');
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample{4,1}=strcat(PATHOUT,rf_file,',1');
    
    % run normalization
    spm_jobman('run',matlabbatch)
    
    delete([PATHOUT rf_file]);
    delete([PATHOUT spm_select('List',PATHOUT,'^wrf*','.nii')]);
    
    %% extract offline time series from MNI swrf images thorugh MNI transformed SS roi
    for r=1:3
        clearvars roi
        mars_img2rois(['w',SUBJ_ID{s},'_roi_',num2str(r),'_SPM.nii'],PATHOUT,['w',SUBJ_ID{s},'_roi_',num2str(r),'_SPM'],'i')
    end
    
    roi_files=[];
    roi_files = spm_select('List',PATHOUT,SUBJ_ID(s),'.mat'); % selection might get difficult with more participants
    roi_files(1:3,:)=[];
    rois = maroi('load_cell', roi_files); % make maroi ROI objects
    
    % get swrf images
    for n=1:length(SESS)
        swrf_img_sess=[];
        PATHIN_SESS = [PATHIN_SUBJ SESS{n}];
        cd(PATHIN_SESS);
        swrf_img_sess=spm_select('List',PATHIN_SESS,'^swrf.*','.nii');
        
        mY = get_marsy(rois{:}, swrf_img_sess, 'mean'); % extract data into marsy data object
        Y  = summary_data(mY); % get summary time course(s)
        cd(PATHOUT)
        save([SUBJ_ID{s} '_time_SESS_' num2str(n) '.mat'],'Y');
    end % Sess
end % SUBJ

