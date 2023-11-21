% 25/06/16
% fmri preprocessing script
% before start: 
% 1. create FCNF folder, containing the following folders: data, DICOM, scripts_FCNF
% 2. copy the folder from the RT PC in DICOM folder
% 3. chane DICOM_ID and SUBJ_ID in line 16 and 17 accordingly
% 4. change line 11 to 14

close all;clear all;clc;

SPMPATH ='F:\FCNF\spm12\';addpath(SPMPATH);
MAINPATH ='F:\FCNF\';
RAWPATH ='F:\FCNF\rawdata\';
SCRIPTPATH ='F:\FCNF\scripts_FCNF\';

DICOM_ID={'20170315.F3T_2015_32_053.F3T_2015_32_053'};
SUBJ_ID={'s111'};

%%
for s=1:length(SUBJ_ID)
    spm fmri
    
        %% dicom
        cd(SCRIPTPATH);
        load dicom_import.mat
        PATHIN_SUBJ=[RAWPATH DICOM_ID{s} '\'];
        
        % data
        matlabbatch{1}.spm.util.import.dicom.data=[];
        rawdata = spm_select('List',PATHIN_SUBJ,'.');
        rawdata = strcat(PATHIN_SUBJ,rawdata);
        
        matlabbatch{1}.spm.util.import.dicom.data{1}=rawdata(1,:);
        for file=1:size(rawdata,1)
            matlabbatch{1}.spm.util.import.dicom.data(file,:)={[rawdata(file,:)]};
        end
        
        %path
        PATHOUT=[MAINPATH 'data\' SUBJ_ID{s} '\'];
        mkdir(PATHOUT);
        matlabbatch{1}.spm.util.import.dicom.outdir{1}=PATHOUT;
        
        %run batc
        spm_jobman('run',matlabbatch)
        
        %% sort folders
        PATHNF=[PATHOUT 'NF\'];
        mkdir(PATHNF);
        
        % copy RT folders and structure
        cd(PATHOUT);
        temp1=dir('*201*');
        copyfile([PATHOUT '\' temp1.name '\bold_mbep2d_2mm_MB6_v2_RT*'],PATHNF,'f');
        copyfile([PATHOUT '\' temp1.name '\t1_mpr_sag_1mm_iso_32ch_v2_0002\s*'],[PATHNF '\structure'],'f');
        
        % keep only full runs (310 volumes)
        cd([PATHNF]);
        temp2=dir('*RT*');
        SESS=[];
        for n=1:length(temp2)
            if length(dir([PATHNF,temp2(n).name,'\*nii']))==310
                SESS{end+1}=[temp2(n).name,'\'];
            else
                rmdir(temp2(n).name,'s');
            end
        end

        %% realignmet
        cd(SCRIPTPATH);
        load realign.mat
        matlabbatch{1}.spm.spatial.realign.estwrite.data=[];
        
        f_img_sess=[];
        
        for n=1:length(SESS)
            PATHIN_SESS= [PATHNF SESS{n}];
            %load functional data for all sessions
            f_img_sess =[spm_select('List',PATHIN_SESS,'^f.*','nii')];
            %get the path of the file before the name of file
            f_img_sess =strcat(PATHIN_SESS,f_img_sess,',1');
            matlabbatch{1}.spm.spatial.realign.estwrite.data{n}(1,:)={f_img_sess(1,:)};
            % write path+name in spm structure as cell
            for file = 1:size(f_img_sess,1)
                matlabbatch{1}.spm.spatial.realign.estwrite.data{n}(file,:) = {f_img_sess(file,:)};
            end
        end
        % run estimate & reslice
        spm_jobman('run',matlabbatch)
        
        %% Coreg
        cd(SCRIPTPATH);
        load coreg.mat
        
        %ref img
        PATHIN_mean=strcat(PATHNF,SESS(1));
        mean_funct_img =spm_select('List',PATHIN_mean,'^mean.*','nii');
        matlabbatch{1}.spm.spatial.coreg.estimate.ref =strcat(PATHIN_mean,mean_funct_img,',1');
        
        %source img
        PATHIN_struct=cellstr([PATHNF 'structure\']);
        struct_img=spm_select('List', PATHIN_struct,'^s.*','nii');
        matlabbatch{1}.spm.spatial.coreg.estimate.source = strcat(PATHIN_struct,struct_img,',1');
        
        %all rf's
        matlabbatch{1}.spm.spatial.coreg.estimate.other=[];
        count=[1 1+310 1+2*310 1+3*310];
        for n=1:length(SESS)
            PATHIN_SESS= cellstr([PATHNF SESS{n}]);
            %load rf's
            rf_img_sess_new =[];
            rf_img_sess_new =spm_select('List',PATHIN_SESS,'^rf.*','nii');
            %get the path of the file before the name of file
            rf_img_sess=[];
            rf_img_sess =[rf_img_sess;strcat(PATHIN_SESS,rf_img_sess_new,',1')];
            % write path+name in spm structure as cell
            for file = count(n):count(n)+size(rf_img_sess,1)-1 
                matlabbatch{1}.spm.spatial.coreg.estimate.other{file} = rf_img_sess{file-(count(n)-1)};
            end
        end
        matlabbatch{1}.spm.spatial.coreg.estimate.other=matlabbatch{1}.spm.spatial.coreg.estimate.other';
        %run coreg
        spm_jobman('run',matlabbatch)
        
        %% Segmentation
        cd(SCRIPTPATH);
        load segment.mat
        
        matlabbatch{1}.spm.spatial.preproc.channel.vols = strcat(PATHIN_struct,struct_img,',1');
        for p=1:6
            matlabbatch{1,1}.spm.spatial.preproc.tissue(1,p).tpm={strcat(SPMPATH, 'tpm\TPM.nii,', num2str(p))};
        end        
        %run segment
        spm_jobman('run',matlabbatch)
        
        %% Normalise 1
        cd(SCRIPTPATH);
        load norm_funct.mat
        
        % parameter file (normalisation parameter = output of segmentation)
        norm_params=spm_select('List', PATHIN_struct,'y_','.nii');
        matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = strcat(PATHIN_struct,norm_params);
        
        % images to write (mean & all rf's) same as in coerg 2 "other images"
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample=[];
        % mean
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1}=strcat(PATHIN_mean,mean_funct_img,',1');        % all rf's
        count=[2 2+310 2+2*310 2+3*310];
        
        for n=1:length(SESS)
            PATHIN_SESS= [PATHNF SESS{n}]
            %load rf's
            rf_img_sess_new =[];
            rf_img_sess_new =[spm_select('List',PATHIN_SESS,'^rf.*','nii')];
            %get the path of the file before the name of file
            rf_img_sess=[]; 
            rf_img_sess =[rf_img_sess;strcat(PATHIN_SESS,rf_img_sess_new,',1')];
            % write path+name in spm structure as cell
            for file = count(n):count(n)+size(rf_img_sess,1)-1
                 matlabbatch{1}.spm.spatial.normalise.write.subj.resample{file} = rf_img_sess(file-(count(n)-1),:);
            end %file
        end %Sess
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample= matlabbatch{1}.spm.spatial.normalise.write.subj.resample';
        
        %run norm1
        spm_jobman('run',matlabbatch)    
        
        %% Normalise 2
        cd(SCRIPTPATH);
        load norm_struct.mat
        %parameter file = same as for norm1
        
        matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = strcat(PATHIN_struct,norm_params);
        
        % images to write structural and mean corrected structural image
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample(1,:)=strcat(PATHIN_struct,struct_img,',1'); %structural one
        mean_corr_struct_img=spm_select('List', PATHIN_struct,'^ms.*','nii'); %mean corrected structural image
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample(2,:)=strcat(PATHIN_struct,mean_corr_struct_img,',1');
        
        % settings changed fromm default:
        %matlabbatch{6}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
        %matlabbatch{6}.spm.spatial.normalise.write.roptions.interp = 4;
        
        % run norm2
        spm_jobman('run',matlabbatch)
        
        %% Smooth
        cd(SCRIPTPATH);
        load smooth.mat
        matlabbatch{1}.spm.spatial.smooth.data=[];
        
        % all functional wrf's
        start_file=0; % because mean img. is the first one
        wrf_img_sess=[]; % at the beginning empty and get filled more in each session iteration
        
        for n=1:length(SESS)
            PATHIN_SESS= [PATHNF SESS{n}];
            %load rf's
            wrf_data_sess_new =[spm_select('List',PATHIN_SESS,'^wrf.*','nii')];
            %get the path of the file before the name of file
            wrf_img_sess =[wrf_img_sess;strcat(PATHIN_SESS,wrf_data_sess_new,',1')];
            % write path+name in spm structure as cell
            matlabbatch{1}.spm.spatial.smooth.data{1}=wrf_img_sess(1);
            for file = start_file+1:start_file+size(wrf_img_sess,1) % mean img., img.'s of sessions
                matlabbatch{1}.spm.spatial.smooth.data(file,:) = {[wrf_img_sess(file-start_file,:)]};
            end
        end
        
        % run smooth
        spm_jobman('run',matlabbatch)
                
end % Subj

%% Validation
for s=1:length(SUBJ_ID)
    PATHOUT=[MAINPATH 'data\' SUBJ_ID{s} '\'];
    PATHNF=[PATHOUT 'NF\'];
    PATHIN_struct=cellstr([PATHNF 'structure\']);
    
    cd([PATHNF]);
    temp2=dir('*RT*');
    SESS=[];
    for n=1:length(temp2)
        if length(dir([PATHNF,temp2(n).name,'\*nii']))>310
            SESS{end+1}=[temp2(n).name,'\'];
        end
    end
    
    %% registration
    
    cd(SCRIPTPATH);
    load check_reg
    % T1 always the same
    % template
    matlabbatch{1,1}.spm.util.checkreg.data(1)={strcat(SPMPATH,'toolbox\DARTEL\icbm152.nii')};
    % structure aligned
    struct_img=spm_select('List', PATHIN_struct,'^wms.*','nii');
    matlabbatch{1,1}.spm.util.checkreg.data(2)=[strcat(PATHIN_struct,struct_img,',1')];
    % function aligned
    PATHIN_SESS= [PATHNF SESS{1}];
    swrf_img=spm_select('List', PATHIN_SESS,'^swrf.*','nii');
    matlabbatch{1,1}.spm.util.checkreg.data(3)={strcat(PATHIN_SESS,swrf_img(200,:),',1')};
    % structure original
    struct_img=spm_select('List', PATHIN_struct,'^s.*','nii');
    matlabbatch{1,1}.spm.util.checkreg.data(4)=[strcat(PATHIN_struct,struct_img,',1')];
    % unciton original
    PATHIN_SESS= [PATHNF SESS{1}];
    f_img=spm_select('List', PATHIN_SESS,'^f.*','nii');
    matlabbatch{1,1}.spm.util.checkreg.data(5)={strcat(PATHIN_SESS,f_img(200,:),',1')};
    
    % run batch
    spm_jobman('run',matlabbatch)
    print('-dpng', [PATHIN_SESS, SUBJ_ID{s},'_check_reg']);
    
    
    %% motion
    cd(SCRIPTPATH);
    for n=1:length(SESS)
        PATHIN_SESS= [PATHNF SESS{n}];
        copyfile([SCRIPTPATH 'get_motion_params.m'],PATHIN_SESS,'f');
        cd(PATHIN_SESS);
        get_motion_params
        print('-dpng', [PATHIN_SESS, SUBJ_ID{s},'_motion_params']);
    end
end % SUBJ

