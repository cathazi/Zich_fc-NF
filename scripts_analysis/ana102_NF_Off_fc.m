% calc offline par corr values

clear all; close all; clc; 

MAINPATH='F:\FCNF\';
PATHIN_OFF=[MAINPATH 'data\ROI_MNI2\'];
PATHIN_ON=[MAINPATH 'Log\'];
PATHOUT=[MAINPATH 'data\ana02\'];
mkdir(PATHOUT);

SUBJ_ID={'s111'};
RUN={'1','2','3','4'};
WINDSIZE=19;

%% Off Vals
for s=1:length(SUBJ_ID)
    for r=1:length(RUN)
        cd(PATHIN_OFF);
        load ([SUBJ_ID{s} '_time_SESS_' num2str(r)]);
        %load ([SUBJ_ID{s} '_time_SESS_' num2str(r)]);
        for q=1:size(Y,2)%roi
            for k=1:size(Y,1)%TR
                off_val(s,r,q,k) = Y(k,q);
            end %TR (k)
        end %ROI (q)
    end %RUN (r)
end %SUBJ (s)


%% Off par Corr
for s=1:length(SUBJ_ID)
    for r=1:length(RUN)
        cd(PATHIN_OFF);
        load ([SUBJ_ID{s} '_time_SESS_' num2str(r)]);
        TR=1;
        for k=1:size(Y,1)-WINDSIZE
            off_corr(s,r,TR) = partialcorr(squeeze(Y(TR:TR+WINDSIZE,1)),squeeze(Y(TR:TR+WINDSIZE,2)),squeeze(Y(TR:TR+WINDSIZE,3)));
            TR=TR+1;
        end %TR
    end % SESS
end

