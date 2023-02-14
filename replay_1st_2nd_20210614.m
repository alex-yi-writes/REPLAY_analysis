%% REPLAY 1st and 2nd

clc;clear

% paths
paths           = [];
paths.parent    = 'D:\work_backups\replay_dat\';
paths.mri       = [paths.parent 'mri\'];
paths.behav     = [paths.parent 'behav\'];
paths.physio    = [paths.parent 'physio\'];
paths.app       = '\Applications\MRIcron.app\Contents\Resources\dcm2niix';
paths.coreg     = [paths.parent 'coreg\'];
paths.analysis  = [paths.parent 'mri\'];
paths.group    = [paths.parent '2nd\'];

% experimental info
% IDs = {'102'; '103'; '104'; '105'; '106'; '107'; '108'; '109'; '110';
%  '111'; '112'; '114'; '115'; '119'; '120'; '122'; '123};
IDs ={'138';'139';'140';'141';'142';'143';'145'};
TR=2.66;
slc=80;


behav=[]; onsets=[];
for id=1:length(IDs)
    
    % load
    %     for enc=1:2
    %         behav.enc{id,enc}  = load([paths.behav IDs{id} '_enc' num2str(enc) '.mat']);
    %     end
    if str2double(IDs{id})==132
        behav.ret1{id,1}=[];
        for ret1=2:4
            behav.ret1{id,ret1}  = load([paths.behav IDs{id} '_ret1' num2str(ret1) '.mat']);
        end
    else
        for ret1=1:4
            behav.ret1{id,ret1}  = load([paths.behav IDs{id} '_ret1' num2str(ret1) '.mat']);
        end
    end
    %     for ret2=1:2
    %         behav.ret2{id,ret2}  = load([paths.behav IDs{id} '_ret2' num2str(ret2) '.mat']);
    %     end
    
    
    % merge retrieval 1 behaviour files for processing fmri data
    if str2double(IDs{id})==132
        eval(['scanstart' num2str(s1) '=NaN;'])
        for s1=2:4
            eval(['clear scanstart' num2str(s1)])
            eval(['scanstart' num2str(s1) '=behav.ret1{id,' num2str(s1) '}.dat.retrieval1.results.SOT.raw.t0_fix0-(TR*5);']);
        end
        
    else
        for s1=1:4
            eval(['clear scanstart' num2str(s1)])
            eval(['scanstart' num2str(s1) '=behav.ret1{id,' num2str(s1) '}.dat.retrieval1.results.SOT.raw.t0_fix0-(TR*5);']);
        end
    end
    clear numvols
    if str2num(IDs{id})==107
        numvols= [210 180 210 210];
    elseif str2num(IDs{id})<=127 && str2num(IDs{id})~=107
        numvols = ones(1,4)*210;
        SessionMarker{id} = [ones(210,1); zeros(210,1); ones(210,1); zeros(210,1)];
    elseif str2num(IDs{id})>127
        numvols= [...
            length(spm_vol([paths.mri IDs{id} '\converted\run1.nii'])),...
            length(spm_vol([paths.mri IDs{id} '\converted\run2.nii'])),...
            length(spm_vol([paths.mri IDs{id} '\converted\run3.nii'])),...
            length(spm_vol([paths.mri IDs{id} '\converted\run4.nii']))];
        SessionMarker{id} = [ones(length(spm_vol([paths.mri IDs{id} '\converted\run1.nii'])),1); zeros(length(spm_vol([paths.mri IDs{id} '\converted\run2.nii'])),1);...
            ones(length(spm_vol([paths.mri IDs{id} '\converted\run3.nii'])),1); zeros(length(spm_vol([paths.mri IDs{id} '\converted\run4.nii'])),1)];
    end
    
    if str2double(IDs{id})==132
        onsets{id,1}.null           = [
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.null-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.null-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.null-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.cue_all        = [
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.cue-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.cue-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.cue-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.control        = [0;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.control-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.control-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.control-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.objects        = [
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.resp-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.resp-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.resp-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.feedback_all   = [
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.feedback-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.feedback-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.feedback-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.fixationX      = [
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.fix-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.fix-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.fix-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.confidence     = [
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.confidence-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.confidence-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.confidence-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        % task related onsets
        clear tmp1 tmp2 indctrl
        tmp1 = [
            cell2mat(behav.ret1{id,2}.dat.retrieval1.results.trl(:,2));...
            cell2mat(behav.ret1{id,3}.dat.retrieval1.results.trl(:,2));...
            cell2mat(behav.ret1{id,4}.dat.retrieval1.results.trl(:,2))];
        indctrl=tmp1==0;
        
        tmp2=[
            behav.ret1{id,2}.dat.retrieval1.results.accuracy;...
            behav.ret1{id,3}.dat.retrieval1.results.accuracy;...
            behav.ret1{id,4}.dat.retrieval1.results.accuracy];
        tmp2(indctrl==1)=[];
        
        tmp3 = tmp1; tmp3(indctrl)=[];
        
        
        onsets{id,1}.cue_correct        = onsets{id,1}.cue_all(tmp2==1);
        onsets{id,1}.cue_wrong          = onsets{id,1}.cue_all(tmp2~=1);
        
        onsets{id,1}.fb_correct         = onsets{id,1}.feedback_all(tmp2==1);
        onsets{id,1}.fb_wrong           = onsets{id,1}.feedback_all(tmp2~=1);
        
        onsets{id,1}.fb_correct_first   = onsets{id,1}.feedback_all(tmp2==1 & tmp3==1);
        onsets{id,1}.fb_wrong_first     = onsets{id,1}.feedback_all(tmp2~=1 & tmp3==1);
        
        onsets{id,1}.fb_correct_second  = onsets{id,1}.feedback_all(tmp2==1 & tmp3==2);
        onsets{id,1}.fb_wrong_second    = onsets{id,1}.feedback_all(tmp2~=1 & tmp3==2);
        
        onsets{id,1}.obj_correct        = onsets{id,1}.cue_correct+8;
        onsets{id,1}.obj_wrong          = onsets{id,1}.cue_wrong+8;
        
        onsets{id,1}.cue_correct_first  = onsets{id,1}.cue_all(tmp2==1 & tmp3==1);
        onsets{id,1}.cue_correct_second = onsets{id,1}.cue_all(tmp2==1 & tmp3==2);
        
        onsets{id,1}.cue_wrong_first    = onsets{id,1}.cue_all(tmp2~=1 & tmp3==1);
        onsets{id,1}.cue_wrong_second   = onsets{id,1}.cue_all(tmp2~=1 & tmp3==2);
        
        onsets{id,1}.obj_correct_first  = onsets{id,1}.cue_correct_first+8;
        onsets{id,1}.obj_correct_second = onsets{id,1}.cue_correct_second+8;
        
        onsets{id,1}.obj_wrong_first    = onsets{id,1}.cue_wrong_first+8;
        onsets{id,1}.obj_wrong_second   = onsets{id,1}.obj_correct_second+8;
        
        
        
    else
        onsets{id,1}.null           = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.null-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.null-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.null-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.null-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.cue_all        = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.cue-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.cue-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.cue-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.cue-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.control        = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.control-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.control-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.control-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.control-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.objects        = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.resp-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.resp-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.resp-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.resp-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.feedback_all   = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.feedback-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.feedback-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.feedback-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.feedback-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.fixationX      = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.fix-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.fix-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.fix-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.fix-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        onsets{id,1}.confidence     = [behav.ret1{id,1}.dat.retrieval1.results.SOT.raw.confidence-scanstart1;...
            (behav.ret1{id,2}.dat.retrieval1.results.SOT.raw.confidence-scanstart2)+(numvols(1)*TR);...
            (behav.ret1{id,3}.dat.retrieval1.results.SOT.raw.confidence-scanstart3)+(numvols(1)*TR+numvols(2)*TR);...
            (behav.ret1{id,4}.dat.retrieval1.results.SOT.raw.confidence-scanstart4)+(numvols(1)*TR+numvols(2)*TR+numvols(3)*TR)];
        
        % task related onsets
        clear tmp1 tmp2 indctrl
        tmp1 = [cell2mat(behav.ret1{id,1}.dat.retrieval1.results.trl(:,2));...
            cell2mat(behav.ret1{id,2}.dat.retrieval1.results.trl(:,2));...
            cell2mat(behav.ret1{id,3}.dat.retrieval1.results.trl(:,2));...
            cell2mat(behav.ret1{id,4}.dat.retrieval1.results.trl(:,2))];
        indctrl=tmp1==0;
        
        tmp2=[behav.ret1{id,1}.dat.retrieval1.results.accuracy;...
            behav.ret1{id,2}.dat.retrieval1.results.accuracy;...
            behav.ret1{id,3}.dat.retrieval1.results.accuracy;...
            behav.ret1{id,4}.dat.retrieval1.results.accuracy];
        tmp2(indctrl==1)=[];
        
        tmp3 = tmp1; tmp3(indctrl)=[];
        
        
        onsets{id,1}.cue_correct        = onsets{id,1}.cue_all(tmp2==1);
        onsets{id,1}.cue_wrong          = onsets{id,1}.cue_all(tmp2~=1);
        
        onsets{id,1}.fb_correct         = onsets{id,1}.feedback_all(tmp2==1);
        onsets{id,1}.fb_wrong           = onsets{id,1}.feedback_all(tmp2~=1);
        
        onsets{id,1}.fb_correct_first   = onsets{id,1}.feedback_all(tmp2==1 & tmp3==1);
        onsets{id,1}.fb_wrong_first     = onsets{id,1}.feedback_all(tmp2~=1 & tmp3==1);
        
        onsets{id,1}.fb_correct_second  = onsets{id,1}.feedback_all(tmp2==1 & tmp3==2);
        onsets{id,1}.fb_wrong_second    = onsets{id,1}.feedback_all(tmp2~=1 & tmp3==2);
        
        onsets{id,1}.obj_correct        = onsets{id,1}.cue_correct+8;
        onsets{id,1}.obj_wrong          = onsets{id,1}.cue_wrong+8;
        
        onsets{id,1}.cue_correct_first  = onsets{id,1}.cue_all(tmp2==1 & tmp3==1);
        onsets{id,1}.cue_correct_second = onsets{id,1}.cue_all(tmp2==1 & tmp3==2);
        
        onsets{id,1}.cue_wrong_first    = onsets{id,1}.cue_all(tmp2~=1 & tmp3==1);
        onsets{id,1}.cue_wrong_second   = onsets{id,1}.cue_all(tmp2~=1 & tmp3==2);
        
        onsets{id,1}.obj_correct_first  = onsets{id,1}.cue_correct_first+8;
        onsets{id,1}.obj_correct_second = onsets{id,1}.cue_correct_second+8;
        
        onsets{id,1}.obj_wrong_first    = onsets{id,1}.cue_wrong_first+8;
        onsets{id,1}.obj_wrong_second   = onsets{id,1}.obj_correct_second+8;
    end
    
end


%% build physio + realignment parameters
% id 102 has no physio data for run 1

for id=1:length(IDs)
    
    % read physio parameters
    if str2double(IDs{id})==132 || str2double(IDs{id})==137
        disp(['no physio for ID ' IDs{id} ', substituting it with CompCor'])
        physioFile     = [paths.physio IDs{id} '\raw\' IDs{id} '.txt'];
        physiodata     = importdata(physioFile);
%         physiofile=[]; clear tmp
%         for r=2:4
%             physioFile     = [paths.physio IDs{id} '\made\' IDs{id} '_physio' num2str(r) '.txt'];
%             tmp{r}         = importdata(physioFile);
%         end
%         tmp{1}= zeros(size(tmp{2}));
%         physiodata=[tmp{1};tmp{2};tmp{3};tmp{4}];
        
    else
        physiofile=[]; clear tmp
        for r=1:4
            physioFile     = [paths.physio IDs{id} '\made\' IDs{id} '_physio' num2str(r) '.txt'];
            tmp{r}         = importdata(physioFile);
        end
        physiodata=[tmp{1};tmp{2};tmp{3};tmp{4}];
        
    end
    
    % read realignment parameters
    realignFile = [paths.mri IDs{id} '\preproc\rp_afunc.txt'];
    realigndata = importdata(realignFile);
    
    % concatenate
    if str2num(IDs{id})==107
        SessionMarker107=[ones(210,1); zeros(180,1); ones(210,1); zeros(210,1)];
        R = [physiodata realigndata SessionMarker107];
    else
        R = [physiodata realigndata SessionMarker{id}];
    end
    
    save([paths.mri IDs{id} '\preproc\reg_all.mat'],'R');
    clear R realigndata realignFile physiodata physioFile
    
end

disp('done')

%% set up a model

for id=1:length(IDs)
    
    mkdir([paths.mri IDs{id} '\analysis\'])
    SPMdir=[paths.mri IDs{id} '\analysis\'];
    
    clear list_scan
    volnum = length(spm_vol([paths.mri IDs{id} '\preproc\s15rafunc.nii']));
    for v1 = 1:volnum
        list_scan{v1,1} = [paths.mri IDs{id} '\preproc\s15rafunc.nii,' num2str(v1)];
    end
    
    %% model specification
    
    clear model1
    
    spm_jobman('initcfg')
    
    model1{1}.spm.stats.fmri_spec.dir = {SPMdir};
    model1{1}.spm.stats.fmri_spec.timing.units = 'secs';
    model1{1}.spm.stats.fmri_spec.timing.RT = TR;
    model1{1}.spm.stats.fmri_spec.timing.fmri_t = slc;
    model1{1}.spm.stats.fmri_spec.timing.fmri_t0 = 40;
    model1{1}.spm.stats.fmri_spec.sess.scans = list_scan;
    
    if str2num(IDs{id})==125 || str2num(IDs{id})==129 % ID 125 and ID 129 does not have wrong items for the 2nd cue
        
    model1{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Cue_1st_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(1).onset = onsets{id,1}.cue_correct_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Cue_2nd_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.cue_correct_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Cue_1st_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(3).onset = onsets{id,1}.cue_wrong_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
    
%     model1{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Cue_2nd_Incorrect';
%     model1{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.cue_wrong_second;
%     model1{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
%     model1{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
%     model1{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
%     model1{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
    
    
    
    % feedbacks
    model1{1}.spm.stats.fmri_spec.sess.cond(4).name = 'FB_1st_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(4).onset = onsets{id,1}.fb_correct_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(5).name = 'FB_2nd_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(5).onset= onsets{id,1}.fb_correct_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(6).name = 'FB_1st_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(6).onset = onsets{id,1}.fb_wrong_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
    
%     model1{1}.spm.stats.fmri_spec.sess.cond(8).name = 'FB_2nd_Incorrect';
%     model1{1}.spm.stats.fmri_spec.sess.cond(8).onset= onsets{id,1}.fb_wrong_second;
%     model1{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
%     model1{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
%     model1{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
%     model1{1}.spm.stats.fmri_spec.sess.cond(8).orth = 1;
    
    
    
    % objects
    model1{1}.spm.stats.fmri_spec.sess.cond(7).name = 'Objects_1st_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(7).onset = onsets{id,1}.obj_correct_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(7).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(7).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(8).name = 'Objects_2nd_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(8).onset= onsets{id,1}.obj_correct_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(8).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(9).name = 'Objects_1st_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(9).onset = onsets{id,1}.obj_wrong_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(9).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(9).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(9).orth = 1;
    
%     model1{1}.spm.stats.fmri_spec.sess.cond(12).name = 'Objects_2nd_Incorrect';
%     model1{1}.spm.stats.fmri_spec.sess.cond(12).onset= onsets{id,1}.obj_wrong_second;
%     model1{1}.spm.stats.fmri_spec.sess.cond(12).duration = 0;
%     model1{1}.spm.stats.fmri_spec.sess.cond(12).tmod = 0;
%     model1{1}.spm.stats.fmri_spec.sess.cond(12).pmod = struct('name', {}, 'param', {}, 'poly', {});
%     model1{1}.spm.stats.fmri_spec.sess.cond(12).orth = 1;
    
    
    
    model1{1}.spm.stats.fmri_spec.sess.cond(10).name = 'Controls';
    model1{1}.spm.stats.fmri_spec.sess.cond(10).onset= onsets{id,1}.control;
    model1{1}.spm.stats.fmri_spec.sess.cond(10).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(10).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(10).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(10).orth = 1;
    
    
    model1{1}.spm.stats.fmri_spec.sess.cond(11).name = 'FixationX';
    model1{1}.spm.stats.fmri_spec.sess.cond(11).onset= onsets{id,1}.fixationX;
    model1{1}.spm.stats.fmri_spec.sess.cond(11).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(11).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(11).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(11).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(12).name = 'confidence';
    model1{1}.spm.stats.fmri_spec.sess.cond(12).onset= onsets{id,1}.confidence;
    model1{1}.spm.stats.fmri_spec.sess.cond(12).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(12).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(12).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(12).orth = 1;
   
        
    else
    
    % cues
    model1{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Cue_1st_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(1).onset = onsets{id,1}.cue_correct_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Cue_2nd_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(2).onset= onsets{id,1}.cue_correct_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(3).name = 'Cue_1st_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(3).onset = onsets{id,1}.cue_wrong_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(4).name = 'Cue_2nd_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(4).onset= onsets{id,1}.cue_wrong_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
    
    
    
    % feedbacks
    model1{1}.spm.stats.fmri_spec.sess.cond(5).name = 'FB_1st_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(5).onset = onsets{id,1}.fb_correct_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(5).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(6).name = 'FB_2nd_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(6).onset= onsets{id,1}.fb_correct_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(7).name = 'FB_1st_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(7).onset = onsets{id,1}.fb_wrong_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(7).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(7).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(7).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(8).name = 'FB_2nd_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(8).onset= onsets{id,1}.fb_wrong_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(8).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(8).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(8).orth = 1;
    
    
    
    % objects
    model1{1}.spm.stats.fmri_spec.sess.cond(9).name = 'Objects_1st_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(9).onset = onsets{id,1}.obj_correct_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(9).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(9).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(9).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(9).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(10).name = 'Objects_2nd_Correct';
    model1{1}.spm.stats.fmri_spec.sess.cond(10).onset= onsets{id,1}.obj_correct_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(10).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(10).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(10).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(10).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(11).name = 'Objects_1st_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(11).onset = onsets{id,1}.obj_wrong_first;
    model1{1}.spm.stats.fmri_spec.sess.cond(11).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(11).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(11).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(11).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(12).name = 'Objects_2nd_Incorrect';
    model1{1}.spm.stats.fmri_spec.sess.cond(12).onset= onsets{id,1}.obj_wrong_second;
    model1{1}.spm.stats.fmri_spec.sess.cond(12).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(12).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(12).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(12).orth = 1;
    
    
    
    model1{1}.spm.stats.fmri_spec.sess.cond(13).name = 'Controls';
    model1{1}.spm.stats.fmri_spec.sess.cond(13).onset= onsets{id,1}.control;
    model1{1}.spm.stats.fmri_spec.sess.cond(13).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(13).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(13).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(13).orth = 1;
    
    
    model1{1}.spm.stats.fmri_spec.sess.cond(14).name = 'FixationX';
    model1{1}.spm.stats.fmri_spec.sess.cond(14).onset= onsets{id,1}.fixationX;
    model1{1}.spm.stats.fmri_spec.sess.cond(14).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(14).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(14).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(14).orth = 1;
    
    model1{1}.spm.stats.fmri_spec.sess.cond(15).name = 'confidence';
    model1{1}.spm.stats.fmri_spec.sess.cond(15).onset= onsets{id,1}.confidence;
    model1{1}.spm.stats.fmri_spec.sess.cond(15).duration = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(15).tmod = 0;
    model1{1}.spm.stats.fmri_spec.sess.cond(15).pmod = struct('name', {}, 'param', {}, 'poly', {});
    model1{1}.spm.stats.fmri_spec.sess.cond(15).orth = 1;
    
    end
    
    model1{1}.spm.stats.fmri_spec.sess.multi = {''};
    model1{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
%     if str2num(IDs{id})==102
%     model1{1}.spm.stats.fmri_spec.sess.multi_reg = {[paths.mri IDs{id} '\preproc\rp_afunc.txt']};
%     else
    model1{1}.spm.stats.fmri_spec.sess.multi_reg = {[paths.mri IDs{id} '\preproc\reg_all.mat']};
%     end
    model1{1}.spm.stats.fmri_spec.sess.hpf = 128;
    model1{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    model1{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    model1{1}.spm.stats.fmri_spec.volt = 1;
    model1{1}.spm.stats.fmri_spec.global = 'None';
    model1{1}.spm.stats.fmri_spec.mthresh = 0.8;
    model1{1}.spm.stats.fmri_spec.mask = {''};
    model1{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    spm_jobman('run',model1)
    
    
    %% estimate: can take some time
    
    clear model1
    
    spm_jobman('initcfg')
    model1{1}.spm.stats.fmri_est.spmmat = {[SPMdir 'SPM.mat']};
    model1{1}.spm.stats.fmri_est.write_residuals = 0;
    model1{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run',model1)
    
    
    %% statistical inference
    
    if str2num(IDs{id})==125 % ID 125 does not have wrong items for the 2nd cue
        
        
        % contrast As
    cueCorrect_v_control         = [1 1 zeros(1,7) -2];
    cueCorrect1st_v_control      = [1 0 zeros(1,7) -1];
    cueCorrect2nd_v_control      = [0 1 zeros(1,7) -1];
    cueCorrect1st_v_cueCorrect2nd= [1 -1];
    
    control_v_cueCorrect         = [-1 -1 zeros(1,7) 2];
    control_v_cueCorrect1st      = [-1 0 zeros(1,7) 1];
    control_v_cueCorrect2nd      = [0 -1 zeros(1,7) 1];
    cueCorrect2nd_v_cueCorrect1st= [-1 1];
    
    
    % contrast Bs
    cueCorrect_v_cueIncorrect    = [1 1 -2];
    cueCorrect1st_v_cueIncorrect1st = [1 0 -1];
    cueCorrect2nd_v_cueIncorrect2nd = [1]; % placeholder
    cue1st_v_cue2nd              = [1 -2 1];
    
    cueIncorrect_v_cueCorrect    = [-1 -1 2];
    cueIncorrect1st_v_cueCorrect1st = [-1 0 1];
    cueIncorrect2nd_v_cueCorrect2nd = [1]; % placeholder
    cue2nd_v_cue1st              = [-1 2 -1];
    
    
    % contrast Cs
    FbCorrect_v_CueCorrect       = [-1 -1 0 1 1];
    FbIncorrect_v_CueIncorrect   = [0 0 -1 0 0 1];
    FbCorrect_v_FbIncorrect      = [0 0 0 1 1 -2];
    FbCorrect_v_control          = [0 0 0 1 1 zeros(1,4) -2];
    
    FbCorrect1st_v_FbCorrect2nd  = [0 0 0 1 -1];
    FbIncorrect1st_v_FbCorrect2nd= [0 0 0 0 -1 1];
    
    CueCorrect_v_FbCorrect       = [1 1 0 -1 -1];
    CueIncorrect_v_FbIncorrect   = [0 0 1 0 0 -1];
    FbIncorrect_v_FbCorrect      = [0 0 0 -1 -1 2];
    control_v_FbCorrect          = [0 0 0 -1 -1 zeros(1,4) 2];
    
    FbCorrect2nd_v_FbCorrect1st  = [0 0 0 -1 1];
    FbIncorrect2nd_v_FbCorrect1st= [1]; % placeholder
    
    
    % contrast Objects
    % contrast AO
    objCorrect_v_control         = [zeros(1,7) 1 1 0 0 -2];
    objCorrect1st_v_control      = [zeros(1,7) 1 0 0 0 -1];
    objCorrect2nd_v_control      = [zeros(1,7) 0 1 0 0 -1];
    objCorrect1st_v_objCorrect2nd= [zeros(1,7) 1 -1];
    
    control_v_objCorrect         = [zeros(1,7) -1 -1 0 0 2];
    control_v_objCorrect1st      = [zeros(1,7) -1 0 0 0 1];
    control_v_objCorrect2nd      = [zeros(1,7) 0 -1 0 0 1];
    objCorrect2nd_v_objCorrect1st= [zeros(1,7) -1 1];
    
    
    % contrast BO
    objCorrect_v_objIncorrect    = [zeros(1,7) 1 1 -2];
    objCorrect1st_v_objIncorrect1st = [zeros(1,7) 1 0 -1];
    objCorrect2nd_v_objIncorrect2nd = [1]; % placeholder
    obj1st_v_obj2nd              = [zeros(1,7) 1 1 -2];
    
    objIncorrect_v_objCorrect    = [zeros(1,7) -1 -1 2];
    objIncorrect1st_v_objCorrect1st = [zeros(1,7) -1 0 1];
    objIncorrect2nd_v_objCorrect2nd = [1]; % placeholder
    obj2nd_v_obj1st              = [zeros(1,7) -1 -1 2];
        
    else
    
    % contrast As
    cueCorrect_v_control         = [1 1 zeros(1,10) -2];
    cueCorrect1st_v_control      = [1 0 zeros(1,10) -1];
    cueCorrect2nd_v_control      = [0 1 zeros(1,10) -1];
    cueCorrect1st_v_cueCorrect2nd= [1 -1];
    
    control_v_cueCorrect         = [-1 -1 zeros(1,10) 2];
    control_v_cueCorrect1st      = [-1 0 zeros(1,10) 1];
    control_v_cueCorrect2nd      = [0 -1 zeros(1,10) 1];
    cueCorrect2nd_v_cueCorrect1st= [-1 1];
    
    
    % contrast Bs
    cueCorrect_v_cueIncorrect    = [1 1 -1 -1];
    cueCorrect1st_v_cueIncorrect1st = [1 0 -1];
    cueCorrect2nd_v_cueIncorrect2nd = [0 1 0 -1];
    cue1st_v_cue2nd              = [1 -1 1 -1];
    
    cueIncorrect_v_cueCorrect    = [-1 -1 1 1];
    cueIncorrect1st_v_cueCorrect1st = [-1 0 1];
    cueIncorrect2nd_v_cueCorrect2nd = [0 -1 0 1];
    cue2nd_v_cue1st              = [-1 1 -1 1];
    
    
    % contrast Cs
    FbCorrect_v_CueCorrect       = [-1 -1 0 0 1 1];
    FbIncorrect_v_CueIncorrect   = [0 0 -1 -1 0 0 1 1];
    FbCorrect_v_FbIncorrect      = [0 0 0 0 1 1 -1 -1];
    FbCorrect_v_control          = [0 0 0 0 1 1 zeros(1,6) -2];
    
    FbCorrect1st_v_FbCorrect2nd  = [0 0 0 0 1 -1];
    FbIncorrect1st_v_FbCorrect2nd= [0 0 0 0 0 -1 1];
    
    CueCorrect_v_FbCorrect       = [1 1 0 0 -1 -1];
    CueIncorrect_v_FbIncorrect   = [0 0 1 1 0 0 -1 -1];
    FbIncorrect_v_FbCorrect      = [0 0 0 0 -1 -1 1 1];
    control_v_FbCorrect          = [0 0 0 0 -1 -1 zeros(1,6) 2];
    
    FbCorrect2nd_v_FbCorrect1st  = [0 0 0 0 -1 1];
    FbIncorrect2nd_v_FbCorrect1st= [0 0 0 0 0 1 -1];
    
    
    % contrast Objects
    % contrast AO
    objCorrect_v_control         = [zeros(1,8) 1 1 0 0 -2];
    objCorrect1st_v_control      = [zeros(1,8) 1 0 0 0 -1];
    objCorrect2nd_v_control      = [zeros(1,8) 0 1 0 0 -1];
    objCorrect1st_v_objCorrect2nd= [zeros(1,8) 1 -1];
    
    control_v_objCorrect         = [zeros(1,8) -1 -1 0 0 2];
    control_v_objCorrect1st      = [zeros(1,8) -1 0 0 0 1];
    control_v_objCorrect2nd      = [zeros(1,8) 0 -1 0 0 1];
    objCorrect2nd_v_objCorrect1st= [zeros(1,8) -1 1];
    
    
    % contrast BO
    objCorrect_v_objIncorrect    = [zeros(1,8) 1 1 -1 -1];
    objCorrect1st_v_objIncorrect1st = [zeros(1,8) 1 0 -1];
    objCorrect2nd_v_objIncorrect2nd = [zeros(1,8) 0 1 0 -1];
    obj1st_v_obj2nd              = [zeros(1,8) 1 1 -1 -1];
    
    objIncorrect_v_objCorrect    = [zeros(1,8) -1 -1 1 1];
    objIncorrect1st_v_objCorrect1st = [zeros(1,8) -1 0 1];
    objIncorrect2nd_v_objCorrect2nd = [zeros(1,8) 0 -1 0 1];
    obj2nd_v_obj1st              = [zeros(1,8) -1 -1 1 1];
    
    end
    
    
    % compute
    clear ReplayStats
    
    spm_jobman('initcfg');      % initiate job manager
    
    ReplayStats{1}.spm.stats.con.spmmat = cellstr([paths.analysis IDs{id} '\analysis\SPM.mat']);
    
    
    % contrast As
    
    ReplayStats{1}.spm.stats.con.consess{1}.tcon.name = 'CueCorrectAll > control';
    ReplayStats{1}.spm.stats.con.consess{1}.tcon.weights = cueCorrect_v_control;
    ReplayStats{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{2}.tcon.name = 'CueCorrect1st > control';
    ReplayStats{1}.spm.stats.con.consess{2}.tcon.weights = cueCorrect1st_v_control;
    ReplayStats{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{3}.tcon.name = 'CueCorrect2nd > control';
    ReplayStats{1}.spm.stats.con.consess{3}.tcon.weights = cueCorrect2nd_v_control;
    ReplayStats{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{4}.tcon.name = 'CueCorrect1st > CueCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{4}.tcon.weights = cueCorrect1st_v_cueCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    % contrast As inv
    
    ReplayStats{1}.spm.stats.con.consess{5}.tcon.name = 'control > CueCorrectAll';
    ReplayStats{1}.spm.stats.con.consess{5}.tcon.weights = control_v_cueCorrect;
    ReplayStats{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{6}.tcon.name = 'control > CueCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{6}.tcon.weights = control_v_cueCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{7}.tcon.name = 'control > CueCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{7}.tcon.weights = control_v_cueCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{8}.tcon.name = 'CueCorrect2nd > CueCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{8}.tcon.weights = cueCorrect2nd_v_cueCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    
    
    % contrast Bs
    
    ReplayStats{1}.spm.stats.con.consess{9}.tcon.name = 'cueCorrectAll > CueIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{9}.tcon.weights = cueCorrect_v_cueIncorrect;
    ReplayStats{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{10}.tcon.name = 'cueCorrect1st > CueIncorrect1st';
    ReplayStats{1}.spm.stats.con.consess{10}.tcon.weights = cueCorrect1st_v_cueIncorrect1st;
    ReplayStats{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
    
    if str2num(IDs{id})==125
    ReplayStats{1}.spm.stats.con.consess{11}.tcon.name = 'null';
    ReplayStats{1}.spm.stats.con.consess{11}.tcon.weights = cueCorrect2nd_v_cueIncorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
    else
    ReplayStats{1}.spm.stats.con.consess{11}.tcon.name = 'cueCorrect2nd > CueIncorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{11}.tcon.weights = cueCorrect2nd_v_cueIncorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
    end
    
    ReplayStats{1}.spm.stats.con.consess{12}.tcon.name = 'Cue1st > Cue2nd';
    ReplayStats{1}.spm.stats.con.consess{12}.tcon.weights = cue1st_v_cue2nd;
    ReplayStats{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
    % contrast Bs inv
    
    ReplayStats{1}.spm.stats.con.consess{13}.tcon.name = 'cueIncorrectAll > CueIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{13}.tcon.weights = cueIncorrect_v_cueCorrect;
    ReplayStats{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{14}.tcon.name = 'cueIncorrect1st > CueCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{14}.tcon.weights = cueIncorrect1st_v_cueCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
    
    if str2num(IDs{id})==125
    ReplayStats{1}.spm.stats.con.consess{15}.tcon.name = 'null';
    ReplayStats{1}.spm.stats.con.consess{15}.tcon.weights = cueIncorrect2nd_v_cueCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    else
    ReplayStats{1}.spm.stats.con.consess{15}.tcon.name = 'cueIncorrect2nd > CueCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{15}.tcon.weights = cueIncorrect2nd_v_cueCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    end
    
    ReplayStats{1}.spm.stats.con.consess{16}.tcon.name = 'Cue2nd > Cue1st';
    ReplayStats{1}.spm.stats.con.consess{16}.tcon.weights = cue2nd_v_cue1st;
    ReplayStats{1}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
        
    % contrast Cs
    
    ReplayStats{1}.spm.stats.con.consess{17}.tcon.name = 'FBCorrectAll > CueCorrectAll';
    ReplayStats{1}.spm.stats.con.consess{17}.tcon.weights = FbCorrect_v_CueCorrect;
    ReplayStats{1}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{18}.tcon.name = 'FBIncorrectAll > CueIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{18}.tcon.weights = FbIncorrect_v_CueIncorrect;
    ReplayStats{1}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{19}.tcon.name = 'FBCorrectall > FBIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{19}.tcon.weights = FbCorrect_v_FbIncorrect;
    ReplayStats{1}.spm.stats.con.consess{19}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{20}.tcon.name = 'FBCorrectAll > Control';
    ReplayStats{1}.spm.stats.con.consess{20}.tcon.weights = FbCorrect_v_control;
    ReplayStats{1}.spm.stats.con.consess{20}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{21}.tcon.name = 'FBCorrect1st > FBCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{21}.tcon.weights = FbCorrect1st_v_FbCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{21}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{22}.tcon.name = 'FBIncorrect1st > FBCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{22}.tcon.weights = FbIncorrect1st_v_FbCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{22}.tcon.sessrep = 'none';
    
    
    % contrast Cs inv
    
    ReplayStats{1}.spm.stats.con.consess{23}.tcon.name = 'CueCorrectAll > FBCorrectAll';
    ReplayStats{1}.spm.stats.con.consess{23}.tcon.weights = CueCorrect_v_FbCorrect;
    ReplayStats{1}.spm.stats.con.consess{23}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{24}.tcon.name = 'CueIncorrectAll > FBIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{24}.tcon.weights = CueIncorrect_v_FbIncorrect;
    ReplayStats{1}.spm.stats.con.consess{24}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{25}.tcon.name = 'FBIncorrectall > FBCorrectAll';
    ReplayStats{1}.spm.stats.con.consess{25}.tcon.weights = FbIncorrect_v_FbCorrect;
    ReplayStats{1}.spm.stats.con.consess{25}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{26}.tcon.name = 'Control > FBCorrectAll';
    ReplayStats{1}.spm.stats.con.consess{26}.tcon.weights = control_v_FbCorrect;
    ReplayStats{1}.spm.stats.con.consess{26}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{27}.tcon.name = 'FBCorrect2nd > FBCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{27}.tcon.weights = FbCorrect2nd_v_FbCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{27}.tcon.sessrep = 'none';
    
    if str2num(IDs{id})==125
    ReplayStats{1}.spm.stats.con.consess{28}.tcon.name = 'null';
    ReplayStats{1}.spm.stats.con.consess{28}.tcon.weights = FbIncorrect2nd_v_FbCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{28}.tcon.sessrep = 'none';
    else
    ReplayStats{1}.spm.stats.con.consess{28}.tcon.name = 'FBCorrect2nd > FBIncorrect1st';
    ReplayStats{1}.spm.stats.con.consess{28}.tcon.weights = FbIncorrect2nd_v_FbCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{28}.tcon.sessrep = 'none';
    end
    
    
    % contrast AO
    
    ReplayStats{1}.spm.stats.con.consess{29}.tcon.name = 'ObjCorrectAll > control';
    ReplayStats{1}.spm.stats.con.consess{29}.tcon.weights = objCorrect_v_control;
    ReplayStats{1}.spm.stats.con.consess{29}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{30}.tcon.name = 'ObjCorrect1st > control';
    ReplayStats{1}.spm.stats.con.consess{30}.tcon.weights = objCorrect1st_v_control;
    ReplayStats{1}.spm.stats.con.consess{30}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{31}.tcon.name = 'ObjCorrect2nd > control';
    ReplayStats{1}.spm.stats.con.consess{31}.tcon.weights = objCorrect2nd_v_control;
    ReplayStats{1}.spm.stats.con.consess{31}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{32}.tcon.name = 'ObjCorrect1st > ObjCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{32}.tcon.weights = objCorrect1st_v_objCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{32}.tcon.sessrep = 'none';
    
    % contrast AO inv
    
    ReplayStats{1}.spm.stats.con.consess{33}.tcon.name = 'control > ObjCorrectAll';
    ReplayStats{1}.spm.stats.con.consess{33}.tcon.weights = control_v_objCorrect;
    ReplayStats{1}.spm.stats.con.consess{33}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{34}.tcon.name = 'control > ObjCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{34}.tcon.weights = control_v_objCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{34}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{35}.tcon.name = 'control > ObjCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{35}.tcon.weights = control_v_objCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{35}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.consess{36}.tcon.name = 'ObjCorrect2nd > ObjCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{36}.tcon.weights = objCorrect2nd_v_objCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{36}.tcon.sessrep = 'none';
    
    
    
    % contrast BO
    
    ReplayStats{1}.spm.stats.con.consess{37}.tcon.name = 'objCorrectAll > ObjIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{37}.tcon.weights = objCorrect_v_objIncorrect;
    ReplayStats{1}.spm.stats.con.consess{37}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{38}.tcon.name = 'objCorrect1st > ObjIncorrect1st';
    ReplayStats{1}.spm.stats.con.consess{38}.tcon.weights = objCorrect1st_v_objIncorrect1st;
    ReplayStats{1}.spm.stats.con.consess{38}.tcon.sessrep = 'none';
    if str2num(IDs{id})==125
    ReplayStats{1}.spm.stats.con.consess{39}.tcon.name = 'null';
    ReplayStats{1}.spm.stats.con.consess{39}.tcon.weights = objCorrect2nd_v_objIncorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{39}.tcon.sessrep = 'none';    
    else
    ReplayStats{1}.spm.stats.con.consess{39}.tcon.name = 'objCorrect2nd > ObjIncorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{39}.tcon.weights = objCorrect2nd_v_objIncorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{39}.tcon.sessrep = 'none';
    end
    
    ReplayStats{1}.spm.stats.con.consess{40}.tcon.name = 'Obj1st > Obj2nd';
    ReplayStats{1}.spm.stats.con.consess{40}.tcon.weights = obj1st_v_obj2nd;
    ReplayStats{1}.spm.stats.con.consess{40}.tcon.sessrep = 'none';
    
    % contrast BO inv
    
    ReplayStats{1}.spm.stats.con.consess{41}.tcon.name = 'objIncorrectAll > ObjIncorrectAll';
    ReplayStats{1}.spm.stats.con.consess{41}.tcon.weights = objIncorrect_v_objCorrect;
    ReplayStats{1}.spm.stats.con.consess{41}.tcon.sessrep = 'none';
   
    ReplayStats{1}.spm.stats.con.consess{42}.tcon.name = 'objIncorrect1st > ObjCorrect1st';
    ReplayStats{1}.spm.stats.con.consess{42}.tcon.weights = objIncorrect1st_v_objCorrect1st;
    ReplayStats{1}.spm.stats.con.consess{42}.tcon.sessrep = 'none';
    
    if str2num(IDs{id})==125
    ReplayStats{1}.spm.stats.con.consess{43}.tcon.name = 'null';
    ReplayStats{1}.spm.stats.con.consess{43}.tcon.weights = objIncorrect2nd_v_objCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{43}.tcon.sessrep = 'none';
    else
    ReplayStats{1}.spm.stats.con.consess{43}.tcon.name = 'objIncorrect2nd > ObjCorrect2nd';
    ReplayStats{1}.spm.stats.con.consess{43}.tcon.weights = objIncorrect2nd_v_objCorrect2nd;
    ReplayStats{1}.spm.stats.con.consess{43}.tcon.sessrep = 'none';    
    end
    
    ReplayStats{1}.spm.stats.con.consess{44}.tcon.name = 'Obj2nd > Obj1st';
    ReplayStats{1}.spm.stats.con.consess{44}.tcon.weights = obj2nd_v_obj1st;
    ReplayStats{1}.spm.stats.con.consess{44}.tcon.sessrep = 'none';
    
    ReplayStats{1}.spm.stats.con.delete = 1;    % 1 means delete, 0 means append
    
    spm_jobman('run', ReplayStats) % run batch
    
    
    %% copy files to a folder for coregistration
    mkdir([paths.coreg IDs{id} '\data\'])
    for ctr=1:44
        copyfile([paths.analysis IDs{id} '\analysis\' strcat('con_',sprintf('%04d',ctr),'.nii')],...
            [paths.coreg IDs{id} '\data\' strcat('con_',sprintf('%04d',ctr),'.nii')])
    end
    
end


%% group level analyses

%% make lists of contrasts


% ------- contrast As ------- %

list_cueCorrect_v_control = [];
for id = 1:length(IDs)
    list_cueCorrect_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0001_mni.nii,1'];
end
clear id

list_cueCorrect1st_v_control = []; 
for id = 1:length(IDs)
    list_cueCorrect1st_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0002_mni.nii,1'];
end
clear id

list_cueCorrect2nd_v_control = []; 
for id = 1:length(IDs)
    list_cueCorrect2nd_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0003_mni.nii,1'];
end
clear id

list_cueCorrect1st_v_cueCorrect2nd = []; 
for id = 1:length(IDs)
    list_cueCorrect1st_v_cueCorrect2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0004_mni.nii,1'];
end
clear id

% --------------------------- %




% ------- contrast Bs ------- %

list_cueCorrect_v_cueIncorrect = []; 
for id = 1:length(IDs)
    list_cueCorrect_v_cueIncorrect{id,1} = [paths.analysis IDs{id} '\analysis\con_0009_mni.nii,1'];
end
clear id

list_cueCorrect1st_v_cueIncorrect1st = []; 
for id = 1:length(IDs)
    list_cueCorrect1st_v_cueIncorrect1st{id,1} = [paths.analysis IDs{id} '\analysis\con_0010_mni.nii,1'];
end
clear id

list_cueCorrect2nd_v_cueIncorrect2nd = []; 
for id = 1:length(IDs)
    list_cueCorrect2nd_v_cueIncorrect2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0011_mni.nii,1'];
end
clear id

list_cue1st_v_cue2nd = []; 
for id = 1:length(IDs)
    list_cue1st_v_cue2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0012_mni.nii,1'];
end
clear id

% --------------------------- %



% ------- contrast Cs ------- %

list_FbCorrect_v_CueCorrect = []; 
for id = 1:length(IDs)
    list_FbCorrect_v_CueCorrect{id,1} = [paths.analysis IDs{id} '\analysis\con_0017_mni.nii,1'];
end
clear id

list_FbIncorrect_v_CueIncorrect = []; 
for id = 1:length(IDs)
    list_FbIncorrect_v_CueIncorrect{id,1} = [paths.analysis IDs{id} '\analysis\con_0018_mni.nii,1'];
end
clear id

list_FbCorrect_v_FbIncorrect = []; 
for id = 1:length(IDs)
    list_FbCorrect_v_FbIncorrect{id,1} = [paths.analysis IDs{id} '\analysis\con_0019_mni.nii,1'];
end
clear id

list_FbCorrect_v_control = []; 
for id = 1:length(IDs)
    list_FbCorrect_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0020_mni.nii,1'];
end
clear id

list_FbCorrect1st_v_FbCorrect2nd = []; 
for id = 1:length(IDs)
    list_FbCorrect1st_v_FbCorrect2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0021_mni.nii,1'];
end
clear id

list_FbIncorrect1st_v_FbCorrect2nd = []; 
for id = 1:length(IDs)
    list_FbIncorrect1st_v_FbCorrect2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0022_mni.nii,1'];
end
clear id

% --------------------------- %



% ------- contrast AOs ------- %

list_objCorrect_v_control = []; 
for id = 1:length(IDs)
    list_objCorrect_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0029_mni.nii,1'];
end
clear id

list_objCorrect1st_v_control = []; 
for id = 1:length(IDs)
    list_objCorrect1st_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0030_mni.nii,1'];
end
clear id

list_objCorrect2nd_v_control = []; 
for id = 1:length(IDs)
    list_objCorrect2nd_v_control{id,1} = [paths.analysis IDs{id} '\analysis\con_0031_mni.nii,1'];
end
clear id

list_objCorrect1st_v_objCorrect2nd = []; 
for id = 1:length(IDs)
    list_objCorrect1st_v_objCorrect2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0032_mni.nii,1'];
end
clear id

% --------------------------- %



% ------- contrast BOs ------- %

list_objCorrect_v_objIncorrect = []; 
for id = 1:length(IDs)
    list_objCorrect_v_objIncorrect{id,1} = [paths.analysis IDs{id} '\analysis\con_0037_mni.nii,1'];
end
clear id

list_objCorrect1st_v_objIncorrect1st = [];
for id = 1:length(IDs)
    list_objCorrect1st_v_objIncorrect1st{id,1} = [paths.analysis IDs{id} '\analysis\con_0038_mni.nii,1'];
end
clear id

list_objCorrect2nd_v_objIncorrect2nd = []; 
for id = 1:length(IDs)
    list_objCorrect2nd_v_objIncorrect2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0039_mni.nii,1'];
end
clear id

list_obj1st_v_obj2nd = []; 
for id = 1:length(IDs)
    list_obj1st_v_obj2nd{id,1} = [paths.analysis IDs{id} '\analysis\con_0040_mni.nii,1'];
end
clear id

% --------------------------- %


%% compute

%% contrast As

% ------- compute: CueCorrectAll > control ------- %

cd(paths.group); mkdir('cueCorrect_v_control'); cd cueCorrect_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'CueCorrectAll > control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: CueCorrect1st > control ------- %

cd(paths.group); mkdir('cueCorrect1st_v_control'); cd cueCorrect1st_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect1st_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'CueCorrect1st > control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: CueCorrect2nd > control ------- %

cd(paths.group); mkdir('cueCorrect2nd_v_control'); cd cueCorrect2nd_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect2nd_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'CueCorrect2nd > control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %




% ------- compute: CueCorrect1st > CueCorrect2nd ------- %

cd(paths.group); mkdir('cueCorrect1st_v_cueCorrect2nd'); cd cueCorrect1st_v_cueCorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect1st_v_cueCorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'CueCorrect1st > CueCorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



%% contrast Bs


% ------- compute: cueCorrectAll > CueIncorrectAll ------- %

cd(paths.group); mkdir('cueCorrect_v_cueIncorrect'); cd cueCorrect_v_cueIncorrect; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect_v_cueIncorrect;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'cueCorrectAll > CueIncorrectAll';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %


% ------- compute: cueCorrect1st > CueIncorrect1st ------- %

cd(paths.group); mkdir('cueCorrect1st_v_cueIncorrect1st'); cd cueCorrect1st_v_cueIncorrect1st; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect1st_v_cueIncorrect1st;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'cueCorrect1st > CueIncorrect1st';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: cueCorrect1st > CueIncorrect1st ------- %

cd(paths.group); mkdir('cueCorrect2nd_v_cueIncorrect2nd'); cd cueCorrect2nd_v_cueIncorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cueCorrect2nd_v_cueIncorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'cueCorrect2nd > CueIncorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: cueCorrect1st > CueIncorrect1st ------- %

cd(paths.group); mkdir('cue1st_v_cue2nd'); cd cue1st_v_cue2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_cue1st_v_cue2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Cue1st > Cue2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %


%% contrast Cs


% ------- compute: FBCorrectAll > CueCorrectAll ------- %

cd(paths.group); mkdir('FbCorrect_v_CueCorrect'); cd FbCorrect_v_CueCorrect; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbCorrect_v_CueCorrect;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBCorrectAll > CueCorrectAll';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: FBIncorrectAll > CueIncorrectAll ------- %

cd(paths.group); mkdir('FbIncorrect_v_CueIncorrect'); cd FbIncorrect_v_CueIncorrect; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbIncorrect_v_CueIncorrect;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBIncorrectAll > CueIncorrectAll';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: FBCorrectall > FBIncorrectAll ------- %

cd(paths.group); mkdir('FbCorrect_v_FbIncorrect'); cd FbCorrect_v_FbIncorrect; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbCorrect_v_FbIncorrect;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBCorrectall > FBIncorrectAll';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: FBCorrectAll > Control ------- %

cd(paths.group); mkdir('FbCorrect_v_control'); cd FbCorrect_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbCorrect_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBCorrectAll > Control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: FBCorrect1st > FBCorrect2nd ------- %

cd(paths.group); mkdir('FbCorrect1st_v_FbCorrect2nd'); cd FbCorrect1st_v_FbCorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbCorrect1st_v_FbCorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBCorrect1st > FBCorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: FBIncorrect1st > FBCorrect2nd ------- %

cd(paths.group); mkdir('FbIncorrect1st_v_FbCorrect2nd'); cd FbIncorrect1st_v_FbCorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbCorrect1st_v_FbCorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBIncorrect1st > FBCorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %


%% contrast AOs


% ------- compute: FBIncorrect1st > FBCorrect2nd ------- %

cd(paths.group); mkdir('FbIncorrect1st_v_FbCorrect2nd'); cd FbIncorrect1st_v_FbCorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_FbCorrect1st_v_FbCorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'FBIncorrect1st > FBCorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %


% ------- compute: ObjCorrectAll > control ------- %

cd(paths.group); mkdir('objCorrect_v_control'); cd objCorrect_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'ObjCorrectAll > control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: ObjCorrect1st > control ------- %

cd(paths.group); mkdir('objCorrect1st_v_control'); cd objCorrect1st_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect1st_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'ObjCorrect1st > control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: ObjCorrect2nd > control ------- %

cd(paths.group); mkdir('objCorrect2nd_v_control'); cd objCorrect2nd_v_control; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect2nd_v_control;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'ObjCorrect2nd > control';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %




% ------- compute: ObjCorrect1st > ObjCorrect2nd ------- %

cd(paths.group); mkdir('objCorrect1st_v_objCorrect2nd'); cd objCorrect1st_v_objCorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect1st_v_objCorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'ObjCorrect1st > ObjCorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



%% contrast BOs


% ------- compute: objCorrectAll > ObjIncorrectAll ------- %

cd(paths.group); mkdir('objCorrect_v_objIncorrect'); cd objCorrect_v_objIncorrect; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect_v_objIncorrect;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'objCorrectAll > ObjIncorrectAll';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %


% ------- compute: objCorrect1st > ObjIncorrect1st ------- %

cd(paths.group); mkdir('objCorrect1st_v_objIncorrect1st'); cd objCorrect1st_v_objIncorrect1st; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect1st_v_objIncorrect1st;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'objCorrect1st > ObjIncorrect1st';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: objCorrect1st > ObjIncorrect1st ------- %

cd(paths.group); mkdir('objCorrect2nd_v_objIncorrect2nd'); cd objCorrect2nd_v_objIncorrect2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_objCorrect2nd_v_objIncorrect2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'objCorrect2nd > ObjIncorrect2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %



% ------- compute: objCorrect1st > ObjIncorrect1st ------- %

cd(paths.group); mkdir('obj1st_v_obj2nd'); cd obj1st_v_obj2nd; dir_spm = pwd;

clear matlabbatch

spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {dir_spm};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = list_obj1st_v_obj2nd;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Obj1st > Obj2nd';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch) % run batch

% ------------------------------------------------ %

