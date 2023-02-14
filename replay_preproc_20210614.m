%% convert 7T MRI data

clc;clear

% paths
paths       = [];
paths.parent= 'D:\work_backups\replay_dat\';
paths.mri   = 'D:\work_backups\replay_dat\mri\';
paths.behav = 'D:\work_backups\replay_dat\behav\';
paths.physio= 'D:\work_backups\replay_dat\physio\';
paths.app   = 'C:\Users\direk\Documents\MRIcron\Resources\dcm2niix';

% experimental info
% IDs = {'102'; '103'; '104'; '105'; '106'; '107'; '108'; '109'; '110'; '111'; '112'};
% IDs = {'129';'130';'131';'132';'133';'134';'135';'136';'137'};
IDs = {'146';'147';'148';'149';'150';'151';'152'};
TR  = 2.66;
slc = 80;

%% convert dcm -> niftii

% t1 t2 epi
for id=1:length(IDs)
    
    
    mkdir([paths.mri IDs{id} '\converted\'])
    
    if str2num(IDs{id}) < 128
        % sort epi
        for r1=1:4
            clear tmp epi
            tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*ep2d_1iso_*_run' num2str(r1) ]);
            
            epi = fullfile(tmp(2).folder,tmp(2).name)
            eval(['!' paths.app ' -f run' num2str(r1) ' -p y -z n -o "' paths.mri IDs{id} '\converted" "' epi '"'])
            
        end
    elseif str2num(IDs{id}) >= 128 && str2num(IDs{id}) < 146
        % sort epi
        for r1=1:4
            clear tmp1 tmp2 epi
            tmp1 = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*ep2d_1iso_*_run' num2str(r1) ]);
            tmp2 = [num2str(str2num(tmp1(1).name(1:2))+1) '-MoCoSeries_DiCo'];
            
            epi = fullfile(tmp1(1).folder,tmp2)
            eval(['!' paths.app ' -f run' num2str(r1) ' -p y -z n -o "' paths.mri IDs{id} '\converted" "' epi '"'])
            
        end
        
    elseif str2num(IDs{id}) >= 146
        
        % sort epi
        for r1=1:4
            clear tmp1 tmp2 epi
            tmp1 = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*ep2d_bold_*_run' num2str(r1) ]);
            tmp2 = [num2str(str2num(tmp1(1).name(1:2))+1) '-MoCoSeries_DiCo'];
            
            epi = fullfile(tmp1(1).folder,tmp2)
            eval(['!' paths.app ' -f run' num2str(r1) ' -p y -z n -o "' paths.mri IDs{id} '\converted" "' epi '"'])
            
        end
        
    end
    
    % whole-brain t1
    clear tmp t1
    try
        tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*t1_mprage_*']);
        t1 = fullfile(tmp.folder,tmp.name);
        eval(['!' paths.app ' -f t1 -p y -z n -o "' paths.mri IDs{id} '\converted" "' t1 '"'])
    catch
        if str2num(IDs{id}) == 147
            tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*mp2rage*RR_UNI*']);
            t1 = fullfile(tmp(end).folder,tmp(end).name);
            eval(['!' paths.app ' -f t1 -p y -z n -o "' paths.mri IDs{id} '\converted" "' t1 '"'])
            
        elseif str2num(IDs{id}) >= 149
            tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*t1_mprage_*iso']);
            t1 = fullfile(tmp(end).folder,tmp(end).name);
            eval(['!' paths.app ' -f t1 -p y -z n -o "' paths.mri IDs{id} '\converted" "' t1 '"'])
            
        end
    end
    
    % t2 hippocampus
    clear tmp t2
    tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*perphippo*0']);
    t2 = fullfile(tmp(2).folder,tmp(2).name);
    eval(['!' paths.app ' -f t2 -p y -z n -o "' paths.mri IDs{id} '\converted" "' t2 '"'])
    
end


% lc scans 

for id=1:length(IDs)
    if  str2num(IDs{id}) < 146
        % maastricht
        clear tmp t1
        tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*multiMTC*']);
        maastricht = fullfile(tmp(1).folder,tmp(1).name);
        eval(['!' paths.app ' -f maastricht -p y -z n -o "' paths.mri IDs{id} '\converted" "' maastricht '"'])
        
        % hackathon
        clear tmp t2
        tmp = dir([paths.mri IDs{id} '\raw\Studies_DZNE\*mtflash3d*']);
        hackathon = fullfile(tmp(1).folder,tmp(1).name);
        eval(['!' paths.app ' -f hackathon -p y -z n -o "' paths.mri IDs{id} '\converted" "' hackathon '"'])
    else
        warning('skipped')
    end
end


% copy files for preprocessing

for id=6:length(IDs)
    
    mkdir([paths.mri IDs{id} '\preproc\'])
    
    copyfile([paths.mri IDs{id} '\converted\t1.nii'],...
        [paths.mri IDs{id} '\preproc\t1.nii'])
    for r1=1:4
    copyfile([paths.mri IDs{id} '\converted\run' num2str(r1) '.nii'],...
            [paths.mri IDs{id} '\preproc\run' num2str(r1) '.nii'])
    end
    
end

%% now preproc

spm fmri

% load slice timing info
preUpdate=load([paths.mri 'SliceTiming1.mat']);
postUpdate=load([paths.mri 'SliceTiming2.mat']);

for id=4:length(IDs)
    
    disp('#####################')
    disp(['#####    ' IDs{id} '    #####'])
    disp('#####################')
    
    for r1=1:4
        clear volnum list_scan
        volnum = length(spm_vol([paths.mri IDs{id} '\preproc\run' num2str(r1) '.nii']));
        for v1 = 1:volnum
            list_scan{v1,1} = [paths.mri IDs{id} '\preproc\run' num2str(r1) '.nii,' num2str(v1)];
        end
        
        
        % slice time correction
        clear sliceTcorrect
        spm_jobman('initcfg')
        sliceTcorrect{1}.spm.temporal.st.scans = {list_scan};
        sliceTcorrect{1}.spm.temporal.st.nslices = slc;
        sliceTcorrect{1}.spm.temporal.st.tr = TR;
        sliceTcorrect{1}.spm.temporal.st.ta = TR-(TR\slc);
        if str2num(IDs{id}) < 146
        sliceTcorrect{1}.spm.temporal.st.so = preUpdate.SliceTiming;
        else
        sliceTcorrect{1}.spm.temporal.st.so = postUpdate.SliceTiming;
        end
        sliceTcorrect{1}.spm.temporal.st.refslice = 1;
        sliceTcorrect{1}.spm.temporal.st.prefix = 'a';
        spm_jobman('run',sliceTcorrect)
        
    end
    
    disp('slice-time corrected')
    
    % realign
    clear matlabbatch
    spm_jobman('initcfg')
    for r1=1:4
    matlabbatch{r1}.spm.util.split.vol = {[paths.mri IDs{id} '\preproc\arun' num2str(r1) '.nii']};
    matlabbatch{r1}.spm.util.split.outdir = {[paths.mri IDs{id} '\preproc\']};
    end
    spm_jobman('run',matlabbatch)
    
    clear files3D
    files3D = dir([paths.mri IDs{id} '\preproc\arun*_00*.nii']); files3D={files3D.name}';
    cd([paths.mri IDs{id} '\preproc\'])
    clear matlabbatch
    spm_jobman('initcfg')
    matlabbatch{1}.spm.util.cat.vols = files3D;
    matlabbatch{1}.spm.util.cat.name = 'afunc.nii';
    matlabbatch{1}.spm.util.cat.dtype = 4;
    matlabbatch{1}.spm.util.cat.RT = TR;
    spm_jobman('run',matlabbatch)
    
    for r1=1:4
    eval(['!del ' paths.mri IDs{id} '\preproc\arun' num2str(r1) '_00*.nii'])
    end
    
    clear volnum list_scan
    volnum = length(spm_vol([paths.mri IDs{id} '\preproc\afunc.nii']));
    for v1 = 1:volnum
        list_scan{v1,1} = [paths.mri IDs{id} '\preproc\afunc.nii,' num2str(v1)];
    end
    
    clear realigned
    spm_jobman('initcfg')
    realigned{1}.spm.spatial.realign.estwrite.data = {list_scan}';
    realigned{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realigned{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realigned{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    realigned{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    realigned{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
    realigned{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realigned{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    realigned{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    realigned{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realigned{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realigned{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realigned{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',realigned)
    
    disp('realigned')
    
    
    % smoothe
    clear volnum list_scan
    volnum = length(spm_vol([paths.mri IDs{id} '\preproc\rafunc.nii']));
    for v1 = 1:volnum
        list_scan{v1,1} = [paths.mri IDs{id} '\preproc\rafunc.nii,' num2str(v1)];
    end
    clear smoothed
    spm_jobman('initcfg')
    smoothed{1}.spm.spatial.smooth.data = list_scan;
    smoothed{1}.spm.spatial.smooth.fwhm = [1.5 1.5 1.5];
    smoothed{1}.spm.spatial.smooth.dtype = 0;
    smoothed{1}.spm.spatial.smooth.im = 0;
    smoothed{1}.spm.spatial.smooth.prefix = 's15';
    spm_jobman('run',smoothed)
    
    disp('smoothed')
    
    
    for r1=1:4
    eval(['!del ' paths.mri IDs{id} '\preproc\arun' num2str(r1) '.nii'])
    eval(['!del ' paths.mri IDs{id} '\preproc\run' num2str(r1) '.nii'])
    end
    eval(['!del ' paths.mri IDs{id} '\preproc\afunc.nii'])
    
    disp('cleaned up')
    
    
    % bias correct t1
    clear biascorrect
    spm_jobman('initcfg')
    biascorrect{1}.spm.spatial.preproc.channel.vols = {[paths.mri IDs{id} '\preproc\t1.nii,1']};
    biascorrect{1}.spm.spatial.preproc.channel.biasreg = 0.01;
    biascorrect{1}.spm.spatial.preproc.channel.biasfwhm = 30;
    biascorrect{1}.spm.spatial.preproc.channel.write = [0 1];
    biascorrect{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\direk\Documents\MATLAB\spm12\tpm\TPM.nii,1'};
    biascorrect{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    biascorrect{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    biascorrect{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    biascorrect{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\direk\Documents\MATLAB\spm12\tpm\TPM.nii,2'};
    biascorrect{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    biascorrect{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    biascorrect{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    biascorrect{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\direk\Documents\MATLAB\spm12\tpm\TPM.nii,3'};
    biascorrect{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    biascorrect{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    biascorrect{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    biascorrect{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\direk\Documents\MATLAB\spm12\tpm\TPM.nii,4'};
    biascorrect{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    biascorrect{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    biascorrect{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    biascorrect{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\direk\Documents\MATLAB\spm12\tpm\TPM.nii,5'};
    biascorrect{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    biascorrect{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    biascorrect{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    biascorrect{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\direk\Documents\MATLAB\spm12\tpm\TPM.nii,6'};
    biascorrect{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    biascorrect{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    biascorrect{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    biascorrect{1}.spm.spatial.preproc.warp.mrf = 1;
    biascorrect{1}.spm.spatial.preproc.warp.cleanup = 1;
    biascorrect{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    biascorrect{1}.spm.spatial.preproc.warp.affreg = 'mni';
    biascorrect{1}.spm.spatial.preproc.warp.fwhm = 0;
    biascorrect{1}.spm.spatial.preproc.warp.samp = 3;
    biascorrect{1}.spm.spatial.preproc.warp.write = [0 0];
    biascorrect{1}.spm.spatial.preproc.warp.vox = NaN;
    biascorrect{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
        NaN NaN NaN];
    
    spm_jobman('run',biascorrect)
    
    eval(['!rm ' paths.mri IDs{id} '\preproc\c*t1.nii'])
    
    disp('done')
    
end
