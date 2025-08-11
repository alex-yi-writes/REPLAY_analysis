%% average LC images

clc;clear;

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin']);
setenv('ANTSPATH','/usr/local/bin')

paths=[];
paths.original='/Users/yeojin/Desktop/E_data/EA_raw/EAB_MRI/EABY_originals/REPLAY/';
paths.converted='/Users/yeojin/Desktop/E_data/EA_raw/EAB_MRI/EABB_preprocessed/REPLAY/';

IDs = {'103'; '104'; '105'; '106'; '107'; '108'; '109'; '110';...
    '111'; '112'; '114'; '115'; '119'; '120'; '122'; '123'; '124'; '125'; '126'; '127'};


%% convert images, use dcm2nii


for id=13:length(IDs)
    
    eval(['!mkdir ' paths.converted IDs{id}])
    clear tmppath tmp tmp2
    
    tmp=dir([paths.original IDs{id}]);
    tmp=tmp(~cellfun(@(x) x==0, {tmp.isdir}));
    tmp=tmp(~ismember({tmp.name} ,{'.','..','.DS_Store'}));
    tmp2=dir([tmp(1).folder '/' tmp(1).name '/Studies_*']);
    tmppath=[tmp2(1).folder '/' tmp2(1).name];
    
    clear tmp3
    tmp3=dir([tmppath '/*mtflash3d*']);
    clear MTwpath
    MTwpath=[tmp3(1).folder '/' tmp3(1).name];
    % convert
    eval(['!/Applications/MRIcron.app/Contents/Resources/dcm2niix -f "hackathon" -p y -z y -ba n "' MTwpath '"'])
    
    % move
    for series=1:4
        eval(['!mv -v ' MTwpath '/hackathon_e' num2str(series) '.nii.gz ' paths.converted IDs{id} '/hackathon' num2str(series) '.nii.gz'])
        eval(['!mv -v ' MTwpath '/hackathon_e' num2str(series) '.json ' paths.converted IDs{id} '/hackathon' num2str(series) '.json'])
    end
    
end


for id=13:length(IDs)
    
    clear tmppath tmp tmp2
    
    tmp=dir([paths.original IDs{id}]);
    tmp=tmp(~cellfun(@(x) x==0, {tmp.isdir}));
    tmp=tmp(~ismember({tmp.name} ,{'.','..','.DS_Store'}));
    tmp2=dir([tmp(1).folder '/' tmp(1).name '/Studies_DZNE*']);
    tmppath=[tmp2(1).folder '/' tmp2(1).name];
    
    clear tmp3
    tmp3=dir([tmppath '/*multiMTC*']);
    clear MTwpath
    MTwpath=[tmp3(1).folder '/' tmp3(1).name];
    % convert
    eval(['!/Applications/MRIcron.app/Contents/Resources/dcm2niix -f "maastricht" -p y -z y -ba n "' MTwpath '"'])
    
    % move
    eval(['!mv -v ' MTwpath '/maastricht.nii.gz ' paths.converted IDs{id} '/maastricht.nii.gz'])
    eval(['!mv -v ' MTwpath '/maastricht.json ' paths.converted IDs{id} '/maastricht.json'])
    
end

%% coregister and average

for id=[4:6 14:length(IDs)]
    
    referenceImage=[paths.converted IDs{id} '/hackathon1.nii.gz'];
    clear regstring
    regstring=[paths.converted IDs{id} '/hackathon1.nii.gz '];
    for s1=2:4
        eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -f ' referenceImage ' -m ' paths.converted IDs{id} '/hackathon' num2str(s1) '.nii.gz -o ' paths.converted IDs{id} '/hackathon'  num2str(s1) '_' ])
        regstring=[regstring paths.converted IDs{id} '/hackathon'  num2str(s1) '_Warped.nii.gz '];
    end
    eval(['!AverageImages 3 ' paths.converted IDs{id} '/hackathon_averaged.nii 2 ' regstring])
    eval(['!rm ' paths.converted IDs{id} '/hackathon*.nii.gz'])
    
    % bias-field correct
    clear matlabbatch
    spm_jobman('initcfg')
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {[paths.converted IDs{id} '/hackathon_averaged.nii,1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01; % medium regularisation
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30; % kernel 30mm cutoff
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Applications/spm12/tpm/TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Applications/spm12/tpm/TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Applications/spm12/tpm/TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/Applications/spm12/tpm/TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Applications/spm12/tpm/TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Applications/spm12/tpm/TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    spm_jobman('run',matlabbatch)
    
    % request cleanup
    eval(['!rm ' paths.converted IDs{id} '/c*.nii'])
    eval(['!rm ' paths.converted IDs{id} '/hackathon*.mat'])

end