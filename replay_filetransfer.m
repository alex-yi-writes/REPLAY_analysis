%% move files for nick

paths           = [];
paths.parent    = 'D:\work_backups\replay_dat\';
paths.mri       = [paths.parent 'mri\'];
paths.behav     = [paths.parent 'behav\'];
paths.physio    = [paths.parent 'physio\'];
paths.app       = '\Applications\MRIcron.app\Contents\Resources\dcm2niix';
paths.coreg     = [paths.parent 'coreg\'];
paths.analysis  = [paths.parent 'mri\'];
paths.group    = [paths.parent '2nd\'];

%%

for id=39:length(IDs)
    
    mkdir(['D:\work_backups\replay_dat\transmat\' IDs{id}])
    %     for ctr=1:44
    %         copyfile([paths.analysis IDs{id} '\analysis\' strcat('con_',sprintf('%04d',ctr),'.nii')],...
    %             ['D:\work_backups\replay_dat\transmat\' IDs{id} '\' strcat('con_',sprintf('%04d',ctr),'.nii')])
    %     end
    copyfile([paths.mri IDs{id} '\preproc\meanafunc.nii'],...
        ['D:\work_backups\replay_dat\transmat\' IDs{id} '\meanEPI.nii'])
    copyfile([paths.mri IDs{id} '\preproc\mt1.nii'],...
        ['D:\work_backups\replay_dat\transmat\' IDs{id} '\T1WB.nii'])
    try
        clear tmp
        tmp=dir([paths.mri IDs{id} '\converted\hackathon*.nii']);
        
        for f1=1:length(tmp)
            copyfile([paths.mri IDs{id} '\converted\' tmp(f1).name],...
                ['D:\work_backups\replay_dat\transmat\' IDs{id} '\' tmp(f1).name])
        end  
    catch
        warning(IDs{id})
    end

    try
        copyfile([paths.mri IDs{id} '\converted\maastricht.nii' ],...
            ['D:\work_backups\replay_dat\transmat\' IDs{id} '\maastricht.nii'])
        
    catch
        warning(IDs{id})
    end

    
end