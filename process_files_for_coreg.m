IDs = {'123'};


for id=1:length(IDs)
    mkdir(['/Users/alex/Documents/replay_dat/coreg/' IDs{id} '/data/'])
    
    
    copyfile(['/Users/alex/Documents/replay_dat/mri/' IDs{id} '/preproc/mt1.nii'],...
        ['/Users/alex/Documents/replay_dat/coreg/' IDs{id} '/data/T1WB_corrected.nii'])
    copyfile(['/Users/alex/Documents/replay_dat/mri/' IDs{id} '/preproc/meanafunc.nii'],...
        ['/Users/alex/Documents/replay_dat/coreg/' IDs{id} '/data/meanEPI.nii'])
    
    
    
end