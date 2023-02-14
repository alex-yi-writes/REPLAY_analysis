%% generate physiological noise matrices for SPM analyses, for the 7T data

%% work log

%   28-02-2020      generated the script

%% set env

clear;clc;warning('off')

% set up search paths
paths = [];
paths.parent    = 'D:\work_backups\replay_dat\';
paths.physio    = [paths.parent 'physio\'];
paths.task      = [paths.parent 'mri\'];
paths.save      = [paths.parent 'physio\'];

% define fixed values
TR          = 2.66;
samplingRate= 1/200; % 200 Hz
numSlices   = 80;
SOI         = 1;

% subjects
% IDs = [102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112];
IDs = {'146';'147';'148'};

for id=1:length(IDs)
    mkdir([paths.physio IDs{id} '\raw\'])
end

%% read the raw data

for id = 2:length(IDs)
    
    disp(IDs{id})
    
    %% extract
    mkdir([paths.physio IDs{id} '\made\'])
    cd([paths.physio IDs{id} '\made\'])
    
    for runs=1:4
        
        if str2num(IDs{id})==102 && runs==1
            
            disp('no physio data for id 102 run 1. making dummy data')
            
        else
            close all
            clear filepath filepath physio numCols tLines breathing_t pulse_t breathing_raw pulse_raw triggers...
                fName_br fName_pul scannertrig filename
            
            try
            filename = dir([paths.physio IDs{id} '\raw\*200Hz_*' num2str(runs) '.txt']); filepath = fullfile(filename.folder,filename.name);
            catch
            filename = dir([paths.physio IDs{id} '\raw\*\*200Hz_*' num2str(runs) '.txt']); filepath = fullfile(filename.folder,filename.name);
            end
            
            delimiter = ' '; % our physio file is delimited by spaces
            fid = fopen(filepath,'rt');
            tLines = fgets(fid);
            numCols = numel(strfind(tLines,delimiter))+1;
            numRows = 4;
            fclose(fid);
            
            [physio]     = importPhysioFile(filepath,numCols);
            
            breathing_raw= physio(1,:);
            pulse_raw    = physio(2,:);
            triggers     = physio(4,:);
            
            scannertrig  = triggers(:)>=triggers(1,1);
            
            % transform files and save
            %     fName = fullfile(filename.folder,[filename.name(1:end-4) '_t.txt']);
            %     physio_BL_t = transpose(physio_BL);
            %     dlmwrite(fName, physio_BL_t, 'delimiter',' ','newline','pc','precision',13);
            fName_br = fullfile(filename.folder,[IDs{id} '_breathing_' num2str(runs) '.txt']);
            breathing_bl_t = transpose(breathing_raw);
            dlmwrite(fName_br, breathing_bl_t, 'delimiter',' ','newline','pc','precision',13);
            fName_pul= fullfile(filename.folder,[IDs{id} '_cardiac_' num2str(runs) '.txt']);
            pulse_bl_t = transpose(pulse_raw);
            dlmwrite(fName_pul, pulse_bl_t, 'delimiter',' ','newline','pc','precision',13);
            fName_trig= fullfile(filename.folder,[IDs{id} '_triggers_' num2str(runs) '.txt']);
            trigg_bl_t = transpose(triggers);
            dlmwrite(fName_trig, trigg_bl_t, 'delimiter',' ','newline','pc','precision',13);
            
            % prepare parameters for the toolbox
            clear nvols
            nvols = length(spm_vol([paths.task IDs{id} '\converted\run' num2str(runs) '.nii']))
            
            % run TAPAS PhysIO
            clear matlabbatch
            spm_jobman('initcfg');      % initiate job manager
            matlabbatch{1}.spm.tools.physio.save_dir = {[paths.save IDs{id} '\made\']};
            matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Custom';
            matlabbatch{1}.spm.tools.physio.log_files.cardiac = {fName_pul};
            matlabbatch{1}.spm.tools.physio.log_files.respiration = {fName_br};
            matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {};%{fName_trig};
            matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = samplingRate;
            matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
            matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'first';
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = numSlices;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = nvols;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = SOI;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.filter.no = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.respiratory.filter.passband = [0.01 2];
            matlabbatch{1}.spm.tools.physio.preproc.respiratory.despike = false;
            matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = [IDs{id} '_physio' num2str(runs) '.txt'];
            matlabbatch{1}.spm.tools.physio.model.output_physio = [IDs{id} '_physio' num2str(runs) '.mat'];
            matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
            matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
            matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
            matlabbatch{1}.spm.tools.physio.verbose.level = 2;
            matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = '';
            matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;
            spm_jobman('run', matlabbatch) % run batch
            
        end
    end
end