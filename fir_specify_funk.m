function [b] = fir_specify_funk(b)

% Block lss one subject and run at a time. Creates one spm.mat for each
% sequence (5 consecutive trials) and across 3 repetitions (3 chunks of 5
% trials). One Beta image for all 15 of those trials. After estimation
% beta images will be passed to RSA toolbox. 
% Author: Walter Reilly

% Usage:
%
%	b = lss_block(b)
%   
%   input arguments:
%
%	b = memolab qa batch structure containing the fields:
%
%       b.dataDir     = fullpath string to the directory where the functional MRI data
%                       is being stored
%
%       b.runs        = cellstring with IDs for each functional time series
%
%       b.auto_accept = a true/false variable denoting whether or not the 
%                       user wants to be prompted
%
%       b.rundir      = a 1 x n structure array, where n is the number of
%                       runs, with the fields:
%
%       b.rundir.rp = motion regressors
% 
%       b.rundir.smfiles  = smoothed volumes

%% get motion regressor rp file names 
% remember to change batch to specify rp if change back to rp from spike
% reg e
for i = 1:length(b.runs)
    b.rundir(i).rp     = spm_select('FPListRec', b.dataDir, ['^rp.*' b.runs{i} '.*bold\.txt']);
    fprintf('%02d:   %s\n', i, b.rundir(i).rp)
end

%% get spike regressor file names
% for i = 1:length(b.runs)
%     b.rundir(i).spike     = spm_select('FPListRec', b.dataDir, ['^spike.*' b.runs{i} '.*\.txt']);
%     fprintf('%02d:   %s\n', i, b.rundir(i).spike)
% end

%% get smoothed nii names
for i = 1:length(b.runs)
    % print success
    b.rundir(i).smfiles = spm_select('ExtFPListRec', b.dataDir, ['^coreg_smoothed.*'  b.runs{i} '.*bold\.nii']);
    fprintf('%02d:   %0.0f smoothed files found.\n', i, length(b.rundir(i).smfiles))
end % end i b.runs


%% get condition files from saved .mat
% cond_dir = '/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_cond_files';
cond_dir   = '/home/wbreilly/sms_scan_crick/cluster_fir_data_10_27_18/fir_cond_files/';

for i = 1:length(b.runs)
    b.rundir(i).cond = cellstr(spm_select('FPList', cond_dir, [ 'fir_condfile_' b.curSubj sprintf('_%s.mat', b.runs{i})]));
end % end i b.runs
        
%%
% hygiene is important
clear matlabbatch
        
fir_dir_spm = fullfile(b.dataDir, 'fir_spm');
mkdir(char(fir_dir_spm));

%% for sessions
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fir_dir_spm);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.22;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 38;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
 
%% run level stuff
for irun = 1:length(b.runs)
    cond_file = b.rundir(irun).cond;
    matlabbatch{1}.spm.stats.fmri_spec.sess(irun).scans = cellstr(b.rundir(irun).smfiles);
    matlabbatch{1}.spm.stats.fmri_spec.sess(irun).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(irun).multi = cond_file;
    matlabbatch{1}.spm.stats.fmri_spec.sess(irun).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(irun).multi_reg = cellstr(b.rundir(irun).rp);
    matlabbatch{1}.spm.stats.fmri_spec.sess(irun).hpf = 128;
end %irun
%% also for all sessions
    
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
% double check this
matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = 39.04;
% seconds or TRs
matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order = 32;
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = .1;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
 %%      

% run
spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);



%% simple enough to estimate here as well
clear matlabbatch



fir_dir_spm = fullfile(b.dataDir, 'fir_spm', 'SPM.mat');

matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(fir_dir_spm);
% changing this to 1 keeps res for every TR
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
      
    
    

end % end function
