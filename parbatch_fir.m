% Author: Walter Reilly 
% FIR model for sms_scan
% Created 10_26_18
% made parallel 10_20_18 not that it took much


%====================================================================================
%			Specify Variables
%====================================================================================

%-- Directory Information
% Paths to relevant directories.
% dataDir   = path to the directory that houses the MRI data
% scriptdir = path to directory housing this script (and auxiliary scripts)
% QAdir     = Name of output QA directory

% dataDir     = '/Users/wbr/walter/fmri/sms_scan_analyses/data_for_spm/fir_data_10_26_18/';
% scriptdir   = '/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_sure'; 
dataDir       = '/home/wbreilly/sms_scan_crick/cluster_fir_data_10_27_18/';
scriptdir     = '/home/wbreilly/sms_scan_crick/cluster_fir_data_10_27_18/fir_sure';

% add spm to path
% this version hasn't been compiled yet
% addpath /home/wbreilly/matlab_toolboxes/spm12/
% this is a version jordan compiled. He edited something to do with
% implicit masking and a .mat default save format
addpath /group/dml/apps/spm12


subjects    = {'s001' 's002' 's003' 's004' 's007' 's008' 's009' 's010' 's011' 's015' 's016' 's018' 's019'  's020'...
               's022' 's023' 's024' 's025' 's027' 's028' 's029' 's030' 's032' 's033' 's034' 's035' 's036' 's037' ...
               's038' 's039' 's040' 's041' 's042' 's043'};
runs        = {'Rifa_1' 'Rifa_2' 'Rifa_3' 'Rifa_4' 'Rifa_5' 'Rifa_6' 'Rifa_7' 'Rifa_8' 'Rifa_9'};  

%-- Auto-accept
% Do you want to run all the way through without asking for user input?
% if 0: will prompt you to take action;
% if 1: skips realignment and ArtRepair if already run, overwrites output files

auto_accept = 0;

fprintf('Initializing and checking paths.\n')


%-- Check for required functions

% SPM
if exist('spm','file') == 0
    error('SPM must be on the path.')
end

fprintf('Running batch lss')

% pool party
pc = parcluster('big_mem'); % or 'single_nose'
poolobj = parpool(pc, 2);

%--Loop over subjects
parfor i = 1:length(subjects)
%     try
        % Define variables for individual subjects - General
        b.curSubj   = subjects{i};
        b.runs      = runs;
        b.dataDir   = fullfile(dataDir, b.curSubj);

        % Define variables for individual subjects - QA General
        b.scriptdir   = scriptdir;
        b.auto_accept = auto_accept;

        % Check whether first level has already been run for a subject

        % Initialize diary for saving output
        diaryname = fullfile(b.dataDir, 'batch_fir_spec_and_est_diary.txt');
        diary(diaryname);

        %======================================================================
        % Run functions (at this point, this could all be in one
        % script/function, but where's the fun in that?
        %======================================================================

        % Run lss script
        fprintf('--FIRing SqUaD--\n')
        [b] = fir_specify_funk(b);
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')

%     catch ME
%         disp(['ID: ' ME.identifier])
%         warning('bad shit happened for %s', subjects{i})
%     end % end catch
    diary off
end % i (subjects)

fprintf('FIR FTW!!\n')


delete(gcp('nocreate'))
exit
