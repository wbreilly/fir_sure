% copy and rename betas to new folder in a way that unifies names across
% subjects for compatibility with rsa toolbox
% Walter Reilly
% created 8_29_18


clear all
clc


dataDir     = '/Users/wbr/walter/fmri/sms_scan_analyses/data_for_spm/fir_data_10_26_18/';
scriptdir   = '/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_sure'; 


% the order of the names for each run in the cond file defined the order of
% the betas in the glm
fir_cond_file_path = '/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_cond_files/';

subjects    = {'s001' 's002' 's003' 's004' 's007' 's008' 's009' 's010' 's011' 's015' 's016' 's018' 's019'  's020'...
               's022' 's023' 's024' 's025' 's027' 's028' 's029' 's030' 's032' 's033' 's034' 's035' 's036' 's037' ...
               's038' 's039' 's040' 's041' 's042' 's043'};
           
% 32 TRs/regressors each for 5 sequences per run, 9 runs, 6 motion reg at
% the end of each run. final unused regressors are 9 constant, 1 for each
% run

%% these indices are the first regressor index of each sequence in each run
first_beta_idx = [1:32:160;...
                  167:32:326;...
                  333:32:492;...
                  499:32:658;...
                  665:32:824;...
                  831:32:990;...
                  997:32:1156;...
                  1163:32:1322;...
                  1329:32:1488];

last_beta_idx  = [32:32:160;...
                  198:32:326;...
                  364:32:492;...
                  530:32:658;...
                  696:32:824;...
                  862:32:990;...
                  1028:32:1156;...
                  1194:32:1322;...
                  1360:32:1488];
%%           
% for icheck = 1:45
%     size(first_beta_idx(icheck):last_beta_idx(icheck),2)
% end


tic
parfor isub = 8:length(subjects) 
    
    b.curSubj   = subjects{isub};
    b.dataDir   = fullfile(dataDir, b.curSubj, 'fir_spm');
    b.rename_dir = fullfile(dataDir, b.curSubj, 'rename_fir_spm');
    mkdir(char(b.rename_dir));
    b.fir_cond_file_path = fir_cond_file_path;
    
    % funk it up
    [b] = rename_beta_funk(b);
%     
%     % grab all the betas 
%     b.allbetas = spm_select('ExtFPListRec', b.dataDir, 'beta_.*.nii');
%     
%     all_names = cell(9,5);
%     for irun = 1:9
%         load(sprintf('%sfir_condfile_%s_Rifa_%d.mat',b.fir_cond_file_path,subjects{isub},irun),'names');
%         all_names(irun,:) = names;
%     end
%     
%     % all unique condition labels
%     all_unique_seq = unique(all_names);
%     
%     for iuniq = 1:15
%         
%         cur_seq = all_unique_seq{iuniq};
%         
%         % make a mask for when current sequence occurred in each run
%         cur_mask = logical(ismember(all_names, cur_seq));
% 
%         % mask the start and end beta indices to know what betas to grab
%         all_start_idxs = first_beta_idx(cur_mask);
%         all_end_idxs = last_beta_idx(cur_mask);
%         
%         % also want to know which run the sequence occured
%         % runs column-wise and so does the masking above so that's nice
%         [seq_runs,notneeded] = find(cur_mask);
%         
%         %% let's do this
%         for irow = 1:3
%             start_idx = all_start_idxs(irow);
%             last_idx  = all_end_idxs(irow);
%             cur_run   = seq_runs(irow);
%             
%             % initiate time counter
%             time_num = 0;
%             
%             % loop through all 32 timepoints 
%             for ibeta = start_idx:last_idx
%                 
%                 % keep track of within sequence time
%                 time_num = time_num + 1;
%                 % source beta
%                 [path,name,ext] = fileparts(b.allbetas(ibeta,:));
%                 % dest name 
%                 dest_name = sprintf('%s_t%02d_run%d',cur_seq,time_num,cur_run);
%                 
%                 % copy with new name to new dir
%                 copyfile(fullfile(b.dataDir, [name, '.nii']), fullfile(b.rename_dir,[dest_name, '.nii']));
%             
%             end % ibeta 32 betas for the FIR window
%         end % irow the 3 runs in which each sequence appears     
%     end % iuniq sequences
%     
%     %% now unify names across subs
%     
%     b.unifynames = spm_select('ExtFPListRec', b.rename_dir, '.*.nii');
%     
%     icur = 1;
%     
%     % VERY IMPORTANT that file names with numbers ranging from singles to
%     % tens place have same length files names i.e. 01 02 using %02d. Mac
%     % sorts these correctly but spm_select does not!!
%     
%     while icur < 577
%         for iseq = 1:6
%             for ipos = 1:32
%                 for irun = 1:3
%                     [path,oldname,ext] = fileparts(b.unifynames(icur,:));
%                     newname = sprintf('intact_seq%d_pos%02d_run%d.nii',iseq,ipos,irun);
%                     movefile(fullfile(path, [oldname, '.nii']), fullfile(path, newname));
%                     icur = icur + 1;
%                 end % end irun
%             end % end ipos
%         end % end iseq     
%     end % end while intact
%     
%     while icur < 865
%         for iseq = 1:3
%             for ipos = 1:32
%                 for irun = 1:3
%                     [path,oldname,ext] = fileparts(b.unifynames(icur,:));
%                     newname = sprintf('random_seq%d_pos%02d_run%d.nii',iseq,ipos,irun);
%                     movefile(fullfile(path, [oldname, '.nii']), fullfile(path, newname));
%                     icur = icur + 1;
%                 end % end irun
%             end % end ipos
%         end % end iseq     
%     end % end while random
%     
%     while icur < 1440
%         for iseq = 1:6
%             for ipos = 1:32
%                 for irun = 1:3
%                     [path,oldname,ext] = fileparts(b.unifynames(icur,:));
%                     newname = sprintf('scrambled_seq%d_pos%02d_run%d.nii',iseq,ipos,irun);
%                     movefile(fullfile(path, [oldname, '.nii']), fullfile(path, newname));
%                     icur = icur + 1;
%                 end % end irun
%             end % end ipos
%         end % end iseq     
%     end % end while scrambled
%     
    
end %isub
toc









