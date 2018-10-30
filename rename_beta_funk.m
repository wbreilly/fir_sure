function [b] = rename_beta_funk(b)

load('beta_indices.mat')
    
    


% grab all the betas 
b.allbetas = spm_select('ExtFPListRec', b.dataDir, 'beta_.*.nii');

all_names = cell(9,5);
for irun = 1:9
    load(sprintf('%sfir_condfile_%s_Rifa_%d.mat',b.fir_cond_file_path,b.curSubj,irun),'names');
    all_names(irun,:) = names;
end

% all unique condition labels
all_unique_seq = unique(all_names);

for iuniq = 1:15

    cur_seq = all_unique_seq{iuniq};

    % make a mask for when current sequence occurred in each run
    cur_mask = logical(ismember(all_names, cur_seq));

    % mask the start and end beta indices to know what betas to grab
    all_start_idxs = first_beta_idx(cur_mask);
    all_end_idxs = last_beta_idx(cur_mask);

    % also want to know which run the sequence occured
    % runs column-wise and so does the masking above so that's nice
    [seq_runs,notneeded] = find(cur_mask);

    %% let's do this
    for irow = 1:3
        start_idx = all_start_idxs(irow);
        last_idx  = all_end_idxs(irow);
        cur_run   = seq_runs(irow);

        % initiate time counter
        time_num = 0;

        % loop through all 32 timepoints 
        for ibeta = start_idx:last_idx

            % keep track of within sequence time
            time_num = time_num + 1;
            % source beta
            [path,name,ext] = fileparts(b.allbetas(ibeta,:));
            % dest name 
            dest_name = sprintf('%s_t%02d_run%d',cur_seq,time_num,cur_run);

            % copy with new name to new dir
            copyfile(fullfile(b.dataDir, [name, '.nii']), fullfile(b.rename_dir,[dest_name, '.nii']));

        end % ibeta 32 betas for the FIR window
    end % irow the 3 runs in which each sequence appears     
end % iuniq sequences

%% now unify names across subs

b.unifynames = spm_select('ExtFPListRec', b.rename_dir, '.*.nii');

icur = 1;

% VERY IMPORTANT that file names with numbers ranging from singles to
% tens place have same length files names i.e. 01 02 using %02d. Mac
% sorts these correctly but spm_select does not!!

while icur < 577
    for iseq = 1:6
        for ipos = 1:32
            for irun = 1:3
                [path,oldname,ext] = fileparts(b.unifynames(icur,:));
                newname = sprintf('intact_seq%d_pos%02d_run%d.nii',iseq,ipos,irun);
                movefile(fullfile(path, [oldname, '.nii']), fullfile(path, newname));
                icur = icur + 1;
            end % end irun
        end % end ipos
    end % end iseq     
end % end while intact

while icur < 865
    for iseq = 1:3
        for ipos = 1:32
            for irun = 1:3
                [path,oldname,ext] = fileparts(b.unifynames(icur,:));
                newname = sprintf('random_seq%d_pos%02d_run%d.nii',iseq,ipos,irun);
                movefile(fullfile(path, [oldname, '.nii']), fullfile(path, newname));
                icur = icur + 1;
            end % end irun
        end % end ipos
    end % end iseq     
end % end while random

while icur < 1440
    for iseq = 1:6
        for ipos = 1:32
            for irun = 1:3
                [path,oldname,ext] = fileparts(b.unifynames(icur,:));
                newname = sprintf('scrambled_seq%d_pos%02d_run%d.nii',iseq,ipos,irun);
                movefile(fullfile(path, [oldname, '.nii']), fullfile(path, newname));
                icur = icur + 1;
            end % end irun
        end % end ipos
    end % end iseq     
end % end while scrambled


end % function