% gettrialinfo_sms_scan
% Walter Reilly
% created 8_22_17 for sms_scan to extract ad create condition names, onsets,
% and durations

clear all
clc

path = '/Volumes/GoogleDrive/My Drive/grad_school/DML_WBR/Sequences_Exp3/sms_scan_drive/sms_scan_fmri_copy/';

subjects    = {'s001' 's002' 's003' 's004' 's007' 's008' 's009' 's010' 's011' 's015' 's016' 's018' 's019'  's020'...
               's022' 's023' 's024' 's025' 's027' 's028' 's029' 's030' 's032' 's033' 's034' 's035' 's036' 's037' ...
               's038' 's039' 's040' 's041' 's042' 's043'};
           
% can't believe I didn't do this sooner..
subs = strrep(subjects, 's00','');
subs = strrep(subs, 's0','');
sub2 = strrep(subjects, 's0','');

for isub = 1:length(subs) 
    for irrb = 1:3
        for iblock = 1:3
            
            load(sprintf('%ss%s_rrb%d_%d.mat',path,sub2{isub},irrb,iblock));

            % cell array of condition names
            % 5 sequences per run, deal with repetitions below
            names = {};
            for iname = 1:5
                names{1,iname} = sprintf('%s_%02d',RETRIEVAL(1).sequence{iname,3},RETRIEVAL(1).sequence{iname,4});
            end
                
            % duration 
            durations{1} = 0;
            durations{2} = 0;
            durations{3} = 0;
            durations{4} = 0;
            durations{5} = 0;
            
            % all the sequences
            nreps = 3;
            allrunseq = {};
            for irep = 1:nreps
                for iname = 1:5
                    allrunseq{irep,iname} = sprintf('%s_%02d',RETRIEVAL(irep).sequence{iname,3},RETRIEVAL(irep).sequence{iname,4});
                end
            end % end irep
            
            % order of first rep is given
            seq1 = [1];
            seq2 = [2];
            seq3 = [3];
            seq4 = [4];
            seq5 = [5];
            
            % now find when those sequences occur in second and third
            % sequences
            for iseq = 1:5
                for irep = 2:3
                    for iname = 1:5
                        if strcmp(allrunseq{1,iseq},allrunseq{irep,iname})
                            eval(sprintf('seq%d = [seq%d; %d];',iseq,iseq,iname));
                        end
                    end
                end
            end

            % convert the idx's of sequence order into onsets for each condition
            % 7 dummy TR's at beginning so everything is +
            onsets{1} = [];
            onsets{2} = [];
            onsets{3} = [];
            onsets{4} = [];
            onsets{5} = [];
            
            % for second and thirs reps, onset is after 25 TR's and 5
            % sequences already 
            reppad = 25*5;

            % write onsets for each condition (aka seq)
            for iseq = 1:5
                for irep = 1:3
                    if iseq == 1 && irep == 1
                        onsets{1} = 8;
                    else
                        eval(sprintf('onsets{%d} = [onsets{%d} (seq%d(%d)-1)*25+8+%d*(%d-1)];',iseq,iseq,iseq,irep,reppad,irep));
                    end % end if  
                end
            end % end icond  
            
            %where to save condition files
            savepath = '/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_cond_files/';
            
            % save with run naming convention (1-9)
            % dumb way
            
            if irrb == 1
                if iblock == 1;
                    run = 1;
                elseif iblock == 2
                    run = 2;
                else 
                    run = 3;
                end
            elseif irrb == 2
                if iblock == 1;
                    run = 4;
                elseif iblock == 2
                    run = 5;
                else 
                    run = 6;
                end
            else
                if iblock == 1;
                    run = 7;
                elseif iblock == 2
                    run = 8;
                else 
                    run = 9;
                end
            end
            
            % save
            save(sprintf('%sfir_condfile_%s_Rifa_%d.mat',savepath,subjects{isub},run),'names', 'durations', 'onsets');
            
            clearvars -EXCEPT isub irrb iblock path subjects sub2
            
        end % end iblock
    end % end irrb
end % end isub




% part of the old way I did this in univariate get trial info
% onsets_tmp(1,icond) = (intact_idx_tmp(icond)-1) *25 + 8;









