% write a 4d nii for showing t values from rsa toolbox > r analyses of
% sliding time window PS between intact and scrambled with glasser ROIs in
% 3 network DMN 

clear all
clc


dataDir     = '/Users/wbr/walter/fmri/rsfc_glasser_ROIs';
subjects    = {'s202','s204','s205','s206','s207','s208','s209','s210','s211',...
            's212','s213','s214','s216','s218','s219','s220','s221',...
            's223','s224','s225','s226','s227','s228','s229','s230','s231',...
            's233','s234','s235','s236','s237','s238','s240','s241','s243',...
            's244','s245','s246','s247','s248','s249','s250'};



% Loop over subjects
isub = 1; %:length(subjects)
b.curSubj   = subjects{isub};
b.dataDir   = fullfile(dataDir, b.curSubj);
b.masks     = cellstr(spm_select('ExtFPListRec', b.dataDir, 'reslice_HCP-MMP1.nii'));

% read in the lookup table
% save txt as windows formatted text in excel!!
firDir = '/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_sure';
fid = fopen(fullfile(firDir,'brainmoviedata_1_3_19.txt'));
% pval tval timepoint roi network
data = textscan(fid,'%f%f%d%s%s');
fclose(fid);



% read in the mask
temp = spm_vol(char(b.masks(1,:)));
brain = spm_read_vols(temp);

for xvox = 1:64
    for yvox = 1:64
        for zvox = 1:38
            if brain(xvox,yvox,zvox) < 999
                brain(xvox,yvox,zvox) = NaN;
            end
        end
    end
end

brain_4d = repmat(brain,[1 1 1 32]);

%intialize counter
counter = 1;
% loop through rois
for itime = 1:32
    for iroi = 1:size(threenetworkLUT,1)
        
        % the data have values starting with first timepoint for every roi,
        % then every roi for the second timepount. Rois are always in the
        % same order, so the LUT is only the length of the number of ROIs,
        % instead of the length of the data
        cur_val = data{1,2}(counter);

        % loop through brain 
        for xvox = 1:64
            for yvox = 1:64
                for zvox = 1:38
                    if brain(xvox,yvox,zvox) == threenetworkLUT(iroi)
                        brain_4d(xvox,yvox,zvox,itime) = cur_val;
                    end
                end
            end
        end
        counter = counter+1;
    end %roi
end %itime
% loop through brain

for itime = 1:32
    for xvox = 1:64
        for yvox = 1:64
            for zvox = 1:38
                if brain_4d(xvox,yvox,zvox,itime) > 999
                    brain_4d(xvox,yvox,zvox,itime) = NaN;
                end
            end
        end
    end
end


% have to delete pinfo scaling from spm_vol struct and fname
for ivol = 26:32
    tempcopy.fname = sprintf('/Users/wbr/walter/fmri/sms_scan_analyses/rsa_singletrial/fir_sure/glasser_rsa_3network_movie_spm_%02d.nii',ivol);
    spm_write_vol(tempcopy,brain_4d(:,:,:,ivol))
end 
% hurray = make_nii(brain_4d);
% save_nii(hurray,'glasser_rsa_3network_movie.nii')

clear brain temp


fclose(fid_roi);