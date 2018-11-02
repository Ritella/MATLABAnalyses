% Performs inter-subject correlations (ISC) on fMRI data.
% See https://www.ncbi.nlm.nih.gov/pubmed/15016991 (Hasson, 2004)
% For each participant, every vertex timecourse is correlated to
% the timecourse of that vertex in the leave-one-out group.
% Pearson R values are then averaged across participants' for each vertex.
% This code can either perform ISC analysis on empirical fMRI values
% or it performs ISC analysis on empirical fMRI values that are  
% first phase shuffled to create null data.
% The latter is used to form a null distribution of R values, against which
% the statistical significance of empirical R is assessed 
% (i.e. a non-parametric, permutation analysis).

% Author: Rita Loiotile
% August 2018

function [] = performISC(server_name, batch_num) %#ok<FNDEF>

% INPUT
% server_name =  'marcc' (compute cluster) or 'arwen'; determines absolute paths set below                           
% batch_num =   only necesseary for permutation analysis, leave empty if
%               performing ISC on empirical data

% OUTPUT
% if empirical data, writes nifti files representing average Fisher
% transformed R for each vertex.
% if permutation analysis, writes csv file of same info on phase-shuffled
% data. also writes max R value for each hemi.

if nargin == 1 % ~exist(batch_num, 'var'); not sure why this doesnt work    
    clearvars -except server_name
    is_perm = 0;
else
    clearvars -except server_name batch_num
    is_perm = 1;    
end

if strcmp(server_name, 'arwen') 
    addpath('/media/BednyDrobo/Tools/matlab/spm12');
    addpath('/opt/fsl/5.0.9/etc/matlab/');
    path_pre = '/media/BednyDrobo/Projects/GNGC';
elseif strcmp(server_name, 'marcc') 
    addpath('/software/apps/fsl/5.0.11/etc/matlab/');
    path_pre = '/home-4/rloioti1@jhu.edu/work/GNGC';
else
    error('Need to set up path structure for that server. Current avalailable options are marcc and arwen.');
end

%% Declare Constants
FMRI_FILENAME_TEMPLATE = [path_pre '/subReplace/crosscorr/preproc-12fwhm/condReplace.feat/NIFTIs/hemiReplace.32k_fs_LR.surfed_data.func.nii.gz'];
GROUPS = {'S', 'CB'};
CONDS = {'rest', 'backward', 'scramble', 'pieman', 'phonecallhome', 'hauntedhouse', 'undercoverwire'};
HEMIS = {'lh', 'rh'};
SUBS = load([path_pre '/CrossCorr/subject_list.mat'], 'subs'); % list of participants to include for each condition
SUBS = SUBS.subs;
switch is_perm
    case 0
        OUT_DIR = [path_pre '/CrossCorr/isc/isc-12fwhm/condReplace.hemiReplace'];
    case 1
        OUT_DIR = [path_pre '/CrossCorr/isc/isc-12fwhm-perm/condReplace.hemiReplace'];
        rng(batch_num); % seed the random number generator so that not all batches are the same
end

tic; % start timer

%% ISC Analysis
%% Get Mean Group Timecourses and LOO Mean Group Timecourses
for g = 1:length(GROUPS)
    group = GROUPS{g};
    for c = 1:length(CONDS)
        cond = CONDS{c};
        for h = 1:length(HEMIS)        
            hemi = HEMIS{h};
            this_out_dir = regexprep(OUT_DIR, {'condReplace', 'hemiReplace'}, {cond, hemi});
            group_subs_by_cond = SUBS.(cond).(group);
            for i = 1:length(group_subs_by_cond)
                sub = group_subs_by_cond{i};
                [sub_tc, sub_tc_mask] = get_cleaned_fmri(sub, cond, hemi, FMRI_FILENAME_TEMPLATE); % get fmri timecourse (tc) for this scan
                % if permutation, shuffle each vertex timecourse
                % only for non-zero vertex or you'll get complex numbers
                if is_perm == 1
                    sub_tc = get_permuted_fmri(sub_tc, sub_tc_mask);
                end
                % add subject data to aggregate struct
                if i == 1  
                    tc.(group).(cond).(hemi).sum = sub_tc;
                    tc.(group).(cond).(hemi).count = sub_tc_mask;            
                else
                    tc.(group).(cond).(hemi).sum = tc.(group).(cond).(hemi).sum + sub_tc;
                    tc.(group).(cond).(hemi).count = tc.(group).(cond).(hemi).count + sub_tc_mask; 
                end
                % store subject data separately
                sub_data.(group).(cond).(hemi).(sub).tc = sub_tc;
                sub_data.(group).(cond).(hemi).(sub).mask = sub_tc_mask;                
                clear sub sub_tc sub_tc_mask
            end
            clear i
            % make mean timecourse file of all subs in a group.
            % no need to mask by vertices that contain all subs
            tc.(group).(cond).(hemi).mean = tc.(group).(cond).(hemi).sum ...
                                                ./ tc.(group).(cond).(hemi).count;
            tc.(group).(cond).(hemi).mask = (tc.(group).(cond).(hemi).count > 0);                               
            % save mean TC if empirical data
            if is_perm == 0
                write_image([this_out_dir '/NIFTIs/' group '_Mean_TC.nii.gz'], ...
                            tc.(group).(cond).(hemi).mean); % mean timecourse
                write_image([this_out_dir '/NIFTIs/' group '_Sub_Count.nii.gz'], ...
                            tc.(group).(cond).(hemi).count); % num of subs per vertex
            end  
            % Now make the leave-one-out group averages 
            % For each subject, subtract subject data from group             
            % Need to divide by number of subjects that contributed to the
            % data. 
            % If you calculate this number by subtracting 1 from the
            % total subject map count, it would be wrong since it assumes
            % that the current subject has contributed a value at every
            % vertex. Instead, you need to subtract the subject specific
            % contribution (either 0 or 1).
            for k = 1:length(group_subs_by_cond)
                sub = group_subs_by_cond{k};
                tc.(group).(cond).(hemi).for_sub.(sub).count = tc.(group).(cond).(hemi).count ...
                                                                - sub_data.(group).(cond).(hemi).(sub).mask;
                tc.(group).(cond).(hemi).for_sub.(sub).sum = tc.(group).(cond).(hemi).sum ...
                                                                - sub_data.(group).(cond).(hemi).(sub).tc; 
                tc.(group).(cond).(hemi).for_sub.(sub).mean = tc.(group).(cond).(hemi).for_sub.(sub).sum ...
                                                                ./ tc.(group).(cond).(hemi).for_sub.(sub).count;
                tc.(group).(cond).(hemi).for_sub.(sub).mask = (tc.(group).(cond).(hemi).for_sub.(sub).count > 0);
                if is_perm == 0
                    write_image([this_out_dir '/' group '_Sub_vs_' group '_Mean/NIFTIs/' group '_Mean_TC_for_' sub '.nii.gz'], ...
                                 tc.(group).(cond).(hemi).for_sub.(sub).mean); % mean LOO timecourse for sub
                    write_image([this_out_dir '/' group '_Sub_vs_' group '_Mean/NIFTIs/' group '_Sub_Count_for_' sub '.nii.gz'], ...
                                 tc.(group).(cond).(hemi).for_sub.(sub).count); % mean LOO timecourse for sub
                end
                clear sub
            end
            clear k hemi this_out_dir group_subs_by_cond
        end
        clear h cond
    end
    clear c group
end
clear g

%% Correlate Subject Timecourses to Group Timecourses
for g = 1:length(GROUPS)
    group1 = GROUPS{g}; 
    for g2 = 1:length(GROUPS) % 2 groups because we correlate both within and across groups
        group2 = GROUPS{g2};
        for c = 1:length(CONDS)
            cond = CONDS{c};
            group_subs_by_cond = SUBS.(cond).(group1); % group1 defines individual subjects
            for h = 1:length(HEMIS)        
                hemi = HEMIS{h};               
                this_out_dir = regexprep(OUT_DIR,{'condReplace','hemiReplace'},{cond,hemi});           
                for k = 1:length(group_subs_by_cond)
                    sub = group_subs_by_cond{k};
                    sub_tc = sub_data.(group1).(cond).(hemi).(sub).tc;
                    sub_mask = sub_data.(group1).(cond).(hemi).(sub).mask;
                    % Correlate each subject with each group's mean
                    % Importantly, use the LOO group if subject is
                    % part of their own group
                    if group1 == group2
                        to_correlate = tc.(group1).(cond).(hemi).for_sub.(sub).mean;
                        to_correlate_mask = tc.(group1).(cond).(hemi).for_sub.(sub).mask;
                    else
                        to_correlate = tc.(group2).(cond).(hemi).mean;
                        to_correlate_mask = tc.(group).(cond).(hemi).mask;
                    end
                    % To avoid this loop, you could use corcoeff and take the diagonal
                    % But that takes longer for this data size
                    sub_corr_to_group.val = zeros(size(sub_tc,1));
                    for a = 1:size(sub_tc,1) % for all vertices
                        sub_corr_to_group.val(a,1) = corr2(sub_tc(a,:), to_correlate(a,:));
                    end
                    clear a
                    sub_corr_to_group.val(isnan(sub_corr_to_group.val)) = 0; % get rid of NaNs now or will zero out for any missing sub vertex
                    sub_corr_to_group.mask = sub_mask(:,1) .* to_correlate_mask(:,1); % don't look for 0 corr since that can happen naturally
                    if k == 1
                        corr.(group1).(group2).(cond).(hemi).sum = sub_corr_to_group.val;
                        corr.(group1).(group2).(cond).(hemi).count = sub_corr_to_group.mask;
                    else
                        corr.(group1).(group2).(cond).(hemi).sum = corr.(group1).(group2).(cond).(hemi).sum ...
                                                                    + sub_corr_to_group.val;
                        corr.(group1).(group2).(cond).(hemi).count = corr.(group1).(group2).(cond).(hemi).count ...
                                                                    + sub_corr_to_group.mask;
                    end
                    if is_perm == 0
                        write_image([this_out_dir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/' sub '_corr.nii.gz'], ...
                                    sub_corr_to_group.val); % corr for sub
                    end
                    clear sub sub_tc sub_mask to_correlate to_correlate_mask sub_corr_to_group
                end
                clear k
                % get the mean corr for all subjects' corrs 
                % in this group and condition
                corr.(group1).(group2).(cond).(hemi).mean = corr.(group1).(group2).(cond).(hemi).sum ...
                                                            ./ corr.(group1).(group2).(cond).(hemi).count;
                corr.(group1).(group2).(cond).(hemi).mask = (corr.(group1).(group2).(cond).(hemi).count > 0);
                corr.(group1).(group2).(cond).(hemi).z_mean = r_to_z(corr.(group1).(group2).(cond).(hemi).mean);
                % REL TO DO: possibly convert R to Z at the individual
                % subject level, before meaning
                if perm == 0
                    write_image([this_out_dir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/Mean_Corr_Z.nii.gz'], corr.(group1).(group2).(cond).(hemi).z_mean);
                else
                    write_perm_info(batch_num, [this_out_dir '/' group1 '_Sub_vs_' group2 '_Mean'], ...
                                    corr.(group1).(group2).(cond).(hemi).z_mean, corr.(group1).(group2).(cond).(hemi).mask);
                end
                clear hemi
            end
            clear h cond this_out_dir group_subs_by_cond
        end
        clear c group2
    end
    clear g2 group1    
end
clear g

toc;
display(toc/60/60, 'hours');

end


function [z_value] = r_to_z(r_value)
    z_value = atanh(r_value);
end


function [r_value] = z_to_r(z_value)
    r_value = tanh(z_value);
end


%% Sub-Functions

function [sub_tc, sub_mask] = get_cleaned_fmri(sub, cond, hemi, fn_temp)  
    % gets EPI values and its associated mask of non-empty vertices
    sub_fn = regexprep(fn_temp, {'subReplace','condReplace','hemiReplace'}, {sub,cond,hemi}); 
    sub_tc = read_avw(sub_fn);
    sub_tc = sub_tc(:,:); % reduce from 4D to 2D (vertex, time)
    % find 0s and nans and mask out all 
    % associated timepoint for that vertex
    has_value = ~isnan(sub_tc(:,:)) .* ~(sub_tc(:,:)==0);
    sub_mask = (sum(has_value, 2) == size(has_value, 2)); % only use vertices that have data for all timepoints
    sub_mask = repmat(sub_mask, [1 size(has_value, 2)]);
    sub_tc = sub_tc .* sub_mask; % mask EPI
end


function permute_tc = get_permuted_fmri(empirical_tc, empirical_tc_mask)
    % get 2 permutation vectors for this subject's fMRI that will be used 
    % for all vertices
    [~, rp1] = sort(rand(round(size(empirical_tc,2)/2),1)); 
    [~, rp2] = sort(rand(round(size(empirical_tc,2)/2),1));
    permute_tc = zeros(size(empirical_tc));
    for a = 1:size(empirical_tc,1)
        if empirical_tc_mask(a) ~= 0 % don't permute null vertices
          permute_tc(a,:) = phase_rand(empirical_tc(a,:), rp1, rp2, true);
        end
    end
end


function [null_x] = phase_rand(x, rp1, rp2, permutation)
    % Permutes a timecourse (x) by shuffling phases of Fourier Transform.
    % This preserves dependencies amongst adjacent timepoints and mimics
    % empirically observed timecourses.
    % rp1 and rp2 are permutation orders for phases and phase signs,
    % respectively.
    % rp1 and rp2 should be round(length of x/2)
    % Author: Chris Honey (edited)
    % See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3787506/

    % param check
    if length(rp1) ~= round(length(x)/2) || length(rp2) ~= round(length(x)/2)
        error('permutation vector length must be half of timecourse length (or ceiling)');
    end

    x = x';

    [Nsamp, ~] = size(x);  %extract number of samples and number of signals

    x_mean = repmat(mean(x,1), Nsamp,1);
    x_std = repmat(std(x,1), Nsamp, 1); % original line assumes x has mean of 0; would have been fine to use if you subtracted x_mean before calculating; original: sqrt(repmat(dot(x,x), Nsamp, 1)/(Nsamp-1));
    x = x - x_mean;  % remove the mean of each column of X
    x = x./x_std; % divide by the standard deviation of each column

    %transform the vectors-to-be-scrambled to the frequency domain
    Fx = fft(x); 

    % identify indices of positive and negative frequency components of the fft
    % we need to know these so that we can symmetrize phase of neg and pos freq
    if mod(Nsamp,2) == 0
        posfreqs = 2:(Nsamp/2);
        negfreqs = Nsamp : -1 : (Nsamp/2)+2;
    else
        posfreqs = 2:(Nsamp+1)/2;
        negfreqs = Nsamp : -1 : (Nsamp+1)/2 + 1;
    end

    x_amp = abs(Fx);  %get the amplitude of the Fourier components

    if permutation
        x_phase = atan2(imag(Fx), real(Fx)); %get the phases of the Fourier components [NB: must use 'atan2', not 'atan' to get the sign of the angle right]
    end

    J = sqrt(-1);  %define the vertical vector in the complex plane

    if permutation  
         x_phase = x_phase(rp1);  % permute phases 
         sign = x_phase./abs(x_phase);
         perm_sign = sign(rp2); % permute the signs of phases
         x_phase = abs(x_phase) .* perm_sign; % change sign so added variability as if permuting second half in
         sym_phase(1) = x_phase(1,:);
         sym_phase(posfreqs,:) = x_phase(2:length(posfreqs)+1,:);
         sym_phase(round(Nsamp/2)+1) = x_phase(1,:); % gets reset if odd length
         sym_phase(negfreqs,:) = -x_phase(2:length(posfreqs)+1,:);
         sym_phase(isnan(sym_phase))=0;
    else
        new_phase = 2*pi*rand(length(posfreqs),K);
        sym_phase(posfreqs,:) = new_phase;
        sym_phase(negfreqs,:) = -new_phase;
    end

    z = x_amp.*exp(J.*sym_phase); % generate (symmetric)-phase-scrambled Fourier components
    null_x = real(ifft(z)); % invert the fft to generate a phase-scrambled version of x

    % Make null have same stdev and mean as original
    null_x_mean = repmat(mean(null_x,1), Nsamp,1);
    null_x_std = repmat(std(null_x,1), Nsamp, 1);
    null_x = null_x - null_x_mean;  %remove the mean of each column of null x
    null_x = null_x./null_x_std; %divide by the standard deviation of each column of null x
    null_x = null_x .* x_std; % now tranform back to the mean and stdev of x
    null_x = null_x + x_mean;
    null_x = null_x';
end


function write_image(fn, vals)
    [fp, ~, ~] = fileparts(fn);
    if ~exist(fp, 'dir')
        mkdir(fp);
    end 
    vals(isnan(vals))=0; % remove NaNs
    nifti(:,1,1,:)=vals(:,:); % map back to 4D structure
    save_avw(nifti, fn,'d', [1 1 1 2]); % 2 sets the TR
end


function [] = write_perm_info(batch_num, fn_root, vals, mask)
    % Because there are 1K+ permutations for each cond + group
    % This saves space and time by writing correlation maps as a csv 
    % file rather than a nifti
    % Also takes the max value across all vertices (in the hemi)
    % for use in multiple correction for vertex-wise stats
    
    % write out all vals 
    fn1 = [fn_root '/CSVs/Mean_corr_z_batch_' mat2str(batch_num) '.csv'];
    fid=fopen(fn1,'w');
    fprintf(fid, '%f \n', vals');
    fclose(fid);
    
    % descriptive stats
    fn2 = [fn_root '/Mean_corr_z_perm_stats.csv'];
    max_corr = max(vals);
    min_corr = min(vals);
    mean_corr = mean(vals(mask==1)); % mean and stdev only for non-zero vertices
    stdev_corr = std(vals(mask==1));
    
    % write stats to file
    fid=fopen(fn2,'a');
    fprintf(fid, '%d %f %f %f %f\n', [batch_num max_corr min_corr mean_corr stdev_corr]');
    fclose(fid);

end
