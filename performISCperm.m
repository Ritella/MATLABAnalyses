% Correlates subjects to each other, according to inter-subject correlation (ISC).

function []=performISCperm(batch_num)
clearvars -except batch_num
addpath('/media/BednyDrobo/Tools/matlab/spm12');
addpath('/opt/fsl/5.0.9/etc/matlab/');

%% DECLARE CONSTANTS
movies={'phonecallhome'};
groups={'S'};
hemis={'lh'};
mainDir='/media/BednyDrobo/Projects/GNGC/subReplace/crosscorr/preproc-12fwhm/movieReplace.feat/NIFTIs/hemiReplace.32k_fs_LR.surfed_data.func.nii.gz';
groupDir='/media/BednyDrobo/Projects/GNGC/CrossCorr/isc/isc-12fwhm-perm-test/movieReplace.hemiReplace';
%movies={'rest', 'backward', 'scramble', 'pieman', 'phonecallhome', 'hauntedhouse', 'undercoverwire'};
%groups={'S', 'CB'};
%hemis={'lh', 'rh'};
%mainDir='/media/BednyDrobo/Projects/GNGC/subReplace/crosscorr/preproc-12fwhm/movieReplace.feat/NIFTIs/hemiReplace.32k_fs_LR.surfed_data.func.nii.gz';
%groupDir='/media/BednyDrobo/Projects/GNGC/CrossCorr/isc/isc-12fwhm-perm/movieReplace.hemiReplace';
subs = load('/media/BednyDrobo/Projects/GNGC/CrossCorr/subject_list.mat', 'subs');
subs = subs.subs;

tic;
rng(batch_num+now); % seed the random number generator so that not all batches are the same
%% START ANALYSIS
for r=1:length(groups)
    group=groups{r};
    for i = 1:length(movies)
        movie=movies{i};
        for j=1:length(hemis)        
            hemi=hemis{j};
            thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{movie,hemi});
            meanSubTC.(group).(movie).(hemi).count=[];
            meanSubTC.(group).(movie).(hemi).value=[];
            meanSubTC.(group).(movie).(hemi).sum=[];
            theseSubs=subs.(movie).(group);
            for k = 1:length(theseSubs)
                sub=theseSubs{k};
                subFile=mainDir;
                subFile=regexprep(subFile,{'subReplace','movieReplace','hemiReplace'},{sub,movie,hemi}); 
                %display(subFile);
                subValue=read_avw(subFile);
                subValue=subValue(:,:);
                % find 0s and nans and mask out all associated timepoint
                % for that vertex
                nonNullCount=~isnan(subValue(:,:)) .* ~(subValue(:,:)==0);
                subMask = (sum(nonNullCount,2) == size(nonNullCount, 2));
                subMask = repmat(subMask, [1 size(nonNullCount,2)]);
                subValue = subValue .* subMask;
                subCount= subMask;
                % permute subject timecourse by vertex
                % only for non-zero vertex or you'll get complex numbers
                for aa = 1:size(subValue,1)
                    if subMask(aa)~=0
                      subValue(aa,:)=phase_rand(subValue(aa,:));
                    end
                end
                if k==1  % change this so that it checks whether it exists
                    meanSubTC.(group).(movie).(hemi).sum=subValue;
                    meanSubTC.(group).(movie).(hemi).count=subCount;            
                else
                    meanSubTC.(group).(movie).(hemi).sum=meanSubTC.(group).(movie).(hemi).sum+subValue;
                    meanSubTC.(group).(movie).(hemi).count=meanSubTC.(group).(movie).(hemi).count + subCount; 
                end
                % Prep the Leave-One-Out by Subtracting out the Subject
                meanSubTC.(group).(movie).(hemi).sub.(sub)= subValue;
                meanSubTC.(group).(movie).(hemi).subCount.(sub)= subCount;                
                rawSubValue.(group).(movie).(hemi).sub.(sub) =  subValue;
                rawSubValue.(group).(movie).(hemi).mask.(sub) =  subMask;
                clear sub subFile subValue subCount subMask nonNullCount
            end
            % Make Mean Timecourse file of all Subs in a Group. No need to
            % mask by vertices that contain all subs
            meanSubTC.(group).(movie).(hemi).value = meanSubTC.(group).(movie).(hemi).sum ./ meanSubTC.(group).(movie).(hemi).count;
            if ~exist([thisDir '/NIFTIs/'], 'dir')
                mkdir([thisDir '/NIFTIs/']);
            end 
            % save mean TC
            meanSubTC.(group).(movie).(hemi).value(isnan(meanSubTC.(group).(movie).(hemi).value))=0;
            %nifti(:,1,1,:)=meanSubTC.(group).(movie).(hemi).value(:,:);
            %save_avw(nifti,[thisDir '/NIFTIs/' group '_Mean_TC.nii.gz'],'d',[1 1 1 2]);
            %fprintf([thisDir '/' group '_Mean_TC.gii\n']);          
            %clear nifti;
            % save subject count mask
            %nifti(:,1,1,:)=meanSubTC.(group).(movie).(hemi).count(:,:);
            %save_avw(nifti,[thisDir '/NIFTIs/' group '_SubCount.nii.gz'],'d',[1 1 1 2]);          
            %clear nifti;
            % Now make the Subject Specific Group Averages for
            % Leave-One-Out
            % If you add the SubValue, it will subtract it from the total
            % sum
            % Need to divide by number of subjects that contributed to the
            % data.  If you calculate this number by subtracting 1 from the
            % total subject map count, it would be wrong.  That assumes
            % that the current subject has contributed a value at every
            % vertex.  Instead, you need to subtract the subject specific
            % contribution (either 0 or 1).
            for k = 1:length(theseSubs)
                sub=theseSubs{k};
                meanSubTC.(group).(movie).(hemi).subCount.(sub) = (meanSubTC.(group).(movie).(hemi).count - meanSubTC.(group).(movie).(hemi).subCount.(sub));
                meanSubTC.(group).(movie).(hemi).sub.(sub)= (meanSubTC.(group).(movie).(hemi).sum - meanSubTC.(group).(movie).(hemi).sub.(sub)) ./ meanSubTC.(group).(movie).(hemi).subCount.(sub);
                meanSubTC.(group).(movie).(hemi).sub.(sub)(isnan(meanSubTC.(group).(movie).(hemi).sub.(sub))) = 0;
                if ~exist([thisDir '/' group '_Sub_vs_' group '_Mean/NIFTIs/'], 'dir')
                    mkdir([thisDir '/' group '_Sub_vs_' group '_Mean/NIFTIs/']);
                end
                %nifti(:,1,1,:)=meanSubTC.(group).(movie).(hemi).sub.(sub)(:,:);
                %nifti(isnan(nifti))=0;
                %save_avw(nifti,[thisDir '/' group '_Sub_vs_' group '_Mean/NIFTIs/' group '_Mean_TC_for_' sub '.nii.gz'],'d',[1 1 1 2]);
                %fprintf([thisDir '/' group '_Sub_vs_' group '_Mean/' group '_Mean_TC_for_' sub '.gii\n']);
                %clear nifti;
                % save subject count mask
                %nifti(:,1,1,:)=meanSubTC.(group).(movie).(hemi).subCount.(sub)(:,:);
                %save_avw(nifti,[thisDir '/' group '_Sub_vs_' group '_Mean/NIFTIs/' group '_Sub_Count_for_' sub '.nii.gz'],'d',[1 1 1 2]);          
                %clear nifti;
            end
            clear k theseSubs thisDir
        end
        clear j hemi
    end
    clear i movie
end
clear r group
%save('ISC-workspace.mat');

for r=1:length(groups)
    group1=groups{r};
    for v=1:length(groups)
        group2=groups{v};
        for i = 1:length(movies)
            movie=movies{i};
            theseSubs=subs.(movie).(group1);
            for j=1:length(hemis)        
                hemi=hemis{j};               
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{movie,hemi});           
                for k = 1:length(theseSubs)
                    sub = theseSubs{k};
                    subValue = rawSubValue.(group1).(movie).(hemi).sub.(sub);
                    subValueMask = rawSubValue.(group1).(movie).(hemi).mask.(sub);
                    % Correlate each subject with each group's mean
                    % Importantly, the procedure is different if subject is
                    % part of their own group
                    if group1==group2
                        toCorrelate = meanSubTC.(group1).(movie).(hemi).sub.(sub);
                        toCorrelateMask=(meanSubTC.(group1).(movie).(hemi).subCount.(sub) > 0);
                    else
                        toCorrelate = meanSubTC.(group2).(movie).(hemi).value;
                        toCorrelateMask=(meanSubTC.(group2).(movie).(hemi).count > 0);
                    end
                    % To avoid this loop, you could calcute corr by hand,
                    % like with fslmaths           
                    for a=1:size(subValue,1)
                        subValueCorrGroup(a,1)=corr2(subValue(a,:),toCorrelate(a,:));
                    end
                    subValueCorrGroup(isnan(subValueCorrGroup))=0;
                    if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs'], 'dir')
                        mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs']);
                    end
                    %nifti(:,1,1,:)=subValueCorrGroup(:,:);           
                    %save_avw(nifti,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/' sub '_corr.nii.gz'],'d',[1 1 1 2]);
                    %fprintf('save nifti\n');
                    %clear nifti;
                    %nifti(:,1,1,:)=1-subValueCorrGroup(:,:).^2;
                    %nifti(isnan(nifti))=0;
                    %save_avw(nifti,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/' sub '_resid.nii.gz'],'d',[1 1 1 2]);
                    %fprintf('save nifti\n');
                    %clear nifti;
                    %nifti(:,1,1,:)=(size(subValue,2)-2) .* ones([size(subValue,1),1]); % This mask should be zeroed-out where there are no values
                    %nifti(isnan(nifti))=0;
                    %save_avw(nifti,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/' sub '_DOFs.nii.gz'],'d',[1 1 1 2]);
                    %fprintf('save nifti\n');
                    %clear nifti toCorrelate  
                    subValueCorrGroup = subValueCorrGroup .* (subValueMask(:,1) .* toCorrelateMask(:,1)); % to be safe. note that it's not good to look for corr values = 0 since those may occur naturally
                    if k==1
                        mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).sum = subValueCorrGroup;
                        mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).count = subValueMask(:,1) .* toCorrelateMask(:,1); % both vectors have to be non-zero
                    else
                        mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).sum = mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).sum + subValueCorrGroup;
                        mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).count = mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).count + (subValueMask(:,1) .* toCorrelateMask(:,1));
                    end
                    clear sub subFile subValue subValueCorrGroup subValueMask toCorrelateMask
                end
                clear k
            end
            clear j hemi thisDir
        end
        clear i movie theseSubs
    end
    clear v group2      
end
clear r group1

% Now aggregate (average) all single-subject to group correlations    
for r=1:length(groups)
    group1=groups{r};
    for v=1:length(groups)
        group2=groups{v};
        for i = 1:length(movies)
            movie=movies{i};
            for j=1:length(hemis)        
                hemi=hemis{j};
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{movie,hemi});  
                mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value = mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).sum ./ mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).count;
                mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value(isnan(mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value))=0;
                mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).mask = (mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).count > 0);
                mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value = mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value .* mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).mask; % to be safe
                nifti1(:,1,1,:)=mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).zvalue = rtoz(mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).value);
                nifti2(:,1,1,:)=mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).zvalue(:,:);             
                write_max(batch_num, nifti1, nifti2,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean'], mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).mask);
                clear nifti1;
                clear nifti2;
                %nifti(:,1,1,:)=mean_corr_sub_by_group.(movie).(hemi).(group1).(group2).count(:,:);
                %save_avw(nifti,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/Sub_Count_for_mean_corr.nii.gz'],'d',[1 1 1 2]);
                %clear nifti;
            end
            clear j hemi 
        end
        clear i movie thisDir
    end
    clear v group2      
end
clear r group1

% Now aggregate across and compare conditions
for r=1:length(groups)
    group1=groups{r};
    for v=1:length(groups)
        group2=groups{v};
            for j=1:length(hemis)        
                hemi=hemis{j};
                % Mean Movies
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).sum  = mean_corr_sub_by_group.('phonecallhome').(hemi).(group1).(group2).value + mean_corr_sub_by_group.('hauntedhouse').(hemi).(group1).(group2).value + mean_corr_sub_by_group.('undercoverwire').(hemi).(group1).(group2).value;
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).count = ...
                    (mean_corr_sub_by_group.('phonecallhome').(hemi).(group1).(group2).mask ...
                    + mean_corr_sub_by_group.('hauntedhouse').(hemi).(group1).(group2).mask ...
                    + mean_corr_sub_by_group.('undercoverwire').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask = (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).count > 0);
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).sum ./ mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).count;
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value =  mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value .*  mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask;
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value(isnan(mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value)) = 0;
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'meanmovies',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zsum  = mean_corr_sub_by_group.('phonecallhome').(hemi).(group1).(group2).zvalue + mean_corr_sub_by_group.('hauntedhouse').(hemi).(group1).(group2).zvalue + mean_corr_sub_by_group.('undercoverwire').(hemi).(group1).(group2).zvalue;
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zsum ./ mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).count;
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue =  mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue .*  mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask;
                mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue(isnan(mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue)) = 0;
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue(:,:);    
                write_max(batch_num, nifti1, nifti2, [thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
                % Mean Movies - Backward
                mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).value  = (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value - mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).value);
                mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).mask = ...
                    (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask .* ...
                     mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).value .* mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).mask;
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'meanmovies_backward',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).zvalue  = (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue - mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).zvalue);
                mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).zvalue .* mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).mask;
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).zvalue(:,:);   
                write_max(batch_num, nifti1, nifti2, [thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('meanmovies_backward').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
                % Pie-Man - Backward
                mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).value  = (mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).value - mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).value);
                mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).mask = ...
                    (mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).mask .* ... 
                     mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).value .* mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).mask;
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'pieman_backward',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).zvalue  = (mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).zvalue - mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).zvalue);
                mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).zvalue .* mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).mask;
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).zvalue(:,:);
                write_max(batch_num, nifti1, nifti2, [thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('pieman_backward').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
                % Scramble - Backward
                mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).value  = (mean_corr_sub_by_group.('scramble').(hemi).(group1).(group2).value - mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).value);
                mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).mask = ...
                    (mean_corr_sub_by_group.('scramble').(hemi).(group1).(group2).mask ...
                    .* mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).value .* mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).mask; 
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'scramble_backward',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).zvalue  = (mean_corr_sub_by_group.('scramble').(hemi).(group1).(group2).zvalue - mean_corr_sub_by_group.('backward').(hemi).(group1).(group2).zvalue);         
                mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).zvalue .* mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).mask; 
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).zvalue(:,:);
                write_max(batch_num, nifti1,nifti2,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('scramble_backward').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
                % Pie-Man - Scramble
                mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).value  = (mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).value - mean_corr_sub_by_group.('scramble').(hemi).(group1).(group2).value);
                mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).mask = ...
                    (mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).mask ...
                    .* mean_corr_sub_by_group.('scramble').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).value .* mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).mask;
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'pieman_scramble',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).zvalue  = (mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).zvalue - mean_corr_sub_by_group.('scramble').(hemi).(group1).(group2).zvalue);         
                mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).zvalue .* mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).mask;
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).zvalue(:,:);
                write_max(batch_num, nifti1, nifti2,[thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('pieman_scramble').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
                % Mean Movies - PieMan
                mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).value  = (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value - mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).value);
                mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).mask = ...
                    (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask ...
                    .* mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).value .* mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).mask;
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'meanmovies_pieman',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).zvalue  = (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue - mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).zvalue);
                mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).zvalue .* mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).mask;
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).zvalue(:,:);
                write_max(batch_num, nifti1, nifti2, [thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('meanmovies_pieman').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
                % PieMan - Mean Movies 
                mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).value  = (-mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).value + mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).value);
                mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).mask = ...
                    (mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).mask ...
                    .* mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).mask);
                mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).value = mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).value .* mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).mask;
                thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{'pieman_meanmovies',hemi});  
                if ~exist([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/'], 'dir')
                    mkdir([thisDir '/' group1 '_Sub_vs_' group2 '_Mean/NIFTIs/']);
                end 
                nifti1(:,1,1,:) = mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).value(:,:);
                mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).zvalue  = (-mean_corr_sub_by_group.('meanmovies').(hemi).(group1).(group2).zvalue + mean_corr_sub_by_group.('pieman').(hemi).(group1).(group2).zvalue);
                mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).zvalue = mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).zvalue .* mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).mask;
                nifti2(:,1,1,:) = mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).zvalue(:,:);
                write_max(batch_num, nifti1, nifti2, [thisDir '/' group1 '_Sub_vs_' group2 '_Mean'],mean_corr_sub_by_group.('pieman_meanmovies').(hemi).(group1).(group2).mask);
                clear nifti1 nifti2;
            end
            clear j hemi 
    end
    clear v group2      
end
clear r group1

% Now average CB-S and S-CB
conditions = {movies{:}, 'meanmovies', 'meanmovies_backward', 'pieman_backward', 'scramble_backward', 'pieman_scramble', 'meanmovies_pieman', 'pieman_meanmovies'};
for i=1:length(conditions)
    cond=conditions{i};
    for j=1:length(hemis)        
        hemi=hemis{j};
        % Mean CB-S and S-CB
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').value  = (mean_corr_sub_by_group.(cond).(hemi).('CB').('S').value + mean_corr_sub_by_group.(cond).(hemi).('S').('CB').value)/2;
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').mask = ...
                    (mean_corr_sub_by_group.(cond).(hemi).('CB').('S').mask .* ...
                    mean_corr_sub_by_group.(cond).(hemi).('S').('CB').mask);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').value = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').value .* mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').mask;
        thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{cond,hemi});  
        if ~exist([thisDir '/' 'CB_and_S' '_Mean/NIFTIs/'], 'dir')
            mkdir([thisDir '/' 'CB_and_S' '_Mean/NIFTIs/']);
        end 
        nifti1(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').value(:,:);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').zvalue  = (mean_corr_sub_by_group.(cond).(hemi).('CB').('S').zvalue + mean_corr_sub_by_group.(cond).(hemi).('S').('CB').zvalue)/2;
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').zvalue = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').zvalue .* mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').mask;        
        nifti2(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').zvalue(:,:);
        write_max(batch_num, nifti1,nifti2,[thisDir '/' 'CB_and_S' '_Mean'],mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').mask);
        clear nifti1 nifti2;
        % CB - S (includes interactions)
        mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').value  = (mean_corr_sub_by_group.(cond).(hemi).('CB').('CB').value - mean_corr_sub_by_group.(cond).(hemi).('S').('S').value);
        mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').mask = ...
                    (mean_corr_sub_by_group.(cond).(hemi).('CB').('CB').mask .* ...
                     mean_corr_sub_by_group.(cond).(hemi).('S').('S').mask);
        mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').value = mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').value .* mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').mask;      
        thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{cond,hemi});  
        if ~exist([thisDir '/' 'CB_minus_S' '_Mean/NIFTIs/'], 'dir')
            mkdir([thisDir '/' 'CB_minus_S' '_Mean/NIFTIs/']);
        end 
        nifti1(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').value(:,:);
        mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').zvalue  = (mean_corr_sub_by_group.(cond).(hemi).('CB').('CB').zvalue - mean_corr_sub_by_group.(cond).(hemi).('S').('S').zvalue);
        mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').zvalue = mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').zvalue .* mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').mask;      
        nifti2(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').zvalue(:,:);       
        write_max(batch_num, nifti1, nifti2, [thisDir '/' 'CB_minus_S' '_Mean'],mean_corr_sub_by_group.(cond).(hemi).('CB_S').('CB_S').mask);
        clear nifti1 nifti2;      
        % CB_and_S - S (includes interactions)
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').value  = (-mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').value + mean_corr_sub_by_group.(cond).(hemi).('S').('S').value);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').mask = ...
                    (mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').mask .* ...
                    mean_corr_sub_by_group.(cond).(hemi).('S').('S').mask);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').value = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').value .* mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').mask;        
        thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{cond,hemi});  
        if ~exist([thisDir '/' 'CB_and_S_minus_S' '_Mean/NIFTIs/'], 'dir')
            mkdir([thisDir '/' 'CB_and_S_minus_S' '_Mean/NIFTIs/']);
        end 
        nifti1(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').value(:,:);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').zvalue  = (-mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').zvalue + mean_corr_sub_by_group.(cond).(hemi).('S').('S').zvalue);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').zvalue = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').zvalue .* mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').mask;        
        nifti2(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').zvalue(:,:);
        write_max(batch_num, nifti1, nifti2, [thisDir '/' 'CB_and_S_minus_S' '_Mean'],mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_S').('CB_and_S_S').mask);
        clear nifti1 nifti2;  
        % CB_and_S - CB (includes interactions)
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').value  = (-mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').value + mean_corr_sub_by_group.(cond).(hemi).('CB').('CB').value);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').mask = ...
                   (mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').mask .* ...
                    mean_corr_sub_by_group.(cond).(hemi).('CB').('CB').mask);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').value = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').value .* mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').mask;
        thisDir=regexprep(groupDir,{'movieReplace','hemiReplace'},{cond,hemi});  
        if ~exist([thisDir '/' 'CB_and_S_minus_CB' '_Mean/NIFTIs/'], 'dir')
            mkdir([thisDir '/' 'CB_and_S_minus_CB' '_Mean/NIFTIs/']);
        end 
        nifti1(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').value(:,:);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').zvalue  = (-mean_corr_sub_by_group.(cond).(hemi).('CB_and_S').('CB_and_S').zvalue + mean_corr_sub_by_group.(cond).(hemi).('CB').('CB').zvalue);
        mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').zvalue = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').zvalue .* mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').mask;
        nifti2(:,1,1,:) = mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').zvalue(:,:);
        write_max(batch_num, nifti1, nifti2,[thisDir '/' 'CB_and_S_minus_CB' '_Mean'],mean_corr_sub_by_group.(cond).(hemi).('CB_and_S_CB').('CB_and_S_CB').mask);
        clear nifti1 nifti2;  
    end
    clear j hemi 
end
clear i cond      

toc;
display(toc/60/60, 'hours');

end


function [zvalue] = rtoz(rvalue)
    zvalue=atanh(rvalue);
end


function [rvalue] = ztor(zvalue)
    rvalue=tanh(zvalue);
end


function [null_x] = phase_rand(x)
x = x';
permutation=true;

[Nsamp K] = size(x);  %extract number of samples and number of signals

x_mean = repmat(mean(x,1), Nsamp,1);
x_std = repmat(std(x,1), Nsamp, 1); % original line assumes x has mean of 0; would have been fine to use if you subtracted x_mean before calculating; original: sqrt(repmat(dot(x,x), Nsamp, 1)/(Nsamp-1));
x = x - x_mean;  %remove the mean of each column of X
x = x./x_std; %divide by the standard deviation of each column

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
    [tmp,rp] = sort(rand(Nsamp,K));
    x_phase=x_phase(rp);
    sym_phase(posfreqs,:) = x_phase(1:length(posfreqs),:);
    sym_phase(negfreqs,:) = -x_phase(1:length(posfreqs),:);
else
    new_phase=2*pi*rand(length(posfreqs),K);
    sym_phase(posfreqs,:) = new_phase;
    sym_phase(negfreqs,:) = -new_phase;
end
    
z = x_amp.*exp(J.*sym_phase); %generate (symmetric)-phase-scrambled Fourier components
null_x = ifft(z); %invert the fft to generate a phase-scrambled version of x

% Make null have same stdev and mean as original
null_x_mean = repmat(mean(null_x,1), Nsamp,1);
null_x_std = repmat(std(null_x,1), Nsamp, 1);
null_x = null_x - null_x_mean;  %remove the mean of each column of Null X
null_x = null_x./null_x_std; %divide by the standard deviation of each column of Null X
null_x = null_x .* x_std; % now tranform back to the mean and stdev of X
null_x = null_x + x_mean;

%display(corr2(x,null_x));

null_x = null_x';
end

function [] = write_max(batch_num, nifti, znifti, file_dir, mask)
% write out max value and other stats
filename = [file_dir '/Mean_corr_perm_stats.csv'];
max_nifti_corr = max(nifti);
min_nifti_corr = min(nifti);
mean_nifti_corr = mean(nifti);
std_nifti_corr = std(nifti);
fid=fopen(filename,'a');
fprintf(fid, '%d %f %f %f %f\n', [batch_num max_nifti_corr min_nifti_corr mean_nifti_corr std_nifti_corr]');
fclose(fid);

% write out max value and other stats
filename = [file_dir '/Mean_corr_z_perm_stats.csv'];
max_znifti_corr = max(znifti);
min_znifti_corr = min(znifti);
mean_znifti_corr = mean(znifti);
std_znifti_corr = std(znifti);
fid=fopen(filename,'a');
fprintf(fid, '%d %f %f %f %f\n', [batch_num max_znifti_corr min_znifti_corr mean_znifti_corr std_znifti_corr]');
fclose(fid);

% write out all vals (masked so that you're not counting non-true zero
% correlations)
mask = ((mask>0) .* ~isnan(nifti));
nifti_masked = nifti(logical(mask));
filename = [file_dir '/NIFTIs/Mean_corr_batch_' mat2str(batch_num) '.csv'];
fid=fopen(filename,'w');
fprintf(fid, '%f \n', [nifti_masked]');
fclose(fid);

% write out all vals (masked so that you're not counting non-true zero
% correlations)
zmask = ((mask>0) .* ~isnan(znifti));
znifti_masked = znifti(logical(zmask));
zfilename = [file_dir '/NIFTIs/Mean_corr_z_batch_' mat2str(batch_num) '.csv'];
fid=fopen(zfilename,'w');
fprintf(fid, '%f \n', [znifti_masked]');
fclose(fid);
end

