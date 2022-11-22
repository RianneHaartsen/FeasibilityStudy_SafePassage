%% Safe Passage: 2) preparing EEG data for FOOOF analyses in python

% This script extracts the frequency spectra for all participants for each
% condition, from each region of interest for 1-30Hz. Note, this script
% includes the code used for the analyses in the publication.
% Input is the data from the script SafePassage_01_Preprocessing.m. 

% Created by Rianne Haartsen, February 2022

%% Prep re-ref to avg data for FOOOF in python

cd XXXX/SafePassage
load XXX/SP_preprocReREFavg_Power.mat
Path_1_30Hz = 'XXX/SP_1_30Hz';

for ss = 1:height(SP_preprocReREFavg_Power)
    if ~isempty(SP_preprocReREFavg_Power.PreprocEEG{ss})
        load(SP_preprocReREFavg_Power.PreprocEEG{ss},'Powerb_soctrls','Powerb_nsoctrls')
        Subj = SP_preprocReREFavg_Power.ID{ss};
        
        % for 1-30 Hz
        cd(Path_1_30Hz)
        % freq indices for 1-30Hz
        Frind_1_30Hz = [find(Powerb_soctrls.Pow_ft_struc.freq == 1):find(Powerb_soctrls.Pow_ft_struc.freq == 30)];
        % get abs power for FoI and save
        % social trials
        Powspec.Abspow = Powerb_soctrls.Abspowspec.Front(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_S_Fr.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_soctrls.Abspowspec.Cent(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_S_Ce.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_soctrls.Abspowspec.Pari(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_S_Pa.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_soctrls.Abspowspec.Occip(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_S_Oc.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_soctrls.Abspowspec.WholeScalp'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab_S_WhSc.csv'))
        clear Tab Powspec
        % non-social trials
        Powspec.Abspow = Powerb_nsoctrls.Abspowspec.Front(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_NS_Fr.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_nsoctrls.Abspowspec.Cent(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_NS_Ce.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_nsoctrls.Abspowspec.Pari(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_NS_Pa.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_nsoctrls.Abspowspec.Occip(1,Frind_1_30Hz)'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab1_30_NS_Oc.csv'))
        clear Abspow_soc_front Tab Powspec
        Powspec.Abspow = Powerb_nsoctrls.Abspowspec.WholeScalp'; 
        Tab = struct2table(Powspec);
        writetable(Tab, strcat(Subj,'_Ab_NS_WhSc.csv'))
        clear Tab Powspec

        clear data and subj
        clear Powerb_soctrls Powerb_nsoctrls Subj 
    end
end


%% Safe Passage: randomly select power spectra for FOOOF
% 2 conditions, 4 regions, 76 subjects: 2*4*76 = 608 spectra
% 608*.2 = 122 power spectra

% find names of power spectra
Path_1_30Hz = 'XXX/SP_1_30Hz';
cd(Path_1_30Hz)
List_allspectra = dir('*.csv');
Ind_spectraoi = zeros(size(List_allspectra,1),1);
for ii = 1:size(List_allspectra,1)
    oi = contains(List_allspectra(ii).name,'WhSc');
    if isequal(oi,0)
        Ind_spectraoi(ii) = 1;
    end
end
Inds_oi = find(Ind_spectraoi == 1);
Random_inds = randperm(size(Inds_oi,1),122);
Random_spectra_inds = Inds_oi(Random_inds);

List_randomspectra = List_allspectra(Random_spectra_inds);


% check for about even number of condition * region combinations
count_S_Fr = 0;
count_NS_Fr = 0;
count_S_Ce = 0;
count_NS_Ce = 0;
count_S_Pa = 0;
count_NS_Pa = 0;
count_S_Oc = 0;
count_NS_Oc = 0;

for ii = 1:size(List_randomspectra,1)
    if contains(List_randomspectra(ii).name,'_S_')
        if contains(List_randomspectra(ii).name,'_Fr.csv')
            count_S_Fr = count_S_Fr+1;
        elseif contains(List_randomspectra(ii).name,'_Ce.csv')
            count_S_Ce = count_S_Ce+1;
        elseif contains(List_randomspectra(ii).name,'_Pa.csv')
            count_S_Pa = count_S_Pa+1;
        elseif contains(List_randomspectra(ii).name,'_Oc.csv')
            count_S_Oc = count_S_Oc+1;
        end
        
    elseif contains(List_randomspectra(ii).name,'_NS_')
        if contains(List_randomspectra(ii).name,'_Fr.csv')
            count_NS_Fr = count_NS_Fr+1;
        elseif contains(List_randomspectra(ii).name,'_Ce.csv')
            count_NS_Ce = count_NS_Ce+1;
        elseif contains(List_randomspectra(ii).name,'_Pa.csv')
            count_NS_Pa = count_NS_Pa+1;
        elseif contains(List_randomspectra(ii).name,'_Oc.csv')
            count_NS_Oc = count_NS_Oc+1;
        end
    else
        disp('unidentified condition')
    end
end

% report numbers
fprintf('Region - Social - NonSocial\n')
fprintf('Frontal - %d - %d\n',count_S_Fr, count_NS_Fr)
fprintf('Central - %d - %d\n',count_S_Ce, count_NS_Ce)
fprintf('Parietal - %d - %d\n',count_S_Pa, count_NS_Pa)
fprintf('Occipital - %d - %d\n',count_S_Oc, count_NS_Oc)

%% Copy random spectra to training folder

for ii = 1:size(List_randomspectra)
    curfile = List_randomspectra(ii).name;
    copyfile(curfile, 'XXX/SP_PowSpectra_trainingdata')
    clear curfile
end

Counts.S_Fr = count_S_Fr;
Counts.NS_Fr = count_NS_Fr;
Counts.S_Ce = count_S_Ce;
Counts.NS_Ce = count_NS_Ce;
Counts.S_Pa = count_S_Pa;
Counts.NS_Pa = count_NS_Pa;
Counts.S_Oc = count_S_Oc;
Counts.NS_Oc = count_NS_Oc;

cd XXX/SP_PowSpectra_trainingdata
save('List_randomspectra.mat','List_randomspectra','Counts');
