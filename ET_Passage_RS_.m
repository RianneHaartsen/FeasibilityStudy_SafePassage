
%% To run this script, you need to install Task Engine 2 v X.X
    % Email/Message Teresa Del Bianco
%% Description

    % This script reads the log and eye-tracking files for 1 participant and
    % generates a table summarising the missing data per trial, the centroid of
    % the fixations per trial and the standard distance deviation. List of
    % variables:

    %id = participant id 
    %event = stimuli file name
    %trial = trial n 1-4
    %onset = timestamp of movie onset
    %movie_duration = time the stimulus is on the screen in seconds (max 1 min)
    %et_looking_prop = proportion of gaze on the screen (et_looking_time) on movie_duration
    %et_looking_time = time the gaze is on the screen in seconds
    %nanP = percentage of lost gaze data (eyes not detected)
    %x_m = centroid of gaze x position
    %y_m = centroid of gaze y position
    %sdd = standard distance deviation of gaze position (root mean squared error)

%% Preliminary ops

warning('off','all')
clear variables

%% Set up directory

BasePath = '/Users/teresa/Documents/MATLAB/data/safepassage_eeg';
DirList  = dir(BasePath);
DirList  = DirList([DirList.isdir]);  % Folders only
    
% Set up local paths to scripts
    %add LM code TaskEngine2
    addpath(genpath('/Users/teresa/Documents/MATLAB/SP_Scripts/TaskEngine2'))
    addpath(genpath('/Users/teresa/Documents/MATLAB/SP_Scripts/lm_tools'));
    addpath(genpath('/Users/teresa/Documents/MATLAB/SP_Scripts/tasks'));  
    
%-------------------------
% participant folder level
%-------------------------
for iDir = 3:numel(DirList)
    
  %% Find participant directories
  aDir = fullfile(BasePath, DirList(iDir).name);
  SDirList = dir(aDir);
  SDirList(strncmp({SDirList.name}, '.', 1)) = []; % remove folders starting with '.'
  SDirList  = SDirList([SDirList.isdir]);
    %-------------------------
    % session folder level
    %-------------------------
  for jDir = 1:numel(SDirList) 
  %% Find session directories
      if numel(SDirList) == 0
        bDir = fullfile(BasePath, DirList(iDir).name); 
        fprintf('Empty folder: %s\n', bDir); % skip if participant folder empty
      elseif strncmp({SDirList.name}, 'SP', 1) 
        bDir = fullfile(BasePath, DirList(iDir).name); 
        fprintf('Missing session folder: %s\n', bDir); % skip if session folder not found
      else
      bDir = fullfile(BasePath, DirList(iDir).name, SDirList(jDir).name);
      fprintf('Processing: %s\n', bDir); 
      %% Read in data
        if isfile(sprintf('%s/tracker_passage_eeg_%s.mat', bDir, DirList(iDir).name)) % skip if et folder not found or object has no tasks property
            
          data = teData(bDir);
          %% Filter Log data
          tab = teLogFilter(data.Log, 'task', 'restingvideos2', 'topic', 'trial_log_data');
          t_n = height(tab); %n of trials
          %% Create empty table
          eyeT_info = table;
          
          %% Loop Through Trials if et folder exist
          if isfolder(sprintf('%s/eyetracking', bDir)) 
              % --------------------------
              % eyetracking folder level
              % --------------------------
              for k = 1:t_n
                  %% Select Onset and Offset of trial *1 & 2*
                  onset = tab.movie_onset{k}; %{trial n}
                  offset = tab.movie_offset{k};
                  s1 = find(data.ExternalData('eyetracking').Buffer(:, 1) >= onset, 1);
                  s2 = find(data.ExternalData('eyetracking').Buffer(:, 1) >= offset, 1);
                  %% Filter 1 trial
                  movieonscreen = data.ExternalData('eyetracking').Buffer(s1:s2, :);
                  %% Calculate mid point of eyes in X and Y
                  x_gaze = (movieonscreen(:, 2)+movieonscreen(:, 17))/2;
                  y_gaze = (movieonscreen(:, 3)+movieonscreen(:, 18))/2;
                  %% Calculate Standard Distance Deviation (SDD - root mean squared distance from the centroid)
                  x_m = nanmean(x_gaze); %centroid x
                  y_m = nanmean(y_gaze); %centroid y
                  x_sumSQ = nansum((x_gaze-x_m).^2); %sum of squares
                  y_sumSQ = nansum((y_gaze-y_m).^2);
                  x_count = sum(not(isnan(x_gaze))); %sum of cases (non NANs)
                  y_count = sum(not(isnan(y_gaze)));
                  x_d = x_sumSQ/x_count; 
                  y_d = y_sumSQ/y_count;
                  sdd = sqrt(x_d+y_d);
                  %% Create 1-row table
                  %t=1:height(tab);
                  id=cellstr(tab.id{k});
                  session=cellstr(SDirList(jDir).name);
                  event=cellstr(tab.key{k});
                  et_looking_prop=tab.et_looking_prop{k};
                  et_looking_time=tab.et_looking_time{k};
                  movie_duration=tab.movie_duration{k};
                  trial=k;
                  nanP=round((sum(movieonscreen(:,4)==0)/length(movieonscreen))*100, 2);
                  eyeT_info_temp = table(id,session,event,trial,onset,movie_duration,et_looking_prop,et_looking_time,nanP,x_m,y_m,sdd);
                  %% Populate table
                  eyeT_info = [eyeT_info;eyeT_info_temp]; 
                  %% Save Table to File
                  save(sprintf('ET_export/safe_passage/RS/matFiles/eyeT_info_%s_%s.mat',tab.id{k}, SDirList(jDir).name), 'eyeT_info')
                  writetable(eyeT_info,sprintf('ET_export/safe_passage/RS/csvFiles/eyeT_info_%s_%s.csv',tab.id{k}, SDirList(jDir).name),'Delimiter',',')
              end %end for k = 1:t_n
          else
              fprintf('Missing/truncated eyetracking folder: %s\n', bDir);
          end %end if isfolder(sprintf('%s/eyetracking', bDir, aDir(end-6:end)))
        else
            fprintf('Missing tracker file: %s\n', bDir);
        end %end if isfile(sprintf('%s/tracker_passage_eeg_%s.mat', bDir, aDir(end-6:end)))
      end %end if numel(SDirList) == 0
  end %end for jDir = numel(SDirList)
end %end for iDir = 4:numel(DirList)

warning('on','all')