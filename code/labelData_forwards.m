function Outputs = labelData_forwards(patientList,dataPath,...
    searchWindowTime,parameterName,parameterThresholdsToTest,alertIfParamGreaterThanThresh,...
    excludeInterventions,alternateInterventionDefinitionFlag,...
    keepDetailedTable)

%%example how to use:
% %---set inputs for labeling
% %which patients data to consider
% [masterPatientList] = helpers.GetListOfPatients();
% patientList = masterPatientList.Total_FT;
% patientGroupName = 'FT';
% 
% dataPath = '..\data'; % where to load the patient data files
% searchWindowTime = 10; %candidates are: 5,10,15.  Length in minutes of search window for forward analysis
% 
% %which parameter to evaluate
% parameterName = 'HPI'; %candidates are: HPI,MAP,CO,SV,PulsePressure,HR,SVV,ShockIndex,dynEa,ShockIndexMap,deltaMap65to75,deltaMap75to85,deltaMap85to95
% 
% parameterThresholdsToTest = [];
% alertIfParamGreaterThanThresh = 0;%1 means Parameter > thresh is an alert, 0 means Parameter < thresh is an alert
% 
% %interventions
% excludeInterventions = 1; % 1 to exclude interventions, 0 to keep
% alternateInterventionDefinitionFlag = 0; %0: intervention is when MAP increased by >5mmHg within 20 seconds or when MAP increased >8mmHg within 2 minutes; 1: uses 10mmHg and 10mmHg
% 
% keepDetailedTable = 0;  % 1 to save the idx and patientID,MAP,&HPI value of each TP,FP,etc.....makes it much slower
% OutputsLabels = labelData_forwards(patientList,dataPath,...
%     searchWindowTime,parameterName,parameterThresholdsToTest,alertIfParamGreaterThanThresh,...
%     excludeInterventions,alternateInterventionDefinitionFlag,...
%     keepDetailedTable);
%--------------------------------------------------------------------------
%--mandatory inputs
%patientList: nx1 cell array, filenames of patient data (without the .mat file extension)
%dataPath: char, file path containing the patient data files
%searchWindowTime: length in minutes of search window for forward analysis
%parameterName: name of the variable to predict hypotension
%
%--handle optional flags
%parameterThresholdsToTest: parameter compared to these thresholds to indicate an alert
if ~exist('parameterThresholdsToTest', 'var') || isempty(parameterThresholdsToTest)
    switch lower(parameterName)
        case 'hpi'%Hypotension Prediction Index
            parameterThresholdsToTest = 0:1:100;
        case 'map'%Mean Arterial Pressure
            parameterThresholdsToTest = [10:49 50:0.5:110 111:150]; %10:1:150;
        case 'co'%Cardiac Output
            parameterThresholdsToTest = 0.1:0.1:20;
        case 'sv'%Stroke Volume
            parameterThresholdsToTest = 10:1:200;
        case 'pulsepressure'%Pulse Pressure
            parameterThresholdsToTest = 10:1:200;
        case 'hr'%Heart Rate
            parameterThresholdsToTest = 10:1:200;
        case 'svv'%Stroke Volume Variation
            parameterThresholdsToTest = 0:1:100;
        case 'shockindex'%Shock Index
            parameterThresholdsToTest = 0.3:0.01:2;
        case 'dynea'%Dynamic Arterial Elastance
            parameterThresholdsToTest = 0.4:0.01:2.4;
        case 'dpdt'%maximum upslope of the arterial pressure waveform
            parameterThresholdsToTest = 150:5:1800;
        case 'svr'%Systemic Vascular Resistance
            parameterThresholdsToTest = 400:10:2600;
        case 'shockindexmap'%Shock Index MAP
            parameterThresholdsToTest = 0.5:0.01:2.5;
        %to-do make below thresholds match what we did for paper
        case 'deltamap65to75'
            parameterThresholdsToTest = [-50:-20 -19.5:0.5:-10 -9.8:0.2:10 10.5:0.5:20 21:50]; %-50:1:50;
        case 'deltamap75to85'
            parameterThresholdsToTest = [-50:-20 -19.5:0.5:-10 -9.8:0.2:10 10.5:0.5:20 21:50]; %-50:1:50;
        case 'deltamap85to95'
            parameterThresholdsToTest = [-50:-20 -19.5:0.5:-10 -9.8:0.2:10 10.5:0.5:20 21:50]; %-50:1:50;
        otherwise
            error('unrecognized parameter')
    end
else
    parameterThresholdsToTest = parameterThresholdsToTest;
end
if ~exist('excludeInterventions', 'var') || isempty(excludeInterventions)
    excludeInterventions = 1; %1 to exclude interventions, 0 to keep
end
if ~exist('alternateInterventionDefinitionFlag', 'var') || isempty(alternateInterventionDefinitionFlag)
    alternateInterventionDefinitionFlag = 0; %0 assumes intervention when MAP increased by more than 5mmHg within 20 seconds or when MAP increased more than 8mmHg within 2 minutes...1 uses 10mmHg and 10mmHg
end
if ~exist('keepDetailedTable', 'var') || isempty(keepDetailedTable)
    keepDetailedTable = 0; %flag where 1 means to save the idx and patientID,MAP,&HPI value of each TP,FP,etc.....makes it much slower
end


%--Outputs: struct named Outputs with the fields
%T_labelCountsPerPatPerThresh: table with threshold,patientid,tp count,fp count, ...
%T_labelCountsPerThresh: table threshold and sum of TP count,fpcount.... at each threshold from all the patients, also has sens,spec,ppv,npv computed from those sum of tpcount,fp count etc.
%T_AllPats: if keepDetailedTable flag is 1, then it is a table with the idx of all TP,FP,FN,TN, along with patientid and HPI and MAP value at those idx
%-------------------------------------------------------------------------

disp(['searchWindowTime=', num2str(searchWindowTime)]);

%preallocate variables to hold results per patient
cellArrOfT_perPatDetails = cell(numel(patientList),1);
cellArrOfT_perPatBasic = cell(numel(patientList),1);

%loop through each patient
parfor i_pts = 1:length(patientList)%parfor
    if rem(i_pts, 100) < 1
        fprintf('%d processed\n', i_pts);
    end

    patientID = patientList{i_pts};

    %---load the data for this patient and extract the parameter used to predict hypotension
    filename = fullfile(dataPath,[patientID '.mat']);
    if exist(filename,'file')
        data = load(filename);
    else
        filename = fullfile(dataPath,[patientID '_v37.mat']);
        data = load(filename);
    end

    % % Calculate delta MAP
    % m = apco.TR_MAP_disp;
    % mb = [repmat(m(1), mapDiffLen, 1); m(1:end-mapDiffLen)];
    % deltaMap = m - mb;
    data.TR_databad = data.TR_databad(:)';
    data.nocal_time = data.nocal_time(:)';
    data.TR_MAP_disp = data.TR_MAP_disp(:)';
    mapDiffLen = 9;
    mb = [repmat(data.TR_MAP_disp(1), 1, mapDiffLen), data.TR_MAP_disp(1:end-mapDiffLen)];
    data.TR_deltaMAP = data.TR_MAP_disp - mb;
    deltaMap = data.TR_deltaMAP(:)';
    switch lower(parameterName)
        case 'hpi'
            param2test = data.TR_HPI_disp(:)';
        case 'map'
            param2test = data.TR_MAP_disp(:)';
            %to-do make data match what we did for paper
        case 'co'
            param2test = data.TR_CO_disp(:)';
        case 'sv'
           param2test = data.TR_SVft(:)';
        case 'pulsepressure'
            param2test = data.TR_pulsepres(:)';
        case 'hr'
            param2test = data.TR_HR_disp(:)';
        case 'svv'
            param2test = data.TR_SVV_disp(:)';
        case 'shockindex'
            param2test = (data.TR_HR_disp(:)'./(data.TR_bp_sys(:)' + 1e-6));
        case 'dynea'
            param2test = data.TR_dynEa_disp(:)';
        case 'dpdt'
            param2test = data.TR_dpdt_disp(:)';
        case 'svr'
            param2test = 80.*(data.TR_MAP_disp(:)'-5) ./ (data.TR_CO_disp(:)'+10^-6);
        case 'shockindexmap'
            param2test = (data.TR_HR_disp(:)'./(data.TR_MAP_disp(:)' + 1e-6));
        case 'deltamap65to75'
            error('to-do')
            param2test = deltaMap(:)';
        case 'deltamap75to85'
            error('to-do')
            param2test = deltaMap(:)';
        case 'deltamap85to95'
            error('to-do')
            param2test = deltaMap(:)';
        otherwise
            error('unrecognized parameter');
    end

    %find the hypotension and intervention indices for this patient
    idxHypotension = findHypotensionIndices(data);
    idxIntervention = findInterventionIndices(data,alternateInterventionDefinitionFlag);


    lenOfPts = length(data.TR_HPI_disp);

    %preallocate variables to hold results per threshold for this patient
    if keepDetailedTable
        idxTested_AllThresholds = [];
        thresholdUsed_AllThresholds = [];
        labels_AllThresholds = [];
        lenDatabadSegment_AllThresholds = [];
    end
    [TP_perThresh,FP_perThresh,FN_perThresh,TN_perThresh] = deal(zeros(numel(parameterThresholdsToTest),1));
    for i_thresh=1:numel(parameterThresholdsToTest)
        paramThreshold = parameterThresholdsToTest(i_thresh);
        TP=0; FP=0; TN=0; FN=0;
        if keepDetailedTable
            label = strings(length(data.TR_HPI_disp),1);
            idxNotFullSearchWindow = max(1,(length(data.TR_HPI_disp)-(searchWindowTime*3) + 1)):length(data.TR_HPI_disp);
            label(idxNotFullSearchWindow) = "exclude_searchWindowEndOfCase";
            %keep track of end of case databad
            idxDatabad = find(data.TR_databad ~= 0);
            idxEndOfCaseAndDatabad = intersect(idxNotFullSearchWindow,idxDatabad);
            label(idxEndOfCaseAndDatabad) = "exclude_searchWindowEndOfCaseAndDatabad";
        end

        %loop through everypoint for this patient and classify alerts as TP,FP,FN,TN for this threshold
        for ii=1:length(data.TR_HPI_disp)-(searchWindowTime*3)

            idxToTest = ii;

            % To search next predTime minutes
            idxInSearchWindow = ii+1:min(ii+1+searchWindowTime*3-1,lenOfPts);

            %skip if databad
            if max(data.TR_databad(idxToTest))~=0
                if keepDetailedTable
                    label(ii) = "exclude_currentDatabad";
                end
                continue;
            end

            %skip if already hypotension
            %...using MAP definition of hypotension
            if any(ismember(idxHypotension,idxToTest))
                if keepDetailedTable
                    label(ii) = "exclude_currentHypotension";
                end
                continue;
            end

            if excludeInterventions
                %skip if already intervention
                %...using MAP definition of intervention
                if any(ismember(idxIntervention,idxToTest))
                    if keepDetailedTable
                        label(ii) = "exclude_currentIntervention";
                    end
                    continue;
                end
            end

            %
            %             % Skip if param is deltaMAP and all MAP are not in the MAP range
            %             if contains(lower(paramFlag), 'dmap')
            %                 if any(apco.TR_MAP_disp(idxToTest) <= mapRange(1)) || ...
            %                         any(apco.TR_MAP_disp(idxToTest) > mapRange(2))
            %                     continue;
            %                 end
            %             end

            %check if parameter is alerting
            if alertIfParamGreaterThanThresh == 1
                if min(param2test(idxToTest)) > paramThreshold % Alert Yes
                    paramIsAlerting = 1;
                else
                    paramIsAlerting = 0;
                end
            else
                if min(param2test(idxToTest)) < paramThreshold % Alert Yes
                    paramIsAlerting = 1;
                else
                    paramIsAlerting = 0;
                end
            end


            %apply the labeling logic to each point
            if paramIsAlerting % Alert Yes
                %check if there is hypotension in this search window
                if any(ismember(idxHypotension,idxInSearchWindow))
                    %Alert Yes,Event Yes
                    TP = TP+1;
                    if keepDetailedTable
                        label(ii) = "TP";
                    end
                else
                    %Alert Yes, Event No
                    %check if there were any interventions in the search window
                    if excludeInterventions && any(ismember(idxIntervention,idxInSearchWindow))
                        %Alert Yes, Event No, Intervention Yes
                        %exclude from analysis
                        if keepDetailedTable
                            label(ii) = "exclude_yesAlertNoHypButSearchWindowIntervention";
                        end
                    else
                        %Alert Yes, Event No, Intervention No
                        FP = FP + 1; %Alert, not event, no intervention, therefore FP
                        if keepDetailedTable
                            label(ii) = "FP";
                        end
                    end
                end
            else
                %Alert No

                %check if there is hypotension in this search window
                if any(ismember(idxHypotension,idxInSearchWindow))
                    %Alert No, Event Yes

                    %check if there was a hypotensive event that the parameter failed to predict at least 2 minutes in advance
                    idxFirstHypoInWindow = min(intersect(idxHypotension,idxInSearchWindow));
                    idxB4FirstHypoInWindow = (ii+1):(idxFirstHypoInWindow - 1);%if hypostarts at ii+1 this will be empty, thats ok

                    if alertIfParamGreaterThanThresh == 1
                        idxB4FirstHypoInWindow_WithParamAlerting = find(param2test(idxB4FirstHypoInWindow) > paramThreshold);
                    else
                        idxB4FirstHypoInWindow_WithParamAlerting = find(param2test(idxB4FirstHypoInWindow) < paramThreshold);
                    end
                    param_alertsLaterOnButStillBeforeFirstHypo = numel(idxB4FirstHypoInWindow_WithParamAlerting) >= 6;%(each point is 20s of data)

                    if ~param_alertsLaterOnButStillBeforeFirstHypo
                        %current Alert No, Event Yes, failed to predict at least 2 mins in advance
                        FN = FN+1; %Non-Alert, Event Yes, therefore FN
                        if keepDetailedTable
                            label(ii) = "FN";
                        end
                    end
                else
                    %Alert No, Event No

                    %check if there were any interventions in the search window
                    if excludeInterventions && any(ismember(idxIntervention,idxInSearchWindow))
                        %Alert No, Event No, Intervention Yes
                        %exclude from analysis
                        if keepDetailedTable
                            label(ii) = "exclude_noAlertButSearchWindowIntervention";
                        end
                    else
                        %Alert No, Event No, Intervention No
                        TN = TN+1; %Non-Alert, Non-Event, therfore TN
                        if keepDetailedTable
                            label(ii) = "TN";
                        end
                    end
                end
            end
        end%end loop through each data point
        %keep track of results for this patient, this threshold
        TP_perThresh(i_thresh) = TP;
        FP_perThresh(i_thresh) = FP;
        FN_perThresh(i_thresh) = FN;
        TN_perThresh(i_thresh) = TN;
        if keepDetailedTable
            labels_AllThresholds = [labels_AllThresholds;label];
            thresholdUsed_AllThresholds = [thresholdUsed_AllThresholds;paramThreshold.*ones(size(label))];
            idxTested_AllThresholds = [idxTested_AllThresholds;(1:numel(data.TR_MAP_disp))'];
        end
    end%end loop through each thresholds
    %keep track of results for this patient, for all thresholds
    cellArrOfT_perPatBasic{i_pts,1} = [parameterThresholdsToTest(:),TP_perThresh,FP_perThresh,FN_perThresh,TN_perThresh];
    if keepDetailedTable
        expersAll = [repelem(string(patientID),numel(labels_AllThresholds))]';
        MAP_All = data.TR_MAP_disp(idxTested_AllThresholds);MAP_All = MAP_All(:);
        HPI_All = data.TR_HPI_disp(idxTested_AllThresholds);HPI_All = HPI_All(:);
        T_thisPat = table(expersAll,idxTested_AllThresholds,thresholdUsed_AllThresholds,labels_AllThresholds,MAP_All,HPI_All);
        cellArrOfT_perPatDetails{i_pts} = T_thisPat;
    end
end%end loop through each patient
%data needed for stat
data =  vertcat(cellArrOfT_perPatBasic{:});
threshold = data(:,1);
TP_all = data(:,2);
FP_all = data(:,3);
FN_all = data(:,4);
TN_all = data(:,5);
tmp = cellfun(@(x) repmat(string(x),numel(parameterThresholdsToTest),1),patientList,'UniformOutput',false);
expersAnalyzed = vertcat(tmp{:});
%extra data for details
T_AllPats =  vertcat(cellArrOfT_perPatDetails{:});


%organize data into a table
T_labelCountsPerPatPerThresh = table(expersAnalyzed, threshold,TP_all,FP_all,FN_all,TN_all);
%get the total number of TP,FP,FN,TN for each threshold
T_labelCountsPerThresh = varfun(@sum,T_labelCountsPerPatPerThresh,'GroupingVariables','threshold','InputVariables',{'TP_all','FP_all','FN_all','TN_all'});
%compute Sens,Spec,PPV,NPV at each threshold
T_labelCountsPerThresh.Sens = T_labelCountsPerThresh.sum_TP_all ./ (T_labelCountsPerThresh.sum_TP_all + T_labelCountsPerThresh.sum_FN_all);
T_labelCountsPerThresh.Spec = T_labelCountsPerThresh.sum_TN_all ./ (T_labelCountsPerThresh.sum_TN_all + T_labelCountsPerThresh.sum_FP_all);
T_labelCountsPerThresh.PPV = T_labelCountsPerThresh.sum_TP_all ./ (T_labelCountsPerThresh.sum_TP_all + T_labelCountsPerThresh.sum_FP_all);
T_labelCountsPerThresh.NPV = T_labelCountsPerThresh.sum_TN_all ./ (T_labelCountsPerThresh.sum_TN_all + T_labelCountsPerThresh.sum_FN_all);


Outputs = struct();
Outputs.T_labelCountsPerPatPerThresh = T_labelCountsPerPatPerThresh;
Outputs.T_labelCountsPerThresh = T_labelCountsPerThresh;
Outputs.T_AllPats = T_AllPats;
end

%%

function idxHypotension = findHypotensionIndices(data)
%hypotension defined as MAP < 65 for at least a minute (each data point represents 20 seconds of data)
hypotensionThreshold= 64.5;

%find index where map < threshold
idxBelowThreshold = [];
for i=1:length(data.TR_databad)-2
    if min(data.TR_MAP_disp(i:i+2) < hypotensionThreshold)==1
        idxBelowThreshold = [idxBelowThreshold i:i+2];
    end
end
idxBelowThreshold = unique(idxBelowThreshold);

%find the start and stop index of each period with MAP < threshold
if isempty(idxBelowThreshold)
    strtidx = [];
    stopidx = [];
else
    diff_idxBelowThreshold = diff(idxBelowThreshold);
    idx = find(diff_idxBelowThreshold>1);
    if isempty(idx)
        strtidx = idxBelowThreshold(1);
        stopidx = idxBelowThreshold(end);
    else
        strtidx = [idxBelowThreshold(1) idxBelowThreshold(idx+1)];
        stopidx = [idxBelowThreshold(idx) idxBelowThreshold(end)];
    end
end

%hypotension starts when MAP < threshold for at least a minute
strtidx = strtidx + 2;
idxHypotension = [];
for i=1:length(strtidx)
    idxHypotension = [idxHypotension strtidx(i):stopidx(i)];
end

%exclude poor arterial waveform signals (e.g., line flushing, significantly damped waveforms, and other waveform artifacts)
idxbad = find(data.TR_databad>=1);
idxHypotension = setdiff(idxHypotension, idxbad);
end

function idxInterventions = findInterventionIndices(data,alternateDefinitionFlag)
%Treatment interventions were assumed when MAP increased more than 5 mmHg within 20 seconds (most probably caused by vasopressor or inotropic injections)
% or when MAP increased more than 8 mmHg within 2 minutes (change in vasopressor or inotropic infusion rate or fluid bolus) when the MAP was <75 mmHg.

% As a sensitivity analysis for the effects of excluding presumed treatment periods, 
% we also conducted a forward analysis with a threshold of 10 mmHg MAP change to define an intervention
if ~alternateDefinitionFlag
    fastThr = 5;
    slowThr = 8;
else
    fastThr = 10;
    slowThr = 10;
end

fastRiseFlag = zeros(1,length(data.TR_MAP_disp));
for ii= 1:length(data.TR_MAP_disp)-5
    if data.TR_MAP_disp(ii+1)-data.TR_MAP_disp(ii)>fastThr && data.TR_MAP_disp(ii)<75.5 && data.TR_MAP_disp(ii+2)-data.TR_MAP_disp(ii+1)>-1
        fastRiseFlag(ii) = 1;
    end
end
slowRiseFlag = zeros(1,length(data.TR_MAP_disp));
for ii= 1:length(data.TR_MAP_disp)-5
    if max(data.TR_MAP_disp(ii:ii+5))-data.TR_MAP_disp(ii)>slowThr && data.TR_MAP_disp(ii)<75.5
        slowRiseFlag(ii) = 1;
    end
end

idxInterventions = find(fastRiseFlag | slowRiseFlag);
end