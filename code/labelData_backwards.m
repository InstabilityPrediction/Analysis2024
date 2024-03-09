function Outputs = labelData_backwards(patientList,dataPath,...
    parameterName,timeBeforeHypotensionOptions)


%%example how to use:
% %---set inputs for labeling
% %which patients data to consider
% [masterPatientList] = helpers.GetListOfPatients();
% patientList = masterPatientList.Total_FT;
% patientGroupName = 'FT';%used to name the results file
% 
% dataPath = '..\data'; %where to load the patient data files
% %time period before to predict hypotension
% timeBeforeHypotensionOptions = [5,10,15];%,5,10,15 minutes used in manuscript
% 
% %which parameter to evaluate to predict hypotension
% parameterName = 'HPI';%candidates are: HPI,MAP,CO,SV,PulsePressure,HR,SVV,ShockIndex,dynEa,ShockIndexMap,deltaMap65to75,deltaMap75to85,deltaMap85to95
% %---set inputs for bootstrapping
% numBootstrapIterations = 2000;
% 
% 
% %--------------------------------------------------------------------------
% %---label certain data points as positives or negatives and get the corresponding value of the desired parameter
% OutputsLabels = labelData_backwards(patientList,dataPath,...
%     parameterName,timeBeforeHypotensionOptions);
%--------------------------------------------------------------------------
%--mandatory inputs
%patientList: nx1 cell array, filenames of patient data (without the .mat file extension)
%dataPath: char, file path containing the patient data files
%parameterName: name of the variable to predict hypotension
%searchWindowTime: length in minutes of search window for forward analysis


%--Outputs: struct named Outputs with the fields
%positives: struct array where each row is table with idx,parameterValues,patientid for a different time before hypotension
%negatives: table with idx,parameterValues,patientid
%-------------------------------------------------------------------------

%preallocate variables to hold results per patient
numTimeOptions = numel(timeBeforeHypotensionOptions);
positivesPerPat = cell(numel(patientList),numTimeOptions);
negativesPerPat = cell(numel(patientList),numTimeOptions);

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
    data.TR_databad = data.TR_databad(:)';
    data.TR_map     = data.TR_MAP_disp(:)';
    switch lower(parameterName)
        case 'hpi'
            param2test = data.TR_HPI_disp(:)';
        case 'map'
            param2test = data.TR_MAP_disp(:)';
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
            % param2test = deltaMap(:)';
        case 'deltamap75to85'
            error('to-do')
            % param2test = deltaMap(:)';
        case 'deltamap85to95'
            error('to-do')
            % param2test = deltaMap(:)';
        otherwise
            error('unrecognized parameter');
    end

    %exclude poor arterial waveform signals (e.g., line flushing, significantly damped waveforms, and other waveform artifacts)
    idxbad = find(data.TR_databad>=1);
    idxgood = setdiff(1:length(data.TR_databad), idxbad);
    if isempty(idxgood)
        continue;
    end
    
    %when there is poor arterial waveform signals, use the previous good value
    for ii=1:length(idxbad)
        data.TR_map(idxbad(ii))    = data.TR_map(max(1,idxbad(ii)-1));
    end
    data.TR_map = round(data.TR_map);
    
    %find the start and stop of hypotension events for this patient
    [idxHypoStarts,idxHypoStops] = findHypotensionStartAndStopIndices(data);

    %---find positive data points(5,10,15 mins before the hypotensive event)
    for i_timeOptions=1:numTimeOptions
        timeBeforeHypotension = timeBeforeHypotensionOptions(i_timeOptions);
        idxPositives_tmp = findPositiveDataPoints(data,idxHypoStarts,idxHypoStops,timeBeforeHypotension,idxbad);
        parameterValuesPositives_tmp = param2test(idxPositives_tmp);
        T_pos = table();
        T_pos.idx = idxPositives_tmp(:);
        T_pos.parameterValues = parameterValuesPositives_tmp(:);
        T_pos.patientID(:,1) = string(patientID);
        positivesPerPat{i_pts,i_timeOptions} = T_pos;
    end
   
    %--find negative data points(mid point of 30-min continous section of data points at least 20 min apart from any hypotensive event and with a MAP > 75 mmHg
    idxNegatives = findNegativeDataPoints(data,idxHypoStarts,idxHypoStops,idxbad);
    parameterValuesNegatives = param2test(idxNegatives);
    T_neg = table();
    T_neg.idx = idxNegatives(:);
    T_neg.parameterValues = parameterValuesNegatives(:);
    T_neg.patientID(:,1) = string(patientID);
    negativesPerPat{i_pts,1} = T_neg;
end
T_idxAndValues_neg = vertcat(negativesPerPat{:});
%organize data for output
Outputs = struct();
Outputs.negatives.T_idxAndValues = T_idxAndValues_neg;
for i_timeOptions=1:numel(timeBeforeHypotensionOptions)
    timeBeforeHypotension = timeBeforeHypotensionOptions(i_timeOptions);
    T_idxAndValues_pos = vertcat(positivesPerPat{:,i_timeOptions});
    Outputs.positives(i_timeOptions).timeBeforeHypotension = timeBeforeHypotension;
    Outputs.positives(i_timeOptions).T_idxAndValues = T_idxAndValues_pos;
end

end

function [strtidx,stopidx] = findHypotensionStartAndStopIndices(data)
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
function idx_x_minsB4HypoStart = findPositiveDataPoints(data,idxHypoStarts,idxHypoStops,timeBeforeHypotension,idxbad)
singleEventDataFlag = 'Single';
predThreshHi = 1000; %[1000,100,90,80,75]; %Count only if map<predThresh in -5,-10,-15 minutes
predThreshLo = 0;
 diffThresh   = -10; % Exclude those sudden drop of pressure whose change is >= this threshold
tmpIndx = [];
for ii=1:length(idxHypoStarts)
    %get data point x mins before hypotension (add little time buffer to account for non-perfect 20s increments)
    tmpstrt = find(data.nocal_time>=data.nocal_time(idxHypoStarts(ii))-timeBeforeHypotension-0.16 & data.nocal_time<=data.nocal_time(idxHypoStarts(ii))-timeBeforeHypotension+0.15);
    tmpstop = find(data.nocal_time>=data.nocal_time(idxHypoStarts(ii))-0.3333-0.16 & data.nocal_time<=data.nocal_time(idxHypoStarts(ii))-0.3333+0.15);
    idxEval = tmpstrt:tmpstop;
    if strcmp(singleEventDataFlag, 'Single')
        idxEval = tmpstrt;
    end
    if isempty(idxEval) || isempty(tmpstrt:tmpstop) || (min(idxEval)<=0 && partialSegmentFlag==0)
        continue;
    end
    %exclude points where MAP drop between point and start of hypotension was not physiological
    if min(data.TR_map(idxEval)<predThreshHi)==1 && min(data.TR_map(idxEval)>predThreshLo)==1 && min(diff(data.TR_map(tmpstrt:tmpstop)))>=diffThresh
        if ii==1 || (ii>1 && (idxEval(1)>idxHypoStops(ii-1)))
            tmpIndx = [tmpIndx idxEval];
        end
    end
end
%exclude if point x mins before hypotension occured when there was poor arterial waveform signal
tmpIndx = setdiff(tmpIndx, idxbad);
%exclude if point occured before start of monitoring
tmpIndx = tmpIndx(tmpIndx>0);
idx_x_minsB4HypoStart = tmpIndx;
end
function idxNegatives = findNegativeDataPoints(data,idxHypoStarts,idxHypoStops,idxbad)
noDatabadInIndex = 1;
singleStableDataFlag = 'Single';
    %% index for stable, >=75mmHg
        stableThresh=75;
        timeMin=20; tmpIndx = [];

        %find segments of data that are not hypotension(outside 'grey zone')
        for ii=1:length(idxHypoStarts)
            if ii==1
                tmpstrt = 1;
                tmpstop = find(data.nocal_time<=data.nocal_time(idxHypoStarts(ii))-timeMin-0.3333+0.15, 1, 'last' );
                idxEval = tmpstrt:tmpstop;
                tmpIndx = [tmpIndx idxEval(data.TR_map(idxEval)>=stableThresh)];
            end
            if ii==length(idxHypoStarts)
                tmpstrt = find(data.nocal_time>=data.nocal_time(idxHypoStops(end))+timeMin+0.3333-0.16, 1 );
                tmpstop = length(data.TR_bp_sys);
                idxEval = tmpstrt:tmpstop;
            else
                tmpstrt = find(data.nocal_time>=data.nocal_time(idxHypoStops(ii))+timeMin+0.3333-0.16, 1 );
                tmpstop = find(data.nocal_time<=data.nocal_time(idxHypoStarts(ii+1))-timeMin-0.3333+0.15, 1, 'last' );
                idxEval = tmpstrt:tmpstop;
            end
            tmpIndx = [tmpIndx idxEval(data.TR_map(idxEval)>=stableThresh)];
        end
        if isempty(idxHypoStarts)
            idxEval = 1:length(data.TR_map);
            tmpIndx = [tmpIndx idxEval(data.TR_map(idxEval)>=stableThresh)];
        end
        %exclude if points occured when there was poor arterial waveform signal
        if noDatabadInIndex==1
            tmpIndx = setdiff(tmpIndx, idxbad);
        end
        %exlcude if point occured before start of monitoring
        tmpIndx = tmpIndx(tmpIndx>0);

        %find the midpoint of 30 minute segments
        if ~isempty(tmpIndx) && strcmp(singleStableDataFlag, 'Single')
            d      = find(diff([tmpIndx 0]) ~= 1);
            strInd = [tmpIndx(1) tmpIndx(d(1:end-1)+1)];
            endInd = tmpIndx(d);
            temp = endInd - strInd;
            idx1 = find(temp>=90);
            idx2 = [];
            for jj = 1:length(idx1)
                len = floor(temp(idx1(jj))/90);
                for kk = 1:len
                    idx2 = [idx2 strInd(idx1(jj))+44+(kk-1)*90];
                end
            end
            tmpIndx = idx2;
            if noDatabadInIndex==1
                tmpIndx = setdiff(tmpIndx, idxbad);
            end
            tmpIndx = tmpIndx(tmpIndx>0);
        end
        idxNegatives = tmpIndx;
end