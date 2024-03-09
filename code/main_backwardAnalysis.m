% This software code is to perform case-control (backwards) analysis. 
%
% Please review LICENSE.md before you use this code. 

%--------------------------------------------------------------------------
% To run this software code:
% 1) "dataPath" is the directory where all the .mat files are located
% 2) patientList" is an array of the names of the .mat files intended to be used in the analysis
% 3) each patient has one .mat file, which is nx1 arrays with the following variable names
%     HPI(Hypotension Prediction Index) = 'TR_HPI_disp'
%     MAP(Mean Arterial Pressure) = 'TR_MAP_disp'
%     CO(Cardiac Ouput) = 'TR_CO_disp'
%     SV(Stroke Volume) = 'TR_SVft'
%     PulsePressure = 'TR_pulsepres'
%     HR(Heart Rate) = 'TR_HR_disp'
%     SVV(Stroke Volume Variation) = 'TR_SVV_disp'
%     Eadyn(Dynamic arterial elsatance) = 'TR_dynEa_disp'
%     Systolic pressure = 'TR_bp_sys'
%     dP/dt             = 'TR_dpdt_disp'
%     data quality indicator = 'TR_databad' % 1 - bad quality; 0 - good quality
%     time = 'nocal_time' %time in minutes, such as 2:20:30PM = 860.50 (=14*60+20+30/60)


clear;

%---traceability
%get name of this script
ScriptToMakeThis = mfilename('fullpath');
TimeScriptRan = datetime('now');
matlabVersion=version;

TraceabilityInfo.ScriptToMakeThis = ScriptToMakeThis;
TraceabilityInfo.TimeScriptRan = TimeScriptRan;
TraceabilityInfo.matlabVersion = matlabVersion;
%--------------------------------------------------------------------------

saveMode = 1;
savePath = '..\results';

%---set inputs for labeling
%which patients data to consider
[masterPatientList] = helpers.GetListOfPatients();
patientList = masterPatientList.Total_FT;
patientGroupName = 'FT';%used to name the results file

dataPath = '..\data'; %where to load the patient data files
%time period before to predict hypotension
timeBeforeHypotensionOptions = [5,10,15];%,5,10,15 minutes used in manuscript

%which parameter to evaluate to predict hypotension
parameterName = 'HPI';%candidates are: HPI,MAP,CO,SV,PulsePressure,HR,SVV,ShockIndex,dynEa,ShockIndexMap,deltaMap65to75,deltaMap75to85,deltaMap85to95

%---set inputs for bootstrapping
numBootstrapIterations = 2000;


%--------------------------------------------------------------------------


%---label certain data points as positives or negatives and get the corresponding value of the desired parameter
OutputsLabels = labelData_backwards(patientList,dataPath,...
    parameterName,timeBeforeHypotensionOptions);

%compare the values to the labels(positives/negatives) to get the ROC statistics for each of the time periods before hypotension
T_idxAndValues_neg = OutputsLabels.negatives.T_idxAndValues;
cellArrOfT_stats = cell(numel(OutputsLabels.positives),1);
allBootstrapOutputs = struct();
for i=1:numel(OutputsLabels.positives)
    timeBeforeHypotension = OutputsLabels.positives(i).timeBeforeHypotension;
    %get the data for this time period
    T_idxAndValues_pos = OutputsLabels.positives(i).T_idxAndValues;
    parameterValues = [T_idxAndValues_pos.parameterValues;T_idxAndValues_neg.parameterValues];
    labels = [ones(height(T_idxAndValues_pos),1);zeros(height(T_idxAndValues_neg),1)];
    patientIDs = [T_idxAndValues_pos.patientID;T_idxAndValues_neg.patientID;];

    %---use bootstrapping to account for repeated measurements from each subject for the calculation of confidence intervals
    OutputsBootstrap = bootstrapBackwardAnalysisOutputs(labels,parameterValues,patientIDs,numBootstrapIterations);
    allBootstrapOutputs(i).timeBeforeHypotension = timeBeforeHypotension;
    allBootstrapOutputs(i).OutputsBootstrap = OutputsBootstrap;

    %---organize the stats
    statsForPPT_ROC = helpers.formatBackwardResultsForPPT(OutputsBootstrap.bootStrapStatsSpecificThresholds,labels);
    [statsForPPT_ROC.timeBeforeHypotension] = deal(timeBeforeHypotension);
    cellArrOfT_stats{i} = struct2table(statsForPPT_ROC);
end
T_stats_ROC = vertcat(cellArrOfT_stats{:});

%---save the data, filename is based on input options
if saveMode
    Inputs = struct();
    Inputs.patientList = patientList;
    Inputs.dataPath = dataPath;
    Inputs.parameterName = parameterName;
    Inputs.timeBeforeHypotensionOptions = timeBeforeHypotensionOptions;
    Inputs.numBoostrapIterations = numBootstrapIterations;
    Outputs = struct();
    Outputs.OutputsLabels = OutputsLabels;
    Outputs.allBootstrapOutputs = allBootstrapOutputs;
    Outputs.T_stats_ROC = T_stats_ROC;
    saveFilename = fullfile(savePath,['backwards_' parameterName '_' patientGroupName '.mat']);
    save(saveFilename,'Inputs','Outputs','TraceabilityInfo')
end