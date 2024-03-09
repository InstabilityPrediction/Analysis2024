% This software code is to perform cohort (forwards) analysis. 
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
searchWindowTime = 10; %candidates are: 5,10,15.  Length in minutes of search window for forward analysis

%which parameter to evaluate to predict hypotension
parameterName = 'HPI'; %candidates are: HPI,MAP,CO,SV,PulsePressure,HR,SVV,ShockIndex,dynEa,ShockIndexMap,deltaMap65to75,deltaMap75to85,deltaMap85to95

parameterThresholdsToTest = [];
alertIfParamGreaterThanThresh = 0;%1 means Parameter > thresh is an alert, 0 means Parameter < thresh is an alert

%interventions
excludeInterventions = 1; % 1 to exclude interventions, 0 to keep
alternateInterventionDefinitionFlag = 0; %0: intervention is when MAP increased by >5mmHg within 20 seconds or when MAP increased >8mmHg within 2 minutes; 1: uses 10mmHg and 10mmHg

keepDetailedTable = 0;  % 1 to save the idx and patientID,MAP,&HPI value of each TP,FP,etc.....makes it much slower

%---set inputs for bootstrapping for confidence interval calculations
numBootstrapIterations = 2000;
randomSeed = 1;
plotMode = 0;


%--------------------------------------------------------------------------
%---label each point as TP,FP,FN,TN or excluded
OutputsLabels = labelData_forwards(patientList,dataPath,...
    searchWindowTime,parameterName,parameterThresholdsToTest,alertIfParamGreaterThanThresh,...
    excludeInterventions,alternateInterventionDefinitionFlag,...
    keepDetailedTable);

%---use bootstrapping to account for repeated measurements from each subject for the calculation of confidence intervals
OutputsBootstrap = bootstrapForwardAnalysisOutputs(OutputsLabels.T_labelCountsPerPatPerThresh,numBootstrapIterations,randomSeed,plotMode);

%---organize the stats
statsForPPT_ROC = helpers.formatForwardResultsForPPT(OutputsLabels.T_labelCountsPerThresh,OutputsBootstrap.bootStrapStatsSpecificThresholds);

%---save the data, filename is based on input options
if saveMode
    if excludeInterventions
        if alternateInterventionDefinitionFlag
            intvString = 'yesExcludeInt10';
        else
            intvString = 'yesExcludeInt';
        end
    else
        intvString = 'noExcludeInt';
    end
    saveFilename = fullfile(savePath,['test_' parameterName '_' num2str(searchWindowTime) 'min_' patientGroupName '_' intvString '.mat']);
end