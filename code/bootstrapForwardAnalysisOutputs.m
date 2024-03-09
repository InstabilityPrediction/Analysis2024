function Outputs = bootstrapForwardAnalysisOutputs(T_labelCountsPerPatPerThresh,numBoostrapIterations,randomSeed,plotMode)
%example how to use:
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
% 
% %---set inputs for bootstrapping for confidence interval calculations
% numBoostrapIterations = 2000;
% randomSeed = 1;
% plotMode = 0;
% 
% 
% %--------------------------------------------------------------------------
% %---label each point as TP,FP,FN,TN or excluded
% OutputsLabels = labelData_forwards(patientList,dataPath,...
%     searchWindowTime,parameterName,parameterThresholdsToTest,alertIfParamGreaterThanThresh,...
%     excludeInterventions,alternateInterventionDefinitionFlag,...
%     keepDetailedTable);
% 
% %---use bootstrapping to account for repeated measurements from each subject for the calculation of confidence intervals
% OutputsBootstrap = bootstrapForwardAnalysisOutputs(OutputsLabels.T_labelCountsPerPatPerThresh,numBoostrapIterations,randomSeed,plotMode);
% 
%--------------------------------------------------------------------------
%--mandatory inputs
%T_labelCountsPerPatPerThresh: table
%numBoostrapIterations: 2000

%--handle optional inputs
if ~exist('randomSeed', 'var') || isempty(randomSeed)
    randomSeed = 1; % used so that results can be replicated since there is random replacement
end
if ~exist('plotMode', 'var') || isempty(plotMode)
    plotMode = 0; % 0 makes no plot, 1 makes an ROC plot with 95% Confidence interval lines as well
end

%set random seed for reproducibility
rng(randomSeed);


useCategorical = 1;
if useCategorical
    %for faster strcmp convert string to categorical
    T_labelCountsPerPatPerThresh.expersCategorical = categorical(T_labelCountsPerPatPerThresh.expersAnalyzed);
    expersForBootStrap = unique(T_labelCountsPerPatPerThresh.expersCategorical);
else
    expersForBootStrap = unique(T_labelCountsPerPatPerThresh.expersAnalyzed);
end

%preallocate containers to hold the bootstrapping results
cellArrOfTStatsPerBootStrap = cell(numBoostrapIterations,1);
resultsAtSpecificThresholdsPerBootstrap = struct();
parfor i=1:numBoostrapIterations%parfor
    rng('default'); rng(i);
    if mod(i,100) == 0
        disp(i);
    end
    %select n random patients with replacement
    expersThisBootStrap = datasample(expersForBootStrap, numel(expersForBootStrap));
    %get the data for these patients(some may be included multiple times)

    cellArrOfIdxToGrabPerPat_thisBootstrap = cell(numel(expersThisBootStrap),1);
    for j=1:numel(expersThisBootStrap)
        if useCategorical
            expnameTmp = expersThisBootStrap(j);
            idxToGrabThisPat = find(T_labelCountsPerPatPerThresh.expersCategorical==expnameTmp);
        else
            expnameTmp = expersThisBootStrap{j};
            idxToGrabThisPat = find(strcmp(T_labelCountsPerPatPerThresh.expersAnalyzed,expnameTmp));
        end
        cellArrOfIdxToGrabPerPat_thisBootstrap{j} = idxToGrabThisPat;
    end
    idxRowsToGrab_thisBootstrap = vertcat(cellArrOfIdxToGrabPerPat_thisBootstrap{:});
    T_labelCountsPerPatPerThresh_thisBootStrap = T_labelCountsPerPatPerThresh(idxRowsToGrab_thisBootstrap,:);
    
    %compute stats for this group of patients at each threshold...will be used for ROC curve
    %TP,FP,FN,TN
    T_labelCountsPerThresh_thisBootStrap = varfun(@sum,T_labelCountsPerPatPerThresh_thisBootStrap,'GroupingVariables','threshold','InputVariables',{'TP_all','FP_all','FN_all','TN_all'});
    %Sens,Spec,PPV,NPV
    T_labelCountsPerThresh_thisBootStrap.Sens = T_labelCountsPerThresh_thisBootStrap.sum_TP_all ./ (T_labelCountsPerThresh_thisBootStrap.sum_TP_all + T_labelCountsPerThresh_thisBootStrap.sum_FN_all);
    T_labelCountsPerThresh_thisBootStrap.Spec = T_labelCountsPerThresh_thisBootStrap.sum_TN_all ./ (T_labelCountsPerThresh_thisBootStrap.sum_TN_all + T_labelCountsPerThresh_thisBootStrap.sum_FP_all);
    T_labelCountsPerThresh_thisBootStrap.PPV = T_labelCountsPerThresh_thisBootStrap.sum_TP_all ./ (T_labelCountsPerThresh_thisBootStrap.sum_TP_all + T_labelCountsPerThresh_thisBootStrap.sum_FP_all);
    T_labelCountsPerThresh_thisBootStrap.NPV = T_labelCountsPerThresh_thisBootStrap.sum_TN_all ./ (T_labelCountsPerThresh_thisBootStrap.sum_TN_all + T_labelCountsPerThresh_thisBootStrap.sum_FN_all);
    T_labelCountsPerThresh_thisBootStrap.Recall = T_labelCountsPerThresh_thisBootStrap.Sens;%recall and sensitvity are the same
    T_labelCountsPerThresh_thisBootStrap.Precision = T_labelCountsPerThresh_thisBootStrap.PPV;%precision and PPV are the same
    T_labelCountsPerThresh_thisBootStrap.bootStrapIteration(:,1) = i;
    
    %compute AUC and stats at specific thresholds that can change each iteration...these stats will be used for PPT slides
    AUC = trapz(T_labelCountsPerThresh_thisBootStrap.Spec, T_labelCountsPerThresh_thisBootStrap.Sens);%we typically do 1-Spec vs Sens for ROC curve but doing trapz with 1-Spec,Sens gives negative number
    % Modification to change AUC if negative and if AUC < 0.5
    if AUC < 0
        AUC = abs(AUC);
    end
    if AUC < 0.5
        AUC = 1-AUC;
    end
    resultsAtSpecificThresholdsPerBootstrap(i).AUC = AUC;
    
    %youden and balanced thresh
    [~, idxYouden] = max(T_labelCountsPerThresh_thisBootStrap.Sens + T_labelCountsPerThresh_thisBootStrap.Spec - 1);
    [~, idxBalanced] = min(abs(T_labelCountsPerThresh_thisBootStrap.Spec - T_labelCountsPerThresh_thisBootStrap.Sens));
    resultsAtSpecificThresholdsPerBootstrap(i).threshold_Youden = T_labelCountsPerThresh_thisBootStrap.threshold(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(i).Sens_Youden = T_labelCountsPerThresh_thisBootStrap.Sens(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(i).Spec_Youden = T_labelCountsPerThresh_thisBootStrap.Spec(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(i).PPV_Youden = T_labelCountsPerThresh_thisBootStrap.PPV(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(i).NPV_Youden = T_labelCountsPerThresh_thisBootStrap.NPV(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(i).threshold_Balanced = T_labelCountsPerThresh_thisBootStrap.threshold(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(i).Sens_Balanced = T_labelCountsPerThresh_thisBootStrap.Sens(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(i).Spec_Balanced = T_labelCountsPerThresh_thisBootStrap.Spec(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(i).PPV_Balanced = T_labelCountsPerThresh_thisBootStrap.PPV(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(i).NPV_Balanced = T_labelCountsPerThresh_thisBootStrap.NPV(idxBalanced);
    
    cellArrOfTStatsPerBootStrap{i} = T_labelCountsPerThresh_thisBootStrap;
end
T_ResultsAllThresholdsAllBootstraps = vertcat(cellArrOfTStatsPerBootStrap{:});
T_ResultsAtSpecificThresholdsAllBootstraps = struct2table(resultsAtSpecificThresholdsPerBootstrap);

%compute stats for each threshold
bootStrapStatsAllThresholds = grpstats(T_ResultsAllThresholdsAllBootstraps,"threshold",{@mean,@std,@median,@(x) quantile(x,0.025),@(x) quantile(x,0.975)},"DataVars",["Sens","Spec","PPV","NPV"]);
varNamesTemp = bootStrapStatsAllThresholds.Properties.VariableNames;
bootStrapStatsAllThresholds = renamevars(bootStrapStatsAllThresholds,varNamesTemp(startsWith(varNamesTemp,'Fun4')),strrep(varNamesTemp(startsWith(varNamesTemp,'Fun4')),'Fun4','prctile2p5'));
bootStrapStatsAllThresholds = renamevars(bootStrapStatsAllThresholds,varNamesTemp(startsWith(varNamesTemp,'Fun5')),strrep(varNamesTemp(startsWith(varNamesTemp,'Fun5')),'Fun5','prctile97p5'));
%compute stats for AUC and specific thresholds
T_ResultsAtSpecificThresholdsAllBootstraps.dummyVar(:,1) = "sameGroup";
bootStrapStatsSpecificThresholds = grpstats(T_ResultsAtSpecificThresholdsAllBootstraps,"dummyVar",{@mean,@std,@median,@(x) quantile(x,0.025),@(x) quantile(x,0.975)});
varNamesTemp = bootStrapStatsSpecificThresholds.Properties.VariableNames;
bootStrapStatsSpecificThresholds = renamevars(bootStrapStatsSpecificThresholds,varNamesTemp(startsWith(varNamesTemp,'Fun4')),strrep(varNamesTemp(startsWith(varNamesTemp,'Fun4')),'Fun4','prctile2p5'));
bootStrapStatsSpecificThresholds = renamevars(bootStrapStatsSpecificThresholds,varNamesTemp(startsWith(varNamesTemp,'Fun5')),strrep(varNamesTemp(startsWith(varNamesTemp,'Fun5')),'Fun5','prctile97p5'));


if plotMode
    % Plot ROC curve: 1-Specificity versus Sensitivity (median with 95% CI)
    figure; hold on; box on;
    plot(1-(bootStrapStatsAllThresholds.median_Spec), bootStrapStatsAllThresholds.median_Sens, '-b', 'LineWidth', 4);
    plot(1-(bootStrapStatsAllThresholds.prctile2p5_Spec), bootStrapStatsAllThresholds.prctile2p5_Sens, '-b', 'LineWidth', 1);
    plot(1-(bootStrapStatsAllThresholds.prctile97p5_Spec), bootStrapStatsAllThresholds.prctile97p5_Sens, '-b', 'LineWidth', 1);
    plot(0:0.1:1, 0:0.1:1, 'm', 'LineWidth', 2);
    legend('HPI');
    xlabel('1-Specificity', 'FontSize', 15); ylabel('Sensitivity', 'FontSize', 15);
end

Outputs = struct();
Outputs.bootStrapStatsAllThresholds = bootStrapStatsAllThresholds;
Outputs.bootStrapStatsSpecificThresholds = bootStrapStatsSpecificThresholds;
Outputs.T_ResultsAllThresholdsAllBootstraps = T_ResultsAllThresholdsAllBootstraps;
Outputs.T_ResultsAtSpecificThresholdsAllBootstraps = T_ResultsAtSpecificThresholdsAllBootstraps;
end
