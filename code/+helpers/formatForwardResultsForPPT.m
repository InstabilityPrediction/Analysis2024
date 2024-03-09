function statsForPPT = formatForwardResultsForPPT(T_labelCountsPerThresh,bootStrapStatsSpecificThresholds)

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
% %---organize the stats
% formatStatsForExcel = 1;
% statsForPPT_ROC = helpers.formatForwardResultsForPPT(OutputsLabels.T_labelCountsPerThresh,OutputsBootstrap.bootStrapStatsSpecificThresholds);


%--------------------------------------------------------------------------
%--mandatory inputs
%T_labelCountsPerThresh = table of nominal results (aka original patients) from labelData_forwards.m
%bootStrapStatsSpecificThresholds = table of bootstrap stats, from bootstrapForwardAnalysisOutputs.m

%--outputs
%formatted results for PPT
%columns: AUC, Sens,Spec,PPV,NPV, Threshold
%rows: Youden, Balanced
%each cell has the boostrap results aka median [2.5%tile, 97.5%tile], except the threshold rows just have the nominal values for the original patients


    statsForPPT = struct('AUC',[],'Sens',[],'Spec',[],'PPV',[],'NPV',[],'Threshold',[]);

    [~, idxYouden] = max(T_labelCountsPerThresh.Sens + T_labelCountsPerThresh.Spec - 1);
    [~, idxBalanced] = min(abs(T_labelCountsPerThresh.Spec - T_labelCountsPerThresh.Sens));
    threshold_Youden = T_labelCountsPerThresh.threshold(idxYouden);
    threshold_Balanced = T_labelCountsPerThresh.threshold(idxBalanced);
    %Youden first
    statsForPPT(1).AUC = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_AUC,...
        bootStrapStatsSpecificThresholds.prctile2p5_AUC,bootStrapStatsSpecificThresholds.prctile97p5_AUC);
    statsForPPT(1).Sens = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_Sens_Youden,...
        bootStrapStatsSpecificThresholds.prctile2p5_Sens_Youden,bootStrapStatsSpecificThresholds.prctile97p5_Sens_Youden);
    statsForPPT(1).Spec = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_Spec_Youden,...
        bootStrapStatsSpecificThresholds.prctile2p5_Spec_Youden,bootStrapStatsSpecificThresholds.prctile97p5_Spec_Youden);
    statsForPPT(1).PPV = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_PPV_Youden,...
        bootStrapStatsSpecificThresholds.prctile2p5_PPV_Youden,bootStrapStatsSpecificThresholds.prctile97p5_PPV_Youden);
    statsForPPT(1).NPV = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_NPV_Youden,...
        bootStrapStatsSpecificThresholds.prctile2p5_NPV_Youden,bootStrapStatsSpecificThresholds.prctile97p5_NPV_Youden);
    %get threshold from nominal values
    statsForPPT(1).Threshold = sprintf('Youden = %d',threshold_Youden);
    
    %Balanced second, AUC is the same
    statsForPPT(2).AUC = statsForPPT(1).AUC;
    statsForPPT(2).Sens = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_Sens_Balanced,...
        bootStrapStatsSpecificThresholds.prctile2p5_Sens_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_Sens_Balanced);
    statsForPPT(2).Spec = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_Spec_Balanced,...
        bootStrapStatsSpecificThresholds.prctile2p5_Spec_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_Spec_Balanced);
    statsForPPT(2).PPV = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_PPV_Balanced,...
        bootStrapStatsSpecificThresholds.prctile2p5_PPV_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_PPV_Balanced);
    statsForPPT(2).NPV = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_NPV_Balanced,...
        bootStrapStatsSpecificThresholds.prctile2p5_NPV_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_NPV_Balanced);
    %get threshold from nominal values
    statsForPPT(2).Threshold = sprintf('Balanced = %d',threshold_Balanced);

    %control the sig figs for threshold
    statsForPPT(1).Threshold = sprintf('Youden = %0.2f', abs(threshold_Youden));
    statsForPPT(2).Threshold = sprintf('Balanced = %0.2f', abs(threshold_Balanced));
    %add number pos and neg
    statsForPPT(1).numPos = T_labelCountsPerThresh.sum_TP_all(idxYouden) + T_labelCountsPerThresh.sum_FN_all(idxYouden);
    statsForPPT(1).numNeg = T_labelCountsPerThresh.sum_TN_all(idxYouden) + T_labelCountsPerThresh.sum_FP_all(idxYouden);
    statsForPPT(2).numPos = T_labelCountsPerThresh.sum_TP_all(idxBalanced) + T_labelCountsPerThresh.sum_FN_all(idxBalanced);
    statsForPPT(2).numNeg = T_labelCountsPerThresh.sum_TN_all(idxBalanced) + T_labelCountsPerThresh.sum_FP_all(idxBalanced);

end