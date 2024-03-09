function Outputs = bootstrapBackwardAnalysisOutputs(labels,scores,patientIds,numBoostrapIterations)
%example how to use:
% scores = nx1 double array of parameterValues intended to predict the label
% labels = nx1 double array of 0 or 1
% patientIds = nx1 string array of patientIDs corresponding to each data point
% numBoostrapIterations = 2000; number of time to do random sampling of patients with replacement
% Outputs = doBootstrappingForBackwardAnalysisOutputs(labels,scores,patientIds,numBoostrapIterations)

labels = labels(:);
scores = scores(:);
patientIds = patientIds(:);


uniquePatientIds = unique(patientIds);
N = length(uniquePatientIds);

resultsAtSpecificThresholdsPerBootstrap = struct();
parfor icr=1:numBoostrapIterations%parfor
%     rng('default');%BS remove default 
    rng(icr,'twister');
    if mod(icr,100) == 0
        disp(icr);
    end

    %randomly sample patients with replacement
    scoresNew=[]; labelsNew=[];
    for icN = 1:N
        try
            LIA = ismember(patientIds, uniquePatientIds{ randperm(N,1) }); idxtmp = find(LIA==1);
            scoresNew  =[scoresNew;  scores(idxtmp)];
            labelsNew =[labelsNew; labels(idxtmp)];
        end
    end

    %compute classification metrics for this sample across all thresholds
    if max(labelsNew)==1 && min(labelsNew)==0
        [Spec_Hack, Sens_Hack, Thresh, AUC] = perfcurve(labelsNew, scoresNew, 0, 'xCrit', 'Spec','yCrit', 'Sens');
        %sens/spec correctly switched
        Spec = Sens_Hack;
        Sens = Spec_Hack;
        [PPV_Hack, NPV_Hack, Thresh2, ~] = perfcurve(labelsNew, scoresNew, 0, 'xCrit', 'PPV','yCrit', 'NPV');
        %sPPV/NPV correctly switched
        PPV = NPV_Hack;
        NPV = PPV_Hack;
        assert(isequaln(Thresh,Thresh2));
        if AUC<0.5
            [Spec, Sens, Thresh, AUC] = perfcurve(labelsNew, scoresNew, 1, 'xCrit', 'Spec','yCrit', 'Sens');
            [PPV, NPV, Thresh2, ~] = perfcurve(labelsNew, scoresNew, 1, 'xCrit', 'PPV','yCrit', 'NPV');
            assert(isequaln(Thresh,Thresh2));
        end
    else
        continue;
    end

    
    Recall = Sens;
    Precision = PPV;
    
    %store the results for this iteration
    resultsAtSpecificThresholdsPerBootstrap(icr).AUC = AUC;

    %---get the classification metric at optimal thresholds
    % for ROC
    %Youden Index
    tmp = Sens + Spec - 1;
    idxYouden = find(tmp==max(tmp)); idxYouden = idxYouden(1);
    resultsAtSpecificThresholdsPerBootstrap(icr).threshold_Youden = Thresh(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(icr).Sens_Youden = Sens(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(icr).Spec_Youden = Spec(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(icr).PPV_Youden = PPV(idxYouden);
    resultsAtSpecificThresholdsPerBootstrap(icr).NPV_Youden = NPV(idxYouden);
    %Balanced Sens,Spec
    tmp = abs(Spec - Sens);
    idxBalanced = find(tmp==min(tmp)); idxBalanced = idxBalanced(1);
    resultsAtSpecificThresholdsPerBootstrap(icr).threshold_Balanced = Thresh(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(icr).Sens_Balanced = Sens(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(icr).Spec_Balanced = Spec(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(icr).PPV_Balanced = PPV(idxBalanced);
    resultsAtSpecificThresholdsPerBootstrap(icr).NPV_Balanced = NPV(idxBalanced);
    
end
T_ResultsAtSpecificThresholdsAllBootstraps = struct2table(resultsAtSpecificThresholdsPerBootstrap);



%compute stats for AUC and specific thresholds
T_ResultsAtSpecificThresholdsAllBootstraps.dummyVar(:,1) = "sameGroup";
bootStrapStatsSpecificThresholds = grpstats(T_ResultsAtSpecificThresholdsAllBootstraps,"dummyVar",{@mean,@std,@median,@(x) quantile(x,0.025),@(x) quantile(x,0.975)});
varNamesTemp = bootStrapStatsSpecificThresholds.Properties.VariableNames;
bootStrapStatsSpecificThresholds = renamevars(bootStrapStatsSpecificThresholds,varNamesTemp(startsWith(varNamesTemp,'Fun4')),strrep(varNamesTemp(startsWith(varNamesTemp,'Fun4')),'Fun4','prctile2p5'));
bootStrapStatsSpecificThresholds = renamevars(bootStrapStatsSpecificThresholds,varNamesTemp(startsWith(varNamesTemp,'Fun5')),strrep(varNamesTemp(startsWith(varNamesTemp,'Fun5')),'Fun5','prctile97p5'));

Outputs = struct();
Outputs.T_ResultsAtSpecificThresholdsAllBootstraps = T_ResultsAtSpecificThresholdsAllBootstraps;
Outputs.bootStrapStatsSpecificThresholds = bootStrapStatsSpecificThresholds;
end