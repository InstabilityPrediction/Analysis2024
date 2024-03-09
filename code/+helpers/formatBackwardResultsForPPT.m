function statsForPPT = formatBackwardResultsForPPT(bootStrapStatsSpecificThresholds,labels,alternateRoundingMethod)
%bootStrapStatsSpecificThresholds = table of bootstrap stats

%example how to use:



%format results for PPT
%columns: AUC, Sens,Spec,PPV,NPV, Threshold
%rows: Youden, Balanced
%each cell has the boostrap results aka median [2.5%tile, 97.5%tile], except the threshold rows just have the nominal values for the original patients


if ~exist('alternateRoundingMethod','var') || isempty(alternateRoundingMethod)
    alternateRoundingMethod = 1;%this was used to print to screen round( (median(newSpec))*1000)/1000;
end 

    statsForPPT = struct('AUC',[],'Sens',[],'Spec',[],'PPV',[],'NPV',[],'Threshold',[]);

    threshold_Youden = bootStrapStatsSpecificThresholds.median_threshold_Youden;
    threshold_Balanced = bootStrapStatsSpecificThresholds.median_threshold_Balanced;
    %Youden first
    if alternateRoundingMethod
        statsForPPT(1).AUC = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_AUC*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_AUC*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_AUC*1000)/1000);
        statsForPPT(1).Sens = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_Sens_Youden*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_Sens_Youden*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_Sens_Youden*1000)/1000);
        statsForPPT(1).Spec = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_Spec_Youden*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_Spec_Youden*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_Spec_Youden*1000)/1000);
        statsForPPT(1).PPV = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_PPV_Youden*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_PPV_Youden*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_PPV_Youden*1000)/1000);
        statsForPPT(1).NPV = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_NPV_Youden*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_NPV_Youden*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_NPV_Youden*1000)/1000);
    else
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
    end
    %get threshold from bootstrapping too
    statsForPPT(1).Threshold = sprintf('Youden = %d',threshold_Youden);
    
    %Balanced second, AUC is the same
    statsForPPT(2).AUC = statsForPPT(1).AUC;
    if alternateRoundingMethod
        statsForPPT(2).AUC = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_AUC*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_AUC*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_AUC*1000)/1000);
        statsForPPT(2).Sens = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_Sens_Balanced*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_Sens_Balanced*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_Sens_Balanced*1000)/1000);
        statsForPPT(2).Spec = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_Spec_Balanced*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_Spec_Balanced*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_Spec_Balanced*1000)/1000);
        statsForPPT(2).PPV = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_PPV_Balanced*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_PPV_Balanced*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_PPV_Balanced*1000)/1000);
        statsForPPT(2).NPV = sprintf('%0.3f [%0.3f, %0.3f]',round(bootStrapStatsSpecificThresholds.median_NPV_Balanced*1000)/1000,...
            round(bootStrapStatsSpecificThresholds.prctile2p5_NPV_Balanced*1000)/1000,round(bootStrapStatsSpecificThresholds.prctile97p5_NPV_Balanced*1000)/1000);
    else
        statsForPPT(2).Sens = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_Sens_Balanced,...
            bootStrapStatsSpecificThresholds.prctile2p5_Sens_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_Sens_Balanced);
        statsForPPT(2).Spec = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_Spec_Balanced,...
            bootStrapStatsSpecificThresholds.prctile2p5_Spec_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_Spec_Balanced);
        statsForPPT(2).PPV = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_PPV_Balanced,...
            bootStrapStatsSpecificThresholds.prctile2p5_PPV_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_PPV_Balanced);
        statsForPPT(2).NPV = sprintf('%0.3f [%0.3f, %0.3f]',bootStrapStatsSpecificThresholds.median_NPV_Balanced,...
            bootStrapStatsSpecificThresholds.prctile2p5_NPV_Balanced,bootStrapStatsSpecificThresholds.prctile97p5_NPV_Balanced);
    end
    %get threshold from bootstrapping too
    statsForPPT(2).Threshold = sprintf('Balanced = %d',threshold_Balanced);
       
    %control the sig figs for threshold
    statsForPPT(1).Threshold = sprintf('Youden = %0.2f', abs(threshold_Youden));
    statsForPPT(2).Threshold = sprintf('Balanced = %0.2f', abs(threshold_Balanced));
    %add number pos and neg
    statsForPPT(1).numPos = numel(find(labels == 1));
    statsForPPT(1).numNeg = numel(find(labels == 0));
    statsForPPT(2).numPos = numel(find(labels == 1));
    statsForPPT(2).numNeg = numel(find(labels == 0));

end
