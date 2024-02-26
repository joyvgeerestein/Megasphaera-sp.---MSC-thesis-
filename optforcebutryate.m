% initialize the toolbox
global TUTORIAL_INIT_CB;
if ~isempty(TUTORIAL_INIT_CB) && TUTORIAL_INIT_CB==1
    initCobraToolbox(false) % false, as we don't want to update
end
 
changeCobraSolver('gurobi','all');
%% 
% load model
model = readCbModel('Megasphaera_sp_MJR8396C.xml');
biomass = 'EX_cpd11416_c0';


growthRate = optimizeCbModel(model); 
fprintf('The maximum growth rate is %1.2f', growthRate.f);

%calculate maximum production for butryate
model = changeObjective(model, 'EX_cpd00211_e0');
maxButr = optimizeCbModel(model);
fprintf('The maximum production rate of butryate is %1.2f', maxButr.f);

%% 
% The WT strain's biomass function ("EX_cpd11416_c0") is constrained to near the maximum growth rate.
constrWT = struct('rxnList', {{'EX_cpd11416_c0'}}, 'rxnValues', 1.11, 'rxnBoundType', 'b')
% The mutant strain's biomass function is set to 10% of the wildtype.
% Butryate is forced to be minimum of 19 (max flux from optknock), max
% production of NH3 and H2s
constrMT = struct('rxnList', {{'EX_cpd11416_c0', 'EX_cpd00211_e0', 'EX_cpd00013_e0', 'EX_cpd00239_e0'}}, 'rxnValues', [0.11, 19, 0.01, 0.01], ...
 'rxnBoundType', 'bbuu')

%% Flux Variability Analysis
[minFluxesW, maxFluxesW, minFluxesM, maxFluxesM] = FVAOptForce(model,constrWT,constrMT)
%% 
disp([minFluxesW, maxFluxesW, minFluxesM, maxFluxesM]);
%% define runID

runID = 'TestOptForceM';

%% constraints

constrOpt = struct('rxnList', {{'EX_cpd11416_c0', 'EX_cpd00211_e0', 'EX_cpd00013_e0', 'EX_cpd00239_e0'}}, 'values', [0.11, 19, 0.01, 0.01]);
%%  finding first order must sets

%MustL set decrease
[mustLSet, pos_mustL] = findMustL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
 'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
 'outputFileName', 'MustL' , 'printExcel', 1, 'printText', 1, ...
 'printReport', 1, 'keepInputs', 1, 'verbose', 0);
 %% 
disp(mustLSet)
%% 

% MustU set increase 
[mustUSet, pos_mustU] = findMustU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
 'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
 'outputFileName', 'MustU' , 'printExcel', 1, 'printText', 1, ...
 'printReport', 1, 'keepInputs', 1, 'verbose', 0);
%% 
disp(mustUSet)

%% finding second order must sets
constrOpt = struct('rxnList', {{'EX_cpd11416_c0', 'EX_cpd00211_e0', 'EX_cpd00013_e0', 'EX_cpd00239_e0'}}, 'values', [1.1, 11.4, 0.01, 0.01]);

exchangeRxns = model.rxns(cellfun(@isempty, strfind(model.rxns, 'EX_')) == 0);
excludedRxns = unique([mustUSet; mustLSet; exchangeRxns]);
%% 

% MustUU:
[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
    findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
 'excludedRxns', excludedRxns,'runID', runID, ...
 'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
 'verbose', 1);
%% 

%MustLL
[mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
    findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
 'excludedRxns', excludedRxns,'runID', runID, ...
 'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
 'verbose', 1);
%% 

%MustUL
[mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = ...
    findMustUL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
 'excludedRxns', excludedRxns,'runID', runID, ...
 'outputFolder', 'OutputsFindMustUL', 'outputFileName', 'MustUL', ...
 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
 'verbose', 1);

%% optforce

mustU = unique(union(mustUSet, mustUU));
mustL = unique(union(mustLSet, mustLL));
targetRxn = 'EX_cpd00211_e0';
biomassRxn = 'EX_cpd11416_c0';
k = 1;
nSets = 1;
constrOpt = struct('rxnList', {{'EX_cpd11416_c0'}}, 'values', [0.11]');
 
[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
 'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
 'runID', runID, 'outputFolder', 'OutputsOptForce', ...
 'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
 'printReport', 1, 'keepInputs', 1, 'verbose', 1);

%% display reactions found 
disp(optForceSets)

%transporter for butryate, so exclude this one from the next
%solution 
%% 
k = 2;
nSets = 10;
runID = 'TestOptForceM2';
excludedRxns = struct('rxnList', {{'rxn05683_c0'}}, 'typeReg','U');
[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
 'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
 'excludedRxns', excludedRxns, ...
 'runID', runID, 'outputFolder', 'OutputsOptForce', ...
 'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
 'printReport', 1, 'keepInputs', 1, 'verbose', 1);
%% 
disp(optForceSets)

%% Sets found by optForce were printed in OptForce.xls  
Sets found by optForce were printed in OptForce.txt  
    {'rxn05683_c0'}    {'rxn14246_c0'   }
    {'rxn05683_c0'}    {'rxn01739_c0'   }
    {'rxn05683_c0'}    {'rxn00870_c0'   }
    {'rxn05683_c0'}    {'rxn00783_c0'   }
    {'rxn00178_c0'}    {'rxn03245_c0'   }
    {'rxn05683_c0'}    {'rxn05242_c0'   }
    {'rxn05683_c0'}    {'EX_cpd00001_e0'}
    {'rxn05683_c0'}    {'rxn00875_c0'   }
    {'rxn05683_c0'}    {'rxn03861_c0'   }
    {'rxn05683_c0'}    {'rxn15962_c0'   }




