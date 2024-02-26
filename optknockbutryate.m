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

% max solutions to find; (i.e., maximum number of remvable reactions that lead to the overproduction of the metabolite of interest)
treshold = 10; 

selectedRxnList = model.rxns

% minimize NH3 and H2S by setting upper bounds to small value 
model = changeRxnBounds(model, 'EX_cpd00013_e0', 0.01, 'u');
model = changeRxnBounds(model, 'EX_cpd00239_e0', 0.01, 'u');

% determine butryate production and growth rate
fbaWT = optimizeCbModel(model);
butrFluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00211_e0'));
nh3FluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00013_e0'));
h2sFluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00239_e0'));
succFluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00036_e0'));
acetFluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00029_e0'));
propiFluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00141_e0'));
formFluxWT = fbaWT.x(strcmp(model.rxns, 'EX_cpd00047_e0'));

growthRateWT = fbaWT.f;
fprintf('The production of butryate before optimization is %.1f \n', butrFluxWT);
fprintf('The growth rate before optimization is %.1f \n', growthRateWT);
fprintf(['The production of other products such as NH3, H2S, succinate, acetate, propionate, formate are %.1f, %.1f, %.1f and %.1f, respectively. \n'], ...
        nh3FluxWT, h2sFluxWT, succFluxWT, acetFluxWT, propiFluxWT, formFluxWT);
%% 
 

% finding optknock solutions of size 2 for increasing production of
% butryate
fprintf('\n...EXAMPLE 1: Finding optKnock sets of size 2 or less...\n\n')


% Set optKnock options
% The exchange of butryate will be the objective of the outer problem
options = struct('targetRxn', 'EX_cpd00211_e0', 'numDel', 2);
treshold = 10;
% We will impose that biomass be at least 10% of the biomass of wild-type
constrOpt = struct('rxnList', {{biomass}},'values', 0.1*fbaWT.f, 'sense', 'G');
% We will try to find 10 optKnock sets of a maximun length of 2
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 1;
while nIter < treshold
    fprintf('...Performing optKnock analysis...\n')
 if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt);
 else
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt, previousSolutions, 1);
 end
 
 % determine butryate production and growth rate after optimization
    butrFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00211_e0'));
    growthRateM1 = optKnockSol.fluxes(strcmp(model.rxns, biomass));
    nh3FluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00013_e0'));
h2sFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00239_e0'));
succFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00036_e0'));
acetFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00029_e0'));
propiFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00141_e0'));
formFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_cpd00047_e0'));
    setM1 = optKnockSol.rxnList;
 
 if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions = contPreviousSolutions + 1;
 %printing results
        fprintf('optKnock found a optKnock set of large %d composed by ', length(setM1));
 for j = 1:length(setM1)
 if j == 1
                fprintf('%s', setM1{j});
 elseif j == length(setM1)
                fprintf(' and %s', setM1{j});
 else
                fprintf(', %s', setM1{j});
 end
 end
        fprintf('\n');
        fprintf('The production of butryate after optimization is %.2f \n', butrFluxM1);
        fprintf('The growth rate after optimization is %.2f \n', growthRateM1);
        fprintf(['The production of other products such as NH3, H2S, succinate, acetate, propionate, formate are' ...
 '%.1f, %.1f, %.1f and %.1f, respectively. \n'], nh3FluxM1, h2sFluxM1, succFluxM1, acetFluxM1, propiFluxM1, formFluxM1);
        fprintf('...Performing coupling analysis...\n');
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, setM1, 'EX_cpd00211_e0');
        fprintf('The solution is of type: %s\n', type);
        fprintf('The maximun growth rate given the optKnock set is %.2f\n', maxGrowth);
        fprintf(['The maximun and minimun production of butryate given the optKnock set is ' ...
 '%.2f and %.2f, respectively \n\n'], minProd, maxProd);
 if strcmp(type, 'growth coupled')
            singleProductionEnvelope(model, setM1, 'EX_cpd00211_e0', biomass, 'savePlot', 1, 'showPlot', 1, ...
 'fileName', ['butryate_ex1_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
 end,
 else
 if nIter == 1
            fprintf('optKnock was not able to found an optKnock set\n');
 else
            fprintf('optKnock was not able to found additional optKnock sets\n');
 end
 break;
 end
    nIter = nIter + 1;
end