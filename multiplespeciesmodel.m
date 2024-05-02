%initialize toolbox
initCobraToolbox

% initialize the toolbox
global TUTORIAL_INIT_CB;
if ~isempty(TUTORIAL_INIT_CB) && TUTORIAL_INIT_CB==1
    initCobraToolbox(false) % false, as we don't want to update
end
 
changeCobraSolver('gurobi','all');
%% 
bifido = readCbModel('Bifidobacterium_adolescentis_ATCC_15703.xml');
megasphaera = readCbModel('Megasphaera_sp_MJR8396C.xml');

models{1,1} = bifido;
models{2,1} = megasphaera;
%% 
biomasses{1,1} = 'EX_cpd11416_c0';
biomasses{2,1} = 'EX_cpd11416_c0';

%% create joint model

[modelJoint] = createMultipleSpeciesModel(models);
%% set objective on maximizing both biomasses

biomassesJoint{1,1} = 'model1_DM_cpd11416[c0]';
biomassesJoint{2,1} = 'model2_DM_cpd11416[c0]';

modelJoint = changeObjective(modelJoint,biomassesJoint);
%% FBA
fbasolutionmodelJoint = optimizeCbModel(modelJoint)