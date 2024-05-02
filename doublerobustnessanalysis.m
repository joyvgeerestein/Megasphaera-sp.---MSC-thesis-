%% double robustness analysis
model = changeRxnBounds(model,'EX_cpd00051_e0',-100,'l');
model = changeRxnBounds(model,'EX_cpd00132_e0',-10,'l');
controlRxn1 = 'EX_cpd00051_e0'; %arginine
controlRxn2 = 'EX_cpd00132_e0'; %asparagine
nPoints = 100;
objRxn = 'EX_cpd00211_e0'; %butyrate
plotResFlag = true;
objType = 'max';

[controlFlux1,controlFlux2,objFlux]=doubleRobustnessAnalysis(model,controlRxn1,controlRxn2,nPoints,plotResFlag,objRxn,objType);
xlabel('Arginine uptake (mmol/gDW*h)')
ylabel('Asparagine uptake (mmol/gDW*h)')
zlabel('Butyrate production (mmol/gDW*h)')