%% productionenvelope
lineColor ='b';
targetRxn = 'EX_cpd00211_e0';
biomassRxn = 'EX_cpd11416_c0';
geneDelFlag = false;
deletions = {};
nPts = 50;
[biomassValues,targetValues,lineHandle]=productionEnvelope(model,deletions,lineColor,targetRxn,biomassRxn,geneDelFlag,nPts);
xlabel('Biomass (mmol/gdw*h)')
ylabel('NH3 (mmol/gDW*h)')