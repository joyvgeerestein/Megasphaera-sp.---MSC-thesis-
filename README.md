Megasphaera MJR8396C - MSc Thesis

This GitHub repository contains the code and data generated as part of the MSc thesis conducted by Joy van Geerestein.

**Initial Validation:**
The tool memote was used to test the viability of the provided genome-scale metabolic model (GEM) for Megasphaera MJR8396C, which yielded a satisfactory memote test score.

**Metabolic Flux Analysis:**
General metabolic fluxes in a gut-like medium were tested using pFBA (code: generalpfba).
The effects of varying amino acids and carbohydrates on the growth of Megasphaera MJR8396C were assessed (code: AAvscarbohydrates).
Attempts to add pectin to the GEM to test growth on this substrate were made, but did not yield feasible results (code: add_pectin_to_model).

**Stickland Metabolism Analysis:**
Stickland metabolism pairs were simulated by pairing two amino acids and conducting pFBA (code: SticklandpairspFBA).
Heatmaps were created using pandas and seaborn in Python to visualize interaction effects (code: interactioneffectheatmap).
Stickland reactions were searched for in all reactions in the GEM (code: reactionstoexcel) 
Specific Stickland metabolites without exchange reactions were added to the model, but this also did not yield feasible results (code: add_stickland_metabolites).

**Butyrate Production Analysis:**
Butyrate production was assessed using production envelope and double robustness analysis (code: productionenvelope and doublerobustnessanalysis).
Optimization algorithms were then employed (code: Optforcebutyrate and OptKnockbutyrate).

**Community Simulations:**
A multi-species model was created for community simulations with Bifidobacterium adolescentis (code: multiplespeciesmodel), resulting in the community model pairedModel_Bifidobacterium_adolescentis_BM.xml_Megasphaera_sp_MJR8396C.xml.mat.
FBA simulations on the community model were performed with initial minimal medium for each microbe (code: minimalmedium).
Pairwise interactions were analyzed (code: pairwiseinteraction), resulting in interactions pairwise.mat.
BacArena was used for community growth simulations over time (BacArenacode).

**Data Availability:**
All generated data from the above simulations can be found in the file "JVG MSc thesis data file"






