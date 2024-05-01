# Script to reproduce Megasphaera sp. and Bifidobacterium simulation in essential  medium
library(BacArena)
library(sybilSBML)
library(data.table)

Mg <- readSBMLmod("/Users/joyvangeerestein/Documents/MSC food tech /MSC thesis /visual code/Megasphaera_sp_MJR8396C.xml")
Bi <- readSBMLmod("/Users/joyvangeerestein/Documents/MSC food tech /MSC thesis /visual code/Bifidobacterium_adolescentis.xml")

Mg@mod_desc <- "Megasphaera"
Bi@mod_desc <- "Bifido"

mega <- Bac(Mg)
bifi <- Bac(Bi)

arena <- Arena(n = 100, m = 100)
arena <- addOrg(arena, mega, amount = 2)
arena <- addOrg(arena, bifi, amount = 2)

#add medium high AA 
arena_subs <- fread("/Users/joyvangeerestein/Downloads/highAANoacetate.csv") 
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]
arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "fmol/cell", addAnyway = T)

# OR for default gut medium
# arena <- addDefaultMed(arena,bac,unit="fmol/cell")

CF_sim <- simEnv(arena, time = 24, sec_obj = "mtf")

# plot growth curves and levels of acetate, butryate, propionate, lactate, succinate 
par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd00211_e0", "cpd00029_e0", "cpd00141_e0", "cpd00159_e0", "cpd00036_e0", "cpd00049_e0"), 
            dict = list(cpd00211_e0 = "Butryate",
                        cpd00029_e0 = "Acetate",
                        cpd00141_e0 = "Propionate",
                        cpd00159_e0 = "Lactate",
                        cpd00036_e0 = "Succinate",))