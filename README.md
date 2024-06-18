Objective: To perform on diagnostic on the algoriitm: detect those particles that wormsim fails to run

Description:

The file contents are exactly the same as in task 2.2. Differences are:

smc_app_dependent_functions_v2_Peri_one_particle_at_a_time.r: smc dependent files linked to distance metric used the difference in the proportion of prevalence between simulations and observations, summed squared, this script now removes any parallelization and tries in sequence on a single core every particle proposed.

It generates two reports (located within 04_wormsim-v2.58Ap27), Report.txt which contains all the particles tested and the failed.runs.log.txt which only includes that no simulated intensity table is generated.

The contents of the failed.runs.log.txt include: i" "seed" "mbr" "exposure.p1
