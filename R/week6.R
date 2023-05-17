usethis::edit_r_profile()

blaseRtemplates::initialize_package(path = "/workspace/brad_workspace/workshop.data")

dir.create("data")
# rename the object for saving
workshop_cds <- vignette_cds

save(workshop_cds, file = "data/workshop_cds.rda", compress = "bzip2")
save(analysis_configs, file = "data/analysis_configs.rda", compress = "bzip2")


analysis_configs
workshop_cds
