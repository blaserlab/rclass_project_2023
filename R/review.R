# Read in and inspect the configuration file.
analysis_configs <- read_csv(system.file("extdata/vignette_config.csv",
                                         package = "blaseRdata"),
                             col_type = list(.default = col_character()))

analysis_configs
readr::write_csv(analysis_configs, file = "~/network/X/Labs/Blaser/network_transfer/analysis_configs.csv")
readr::read_csv("~/network/X/Labs/Blaser/network_transfer/analysis_configs.csv", col_type = list(.default = col_character()))
# Fix the file path

cds_list <- map(.x = analysis_configs$sample,
                .f = \(x, conf = analysis_configs) {
                  conf_filtered <- conf |>
                    filter(sample == x)
                  targz <- list.files(
                    conf_filtered$targz_path,
                    pattern = "filtered_feature_bc_matrix.tar.gz",
                    recursive = T,
                    full.names = T
                  )
                  cds <- bb_load_tenx_targz(targz_file = targz,
                                            sample_metadata_tbl = conf_filtered |>
                                              select(-c(sample, targz_path)))
                  if ("Antibody Capture" %in% bb_rowmeta(cds)$data_type) {
                    cds <- bb_split_citeseq(cds)
                  }
                  return(cds)
                }) %>%
  set_names(nm = analysis_configs$sample)
