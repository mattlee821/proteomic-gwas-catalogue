# Quarto pre-render hook: validates catalogue/*.tsv and compiles site/data/catalogue.json
# Quarto executes pre-render scripts from the project directory (site/).
old_wd <- setwd("..")
on.exit(setwd(old_wd), add = TRUE)
status <- system2("Rscript", "src/build_site_data.R")
if (status != 0L) stop("src/build_site_data.R failed during Quarto pre-render", call. = FALSE)
