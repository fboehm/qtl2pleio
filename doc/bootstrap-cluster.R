## -----------------------------------------------------------------------------
fn <- file.path("https://raw.githubusercontent.com", 
                "fboehm/qtl2pleio-manuscript/master/chtc", 
                "Recla-bootstrap/Rscript/boot-Recla-10-22.R"
                )
foo <- readLines(fn)

## ---- echo = FALSE, comment = ""----------------------------------------------
cat(foo, sep = "\n")

## -----------------------------------------------------------------------------
devtools::session_info()

