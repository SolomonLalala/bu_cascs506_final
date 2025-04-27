# Makefile for bu_cascs506_final project

# Define file names and paths
RMD = final_report.Rmd
HTML = final_report.html
CV_SCRIPT = scripts/cv.R
BOOTSTRAP_SCRIPT = scripts/bootstrap.R

# Default target
all: run_cv run_bootstrap build_report

# Run cross-validation
run_cv:
	Rscript $(CV_SCRIPT)

# Run bootstrap feature selection
run_bootstrap:
	Rscript $(BOOTSTRAP_SCRIPT)

# Knit the RMarkdown report
build_report:
	Rscript -e "rmarkdown::render('$(RMD)', output_format = 'html_document')"

# Clean up generated files
clean:
	rm -f data/*.rds
	rm -f $(HTML)

# Save session information to a file (optional)
session_info:
	Rscript -e "writeLines(capture.output(sessionInfo()), 'session_info.txt')"

.PHONY: all run_cv run_bootstrap build_report clean session_info
