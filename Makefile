OUTPUT = build
JOB = qce26_claim_against_measurement
.PHONY: all clean compile_plots dev repro repro_docker

# Programs
COMPOSE=docker/docker-compose.yml
DC=docker compose
PY=python
R=Rscript

# Directories
D_RESULTS=$(OUTPUT)/results
D_PLOTS=$(OUTPUT)/plots
D_PAPER=$(OUTPUT)/paper

D_REPRODUCTION=./reproduction
D_R=$(D_REPRODUCTION)/R
D_SCRIPTS=$(D_REPRODUCTION)/scripts
D_DATA=$(D_REPRODUCTION)/data


OUTDIRS=$(D_RESULTS) $(D_PLOTS) $(D_PAPER)

# Dependencies
RAW_PLOTS=review_compliance khan_sensitivity desdentado_sweep drift_combined_all drift_week drift_ibm_brussels
PLOTS=$(addprefix $(D_PLOTS)/,$(addsuffix .tex,$(RAW_PLOTS)))
DESDENTADO_RESULT=$(D_RESULTS)/desdentado_original_boxplot.csv
DESDENTADO_INPUT=$(wildcard $(D_DATA)/desdentado_original/summary_*shots.csv)
KHAN_RESULTS=$(D_RESULTS)/khan_detail.csv $(D_RESULTS)/khan_power.csv $(D_RESULTS)/khan_shots.csv $(D_RESULTS)/khan_summary.csv
MC_RESULT=$(D_RESULTS)/khan_multiple_comparisons.csv

# Drift statistics (Table tab:drift-summary: eta^2/ICC, r1, n_eff, Cohen's d range)
DRIFT_SESSIONS=first_run day2_full weekend
D_DRIFT_STATS=$(D_RESULTS)/drift_stats
DRIFT_STATS_RESULT=$(addsuffix /drift_stats_summary.csv,$(addprefix $(D_DRIFT_STATS)/,$(DRIFT_SESSIONS)))

# Rules
all: build/$(JOB).pdf

build/$(JOB).pdf: paper/main.tex paper/$(JOB).bbl | $(OUTDIRS)
	cp paper/$(JOB).bbl $(D_PAPER)/$(JOB).bbl
	BIBINPUTS=paper:$$BIBINPUTS latexmk -shell-escape -lualatex -e '$$bibtex_use=0' -output-directory=$(D_PAPER) -jobname=$(JOB) $<

dev: $(COMPOSE)
	$(DC) -f $^ run --rm $@

repro_docker: $(COMPOSE)
	# ensure existence of container
	$(DC) -f $^ build
	$(DC) -f $^ run --rm repro
	make

repro: plots $(MC_RESULT) $(DRIFT_STATS_RESULT)
	echo "Up to date!"

$(OUTDIRS):
	mkdir -p $@

plots: $(PLOTS) compile_plots

compile_plots: $(PLOTS) | $(D_PLOTS)
	$(D_REPRODUCTION)/gen_img.sh $(D_PLOTS)

$(D_PLOTS)/review_compliance.tex &: $(D_R)/plot_review_compliance.R $(D_DATA)/review_criteria.csv | $(OUTDIRS)
	$(R) $<

$(KHAN_RESULTS) &: $(D_SCRIPTS)/khan_backend.py | $(OUTDIRS)
	$(PY) $< --outdir $(OUTPUT)/results/

$(D_PLOTS)/khan_sensitivity.tex &: $(D_R)/plot_khan_backend.R $(KHAN_RESULTS) | $(OUTDIRS)
	$(R) $<
	sed -i 's/khan_heatmap_ras1/build\/plots\/khan_heatmap_ras1/' $(D_PLOTS)/khan_heatmap.tex
	sed -i 's/khan_heatmap_ras2/build\/plots\/khan_heatmap_ras2/' $(D_PLOTS)/khan_heatmap.tex

$(DESDENTADO_RESULT): $(D_SCRIPTS)/generate_desdentado_original_boxplot.py $(DESDENTADO_INPUT) | $(OUTDIRS)
	$(PY) $< --output $@

$(D_PLOTS)/desdentado_sweep.tex &: $(D_R)/plot_desdentado_original.R $(DESDENTADO_RESULT) | $(OUTDIRS)
	$(R) $<

$(D_PLOTS)/drift_combined_all.tex &: $(D_R)/plot_drift.R $(D_DATA)/qexa_drift/raw_data_first_run.csv $(D_DATA)/qexa_drift/raw_data_day2_full.csv $(D_DATA)/qexa_drift/raw_data_weekend.csv | $(OUTDIRS)
	$(R) $<

# Drift power/ICC/autocorrelation statistics, one run per session, logged separately
$(D_DRIFT_STATS)/%/drift_stats_summary.csv $(D_DRIFT_STATS)/%/drift_power_summary.txt &: $(D_SCRIPTS)/full_statistics.py $(D_SCRIPTS)/power_analysis.py $(D_DATA)/qexa_drift/raw_data_%.csv | $(OUTDIRS)
	mkdir -p $(D_DRIFT_STATS)/$*
	$(PY) $(D_SCRIPTS)/full_statistics.py --csv $(D_DATA)/qexa_drift/raw_data_$*.csv --outdir $(D_DRIFT_STATS)/$*
	$(PY) $(D_SCRIPTS)/power_analysis.py --csv $(D_DATA)/qexa_drift/raw_data_$*.csv --outdir $(D_DRIFT_STATS)/$*

$(D_PLOTS)/drift_week.tex &: $(D_R)/plot_drift_week.R $(D_DATA)/qexa_drift/raw_data_week.csv | $(OUTDIRS)
	$(R) $<

$(D_PLOTS)/drift_ibm_brussels.tex &: $(D_R)/plot_ibm_drift.R $(D_DATA)/ibm_drift_results/raw_data_ibm_drift.csv | $(OUTDIRS)
	$(R) $<

$(MC_RESULT): $(D_SCRIPTS)/multiple_comparisons.py $(D_RESULTS)/khan_summary.csv | $(OUTDIRS)
	$(PY) $<

clean:
	rm -rf build
