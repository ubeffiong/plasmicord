.PHONY: help demo full-demo test real-cohort-test setup check run clean
help:
	@echo "make demo             - offline end-to-end on synthetic sequences (no tools)"
	@echo "make full-demo        - offline synthetic demo exercising every dashboard feature"
	@echo "make test             - unit tests (standalone + gap workflows + detailed report + CLI end-to-end + edge cases + documentation links + clustering + network + discordance)"
	@echo "make real-cohort-test - opt-in: download a small real cohort from NCBI and run the pipeline (network required, not part of make test/CI)"
	@echo "make setup            - create conda env"
	@echo "make check            - verify dependencies (stage 0)"
	@echo "make run              - full pipeline (stages 0-5)"
	@echo "make clean            - remove demo/tmp/pycache"
demo:  ; python3 plasmicord.py demo --out results_demo
full-demo: ; python3 scripts/full_feature_demo.py --out results_full_demo
test:
	python3 -m unittest discover -s test -p test_standalone.py -v
	python3 -m unittest discover -s test -p test_gap_workflows.py -v
	python3 -m unittest discover -s test -p test_detailed_report.py -v
	python3 -m unittest discover -s test -p test_cli_end_to_end.py -v
	python3 -m unittest discover -s test -p test_edge_cases.py -v
	python3 -m unittest discover -s test -p test_documentation_links.py -v
	python3 test/test_clustering.py
	python3 test/test_network_discordance.py
real-cohort-test:
	PLASMICORD_REAL_COHORT_TEST=1 python3 -m unittest discover -s test -p test_real_cohort.py -v
setup: ; bash env/setup_conda.sh
check: ; bash scripts/run_all.sh 0
run:   ; bash scripts/run_all.sh
clean:
	rm -rf results_demo results_full_demo tmp
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
