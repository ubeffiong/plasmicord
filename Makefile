.PHONY: help demo test setup check run clean
help:
	@echo "make demo   - offline end-to-end on synthetic sequences (no tools)"
	@echo "make test   - unit tests (standalone + gap workflows + detailed report + clustering + network + discordance)"
	@echo "make setup  - create conda env"
	@echo "make check  - verify dependencies (stage 0)"
	@echo "make run    - full pipeline (stages 0-5)"
	@echo "make clean  - remove demo/tmp/pycache"
demo:  ; python3 plasmicord.py demo --out results_demo
test:
	python3 -m unittest discover -s test -p test_standalone.py -v
	python3 -m unittest discover -s test -p test_gap_workflows.py -v
	python3 -m unittest discover -s test -p test_detailed_report.py -v
	python3 test/test_clustering.py
	python3 test/test_network_discordance.py
setup: ; bash env/setup_conda.sh
check: ; bash scripts/run_all.sh 0
run:   ; bash scripts/run_all.sh
clean:
	rm -rf results_demo tmp
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
