cd tests
COVERAGE_FILE=".coverage"
python -m pytest --cov --cov-config=../.coveragerc -vv
python -m coverage combine
python -m coverage json
python -m coverage report --fail-under=100
unset COVERAGE_FILE
cd ..
