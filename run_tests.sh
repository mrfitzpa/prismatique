cd tests
pytest --cov --cov-config=../.coveragerc -vv
coverage combine
coverage json
coverage report --fail-under=100
cd ..
