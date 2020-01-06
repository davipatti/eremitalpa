.PHONY: test check

check:
	pyflakes **/*.py
	pycodestyle **/*.py

test:
	pytest
