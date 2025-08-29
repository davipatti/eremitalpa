# eremitalpa

[![tests](https://github.com/davipatti/eremitalpa/actions/workflows/run-tests.yml/badge.svg)](https://github.com/davipatti/eremitalpa/actions/workflows/run-tests.yml)

Plot dendropy trees.

## Installation

```bash
pip install eremitalpa
```

## Basic usage:

```python
import eremitalpa as ere
tree = ere.Tree.from_disk("tree.tre")
ere.plot_dendropy_tree(tree)
```

Check out the [docs](https://davipatti.github.io/eremitalpa/).