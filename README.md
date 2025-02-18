# eremitalpa

Plot dendropy trees.

## Installation

Dependencies

- [matplotlib](https://matplotlib.org/)
- [dendropy](https://dendropy.org/)

```bash
pip install git+https://github.com/davipatti/eremitalpa.git
```

## Basic usage:

```python
import eremitalpa as ere
tree = ere.Tree.from_disk("tree.tre")
ere.plot_dendropy_tree(tree)
```


Check out the [docs](https://davipatti.github.io/eremitalpa/).