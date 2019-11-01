# eremitalpa

Plot dendropy trees.

#### Installation

Dependencies
- [matplotlib](https://matplotlib.org/)
- [dendropy](https://dendropy.org/)

```bash
pip3 install git+https://github.com/davipatti/eremitalpa.git
```

#### Basic usage:

```python
import dendropy
import eremitalpa
tree = dendropy.Tree.get("tree.tre", schema="newick")
eremitalpa.plot_dendropy_tree(tree)
```
