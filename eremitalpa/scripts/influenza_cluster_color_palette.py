import eremitalpa as ere
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main():
    _, ax = plt.subplots(figsize=(10, 1))
    for i, cluster in enumerate(ere.influenza.clusters):
        ax.text(i + 0.5, 1.05, cluster, va="bottom", fontsize=8, ha="center")

        try:
            fc = cluster.color
        except KeyError:
            ax.text(i + 0.5, 0.5, "Unknown", fontsize=4, ha="center")
            continue

        ax.add_artist(Rectangle((i, 0), width=1, height=1, facecolor=fc))
        ax.text(i + 0.5, -0.05, cluster.color, va="top", fontsize=4, ha="center")

    ax.set_xlim(0, i + 1)
    ax.set_ylim(0, 1.5)
    ax.axis("off")
    plt.savefig("cluster-colors.pdf", bbox_inches="tight")
