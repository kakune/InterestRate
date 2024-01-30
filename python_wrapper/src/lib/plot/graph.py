import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plotGraph(
    xs: np.array,
    ys: np.array,
    fig: matplotlib.figure.Figure = None,
    ax: matplotlib.axes._axes.Axes = None,
    title: str = None,
    title_fontsize: int = 12,
    xlabel: str = None,
    ylabel: str = None,
    label_fontsize: int = 12,
    x_reg: tuple = None,
    y_reg: tuple = None,
    label: str = None,
    legend_fontsize: int = 6,
    color: str = None,
    linestyle: str = 'solid',
    marker: str = None,
    legend_loc: str = None,
    figsize: tuple[float, float] = (3.375, 2.08586),
    dpi: int = 600
):
    plt.rcParams['font.family'] = 'Arial'
    if (fig is None or ax is None):
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    ax.plot(xs, ys, linestyle=linestyle,
            marker=marker, label=label, color=color)

    if (title is not None):
        ax.set_title(title, fontsize=title_fontsize)
    if (xlabel is not None):
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if (ylabel is not None):
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if (x_reg is not None):
        ax.set_xlim(x_reg[0], x_reg[1])
    if (y_reg is not None):
        ax.set_ylim(y_reg[0], y_reg[1])

    if (label is not None):
        ax.legend(loc=legend_loc, fontsize=legend_fontsize)
        plt.legend(fontsize=legend_fontsize, bbox_to_anchor=(
            1, 1), loc='upper right', borderaxespad=0)

    return fig, ax
