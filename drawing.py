import matplotlib.pyplot as plt
import numpy as np

from seqhandling import pattern_count_kprototype, pattern_count_prototype


def lines(data, xlabel="X", ylabel="Y", main="", markers=None, show=True, new=True):
    """
    Generic method to draw an array of values
    """
    if markers is None:
        markers = []
    x = [float(i) for i in range(0, len(data))]
    y = [float(data[i]) for i in range(0, len(data))]

    if new:
        plt.close('all')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(main)
    plt.plot(x, y, '-o', markevery=markers)
    plt.xticks(np.arange(0, len(data) + 1, 5.0))
    plt.grid(axis="x")
    if show:
        plt.show()


def plot_clustering(clustering, Kprototypes, f0, f1, elbow):
    # data
    x = range(len(clustering))
    # norm_rates = [Kprototypes[i].norm_rate for i in range(len(Kprototypes))]
    rates = [Kprototypes[i].mean_rate for i in range(len(Kprototypes))]

    #  Plot
    plt.close('all')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))

    # Fig #0
    ax[0].plot(x, [clustering[i].mean for i in range(len(clustering))])
    ax[0].plot(x, [clustering[i].median for i in range(len(clustering))])
    ax[0].plot(x, [clustering[i].stdev for i in range(len(clustering))])
    ax[0].axvline(x=elbow, color='0.75', linestyle='--')
    ax[0].axvline(x=f0, color='0.75', linestyle='--')
    ax[0].axvline(x=f1, color='0.75', linestyle='--')
    ax[0].legend(['Mean', 'Median', 'Stdev'], loc='upper right')
    ax[0].set_xticks(np.arange(0, len(x) + 1, 5.0))
    ax[0].set_title("Intra distance metrics at each clustering.")

    # Fig #1
    # ax[1].plot(x, norm_rates, '.', color="blue", alpha=0.75)
    ax[1].plot(x, rates, '.', color="lightblue", alpha=0.75)
    ax[1].axvline(x=elbow, color='0.75', linestyle='--')
    ax[1].axvline(x=f0, color='0.75', linestyle='--')
    ax[1].axvline(x=f1, color='0.75', linestyle='--')
    ax[1].set_xticks(np.arange(0, len(x) + 1, 5.0))
    ax[1].set_title("Rate weighted by nr of clusters SCORE")
    fig.tight_layout()
    plt.show()


def hbar_prototype_count(data, prototypes, threshold=0):
    plt.close('all')
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))

    if threshold == 0:
        threshold = len(data) * 0.2

    pc = pattern_count_prototype(data, prototypes)
    pc = list(filter(lambda x: x[1] > threshold, pc))
    ax.barh(range(0, len(pc)), [pc[i][1] for i in range(len(pc))], color="lightblue")
    ax.set_yticks(range(len(pc)))
    ax.set_yticklabels([pc[i][0] for i in range(len(pc))], fontsize='small', family='monospace')
    ax.set_title("Prototypes of Best Clustering")

    plt.show()


def hbar_prototypes_count(data, Kprototype, f0, f1, elbow, threshold=0):
    plt.close('all')

    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    if threshold == 0:
        threshold = len(data) * 0.2

    pc = pattern_count_kprototype(data, Kprototype[elbow])
    pc = list(filter(lambda x: x[1] > threshold, pc))
    ax[0].barh(range(0, len(pc)), [pc[i][1] for i in range(len(pc))], color="lightblue")
    ax[0].set_yticks(range(len(pc)))
    ax[0].set_yticklabels([pc[i][0] for i in range(len(pc))], fontsize='small', family='monospace')
    ax[0].set_title("Elbow")

    pc = pattern_count_kprototype(data, Kprototype[f0])
    pc = list(filter(lambda x: x[1] > threshold, pc))
    ax[1].barh(range(0, len(pc)), [pc[i][1] for i in range(len(pc))], color="lightblue")
    ax[1].set_yticks(range(len(pc)))
    ax[1].set_yticklabels([pc[i][0] for i in range(len(pc))], fontsize='small')
    ax[1].set_title("F0")

    pc = pattern_count_kprototype(data, Kprototype[f1])
    pc = list(filter(lambda x: x[1] > threshold, pc))
    ax[2].barh(range(0, len(pc)), [pc[i][1] for i in range(len(pc))], color="lightblue")
    ax[2].set_yticks(range(len(pc)))
    ax[2].set_yticklabels([pc[i][0] for i in range(len(pc))], fontsize='small')
    ax[2].set_title("Fscore")

    plt.tight_layout()
    plt.show()


def histogram_fscores(elbs, f0s, f1s, nbins=25):
    if len(elbs) < 10: return
    plt.close('all')
    if len(elbs) < 25:
        nbins = len(elbs)
    else:
        nbins = int(len(elbs))
    hist, bins = np.histogram(elbs, bins=nbins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width, color="green", alpha=0.7)

    hist, bins = np.histogram(f0s, bins=nbins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width, color="blue", alpha=0.7)

    hist, bins = np.histogram(f1s, bins=nbins)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width, color="red", alpha=.7)
    plt.legend(["Elbow", "F0", "F1"])
    plt.show()


def plot_score_rate(k_prototypes):
    """
    Rate is the ratio between the number of sequences matched by a given prototype extracted from a cluster,
    and the total number of sequences in that cluster. A low rate is a bad prototyping or clustering result.
    The Score is the number of matches of each prototype within the cluster.
    The graph shows the average values for Mean and Rate for each clustering. Each clustering tries to group
    with a different value of 'k' (X-Axis).
    """
    x_width = len(k_prototypes)
    lines([k_prototypes[i].score for i in range(x_width)], show=False)
    lines([k_prototypes[i].rate for i in range(x_width)],
          ylabel="Avg. score & rate", xlabel="# clusters", show=False, new=False)
    plt.legend(['Score', 'Rate'], loc='upper right')
    plt.axhline(1.0, color='0.75', linestyle='--')
    plt.xticks(np.arange(0, x_width + 1, 5.0))
    plt.show()
