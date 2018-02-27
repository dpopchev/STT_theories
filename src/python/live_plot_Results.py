#!/usr/bin/env python

# Based on
#   https://pythonprogramming.net/live-graphs-matplotlib-tutorial/

file_path = "/home/dimitar/projects/STT_theories/results/"
file_name_results = "STT_phiScal_"

system_names = [ "r", "phiScal", "Q", "p", "LambdaMetr", "m"]

def my_ploting(ax, label_x, x, label_y, y):

    label_fontsize = 12
    ticks_label_size = 10

    ax.clear()

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=label_fontsize)
    ax.xaxis.set_tick_params(labelsize=ticks_label_size)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=label_fontsize)
    ax.yaxis.set_tick_params(labelsize=ticks_label_size)

    try:
        ax.plot(
          x, y,
          marker="o", markersize=7,
          linewidth=1.5
        )
    except ValueError:
        x = []
        y = []
        pass

def results_live(arg):

    import glob
    import os

    all_files = glob.glob(file_path + file_name_results + "*")
    while not len(all_files):
        all_files = glob.glob(file_path + file_name_results + "*")

    file_to_use = max(all_files, key=os.path.getctime)

    with open(file_to_use,'r') as f:
        graph_data = f.read()

    graph_data = [ k.split("\n") for k in graph_data.split("# ") if len(k) > 1 ][1]

    plot_labels = graph_data.pop(0).split(" ")

    plots_all = [ [] for k in range(len(plot_labels)) ]

    try:
        for single_line in graph_data:
            if len(single_line) > 1:
                tmp = [ j for j in single_line.split(" ") if len(j) > 1 ]
                for m,n in zip(plots_all, tmp):
                    m.append(float(n))
    except ValueError:
        graph_data = []
        plot_labels = []
        plots_all = []
        pass

    my_ploting(
      ax_phiScal,
      plot_labels[0], plots_all[0],
      plot_labels[1], plots_all[1],
    )

    my_ploting(
      ax_mR,
      plot_labels[3], plots_all[3],
      plot_labels[2], plots_all[2],
    )

    #ax_phiScal.set_title(file_to_use.split("_")[3:])
    plt.suptitle(file_to_use.split("_")[3:], fontsize=16, y=1.001)

    graph_data = []
    plot_labels = []
    plots_all = []

if __name__ == "__main__":

    from time import sleep
    from matplotlib import pyplot as plt
    from matplotlib import animation as animation
    from matplotlib import style
    from matplotlib.ticker import FormatStrFormatter
    from matplotlib import gridspec

    #sleep(30)

    style.use("seaborn-poster")
    #style.use("fivethirtyeight")

    gs = gridspec.GridSpec(nrows=1,ncols=2)

    fig = plt.figure()

    fig.set_tight_layout(True)

    ax_phiScal = fig.add_subplot(gs[0,0])
    ax_mR = fig.add_subplot(gs[0,1])

    ani_live = animation.FuncAnimation(
      fig = fig, func = results_live, interval = 1000, save_count = None
    )

    plt.show()
