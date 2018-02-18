#!/usr/bin/env python

# Based on
#   https://pythonprogramming.net/live-graphs-matplotlib-tutorial/

file_path = "/home/dimitar/projects/STT_theories/results/"
file_name_results = "STT_phiScal_"
file_name_live = "live_plot_solver"

system_names = ["r", "phiScal", "Q", "p", "LambdaMetr", "m", "rho"]
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

def my_plotting_phiScal(
  ax,
  label_x, plots_x,
  label_y, plots_y,
  plots_R,
  plots_phiScal_inf
):

    from random import randint

    label_fontsize = 12
    ticks_label_size = 10

    ax.clear()

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=label_fontsize)
    ax.xaxis.set_tick_params(labelsize=ticks_label_size)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=label_fontsize)
    ax.yaxis.set_tick_params(labelsize=ticks_label_size)

    colors_iter = iter(colors)
    c = "w"

    try:
        for x, y, R, phiScal_inf in zip(plots_x, plots_y, plots_R, plots_phiScal_inf):

            if R:
                c = next(colors_iter)

            ax.plot(
                x, y[1],
                marker="o", markersize = 5,
                linewidth = 1.5,
                color = c
            )

            if R:
                ax.axvline(x=R, alpha = 0.4, linestyle = "--", linewidth = 1.5)

            ax.axvline(x=phiScal_inf, alpha = 0.4, linestyle = "-", linewidth = 1.5)
    except ValueError:
        pass

def my_plotting_rho(
  ax,
  label_x, plots_x,
  label_y, plots_y,
):
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
        for x, y in zip(plots_x, plots_y):

            ax.plot(
                x[:len(y[-1])], y[-1],
                marker="o", markersize = 5,
                linewidth = 1.5
            )
    except ValueError:
        pass

def animate(something):

    import glob
    import os

    all_files = glob.glob(file_path + file_name_live + "*")
    while not len(all_files):
        all_files = glob.glob(file_path + file_name_live + "*")

    file_to_use = max(all_files, key=os.path.getctime)

    with open(file_to_use,"r+") as f:
        graph_data = f.read()

    graph_data = [ k for k in graph_data.split("# ") if len(k) ]

    print("asdasda " + str(len(graph_data)))

    if len(graph_data) > 2:

        if not "R = 0.000000e+00" in graph_data[-1] and not "R = 0.000000e+00" in graph_data[-2]:
            graph_data = graph_data[-len(colors):]
        elif "R = 0.000000e+00" in graph_data[-1]:
            graph_data = graph_data[-2*len(colors):]
        else:
            del graph_data[-1]
            graph_data = graph_data[-2*len(colors):]

    print(str(len(graph_data)) + "\n\n")

    graph_data = [ k.split("\n") for k in graph_data ]

    plot_titles = [ k.pop(0) for k in graph_data ]

    plot_R = [ float(k.split(",")[1].split(" = ")[1]) for k in plot_titles ]
    plot_phiScal_inf = [ float(k.split(",")[2].split(" = ")[1]) for k in plot_titles ]

    plot_titles = [ k.split(",")[0] for k in plot_titles ]

    plot_x, plot_y = \
        [ [] for k in range(0,len(graph_data)) ], \
        [ [ [] for l in range(0, len(system_names)) ] for k in range(0,len(graph_data)) ]

    try:
        for single_set, single_set_x, single_set_y \
          in zip(graph_data, plot_x, plot_y):

            for single_line in single_set:
                if len(single_line) > 1:
                    tmp = [ j for j in single_line.split(" ") if len(j) > 1 ]

                    single_set_x.append(float(tmp[0]))

                    for m,n in zip(single_set_y, tmp):
                        m.append(float(n))

                    single_set_y[-1] = [ m for m in single_set_y[-1] if m ]
    except ValueError:
        pass

    my_plotting_phiScal(
      ax_phiScal,
      "r", plot_x,
      "phiScal", plot_y,
      plot_R,
      plot_phiScal_inf
    );

    my_plotting_rho(
      ax_rho,
      "r", plot_x,
      "rho", plot_y
    );

    plt.suptitle(file_to_use.split("_")[6:10], fontsize=16, y=1.001)

if __name__ == "__main__":

    from time import sleep
    from sys import argv
    from matplotlib import pyplot as plt
    from matplotlib import animation as animation
    from matplotlib import style
    from matplotlib.ticker import FormatStrFormatter

    #sleep(1)

    style.use("seaborn-poster")
    #style.use("fivethirtyeight")

    fig, axes = plt.subplots(1,2)

    fig.set_tight_layout(True)

    ax_phiScal = axes[0]
    ax_rho = axes[1]

    ani_live = animation.FuncAnimation(fig,animate,interval=100)

    plt.show()
