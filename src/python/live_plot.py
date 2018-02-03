#!/usr/bin/env python

# Based on
#   https://pythonprogramming.net/live-graphs-matplotlib-tutorial/

file_path = "/home/dimitar/projects/STT_theories/results/"
file_name_results = "LogisticEq_"
file_name_live = "live_plot_"

def animate_live(i):

    import glob
    import os

    all_files = glob.glob(file_path + file_name_live + "*")

    file_to_use = max(all_files, key=os.path.getctime)

    with open(file_to_use,'r') as f:
        graph_data = f.read()

    graph_data = [ k.split("\n") for k in graph_data.split("# ") if len(k) > 1 ]

    plot_x, plot_y = \
        [ [] for k in range(0,len(graph_data)) ], \
        [ [] for k in range(0,len(graph_data)) ]

    plot_titles = [ k.pop(0) for k in graph_data ]

    for single_set, single_set_x, single_set_y \
      in zip(graph_data, plot_x, plot_y):

        for single_line in single_set:
            if len(single_line) > 1:
                x,y = single_line.split(",")
                single_set_x.append(float(x))
                single_set_y.append(float(y))

    ax_live.clear()
    ax_live.set_title(file_to_use.split("_")[-1])
    ax_live.xaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    ax_live.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    for single_set_x, single_set_y , single_set_title \
      in zip(plot_x, plot_y, plot_titles):
        ax_live.plot(\
          single_set_x, single_set_y, \
          marker="o", markersize=7, \
          linewidth=1.5, \
          label=single_set_title \
        )

    #ax_live.legend(loc="best")
    #ax_live.legend(loc="right")

def animate_result(i):

    import glob
    import os

    all_files = glob.glob(file_path + file_name_results + "*")

    file_to_use = max(all_files, key=os.path.getctime)

    with open(file_to_use,'r') as f:
        graph_data = f.read()

    only_test_const = 87
    graph_data = graph_data.split("\n")[only_test_const:]

    plot_x, plot_y = [], []

    for single_set_points in graph_data:
        if len(single_set_points) > 1:
            x,y = single_set_points.split(",")
            plot_x.append(float(x))
            plot_y.append(float(y))

    ax_result.clear()
    ax_result.set_title(file_to_use.split("_")[-1])
    ax_result.xaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    ax_result.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    ax_result.plot( plot_x, plot_y, marker="o", markersize=7, linewidth=1.5 )

    #ax_result.legend(loc="best")

if __name__ == "__main__":

    from time import sleep
    from sys import argv
    from matplotlib import pyplot as plt
    from matplotlib import animation as animation
    from matplotlib import style
    from matplotlib.ticker import FormatStrFormatter

    sleep(1)

    style.use("seaborn-poster")
    #style.use("fivethirtyeight")

    fig_live = plt.figure()
    ax_live = fig_live.add_subplot(1,1,1)

    fig_result = plt.figure()
    ax_result = fig_result.add_subplot(1,1,1)

    ani_live = animation.FuncAnimation(fig_live,animate_live,interval=1000)
    ani_result = animation.FuncAnimation(fig_result,animate_result,interval=1000)

    plt.show()
