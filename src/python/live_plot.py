#!/usr/bin/env python

# Based on
#   https://pythonprogramming.net/live-graphs-matplotlib-tutorial/

file_path = "/home/dimitar/projects/STT_theories/results/"
file_name = "LogisticEq_live_plot"

def animate(i):

    with open(file_path + file_name,'r') as f:
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

    ax.clear()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    for single_set_x, single_set_y , single_set_title \
      in zip(plot_x, plot_y, plot_titles):
        ax.plot(\
          single_set_x, single_set_y, \
          marker="o", markersize=7, \
          linewidth=1.5, \
          label=single_set_title \
        )

    ax.legend(loc="best")

if __name__ == "__main__":

    from sys import argv

    from matplotlib import pyplot as plt
    from matplotlib import animation as animation
    from matplotlib import style
    from matplotlib.ticker import FormatStrFormatter

    style.use("seaborn-poster")
    #style.use("fivethirtyeight")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    ani = animation.FuncAnimation(fig,animate,interval=1000)
    plt.show()






