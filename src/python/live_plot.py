#!/usr/bin/env python

# Based on
#   https://pythonprogramming.net/live-graphs-matplotlib-tutorial/

def animate(i):

    file_path = "../../results/"
    file_name = "live_plot_data"

    with open(file_path + file_name,'r') as f:
        graph_data = f.read()

    plot_x, plot_y = \
        [ [] for k in range(0,len(graph_data.split("\n\n"))) ], \
        [ [] for k in range(0,len(graph_data.split("\n\n"))) ]

    for single_set, single_set_x, single_set_y \
      in zip(graph_data.split("\n\n"),plot_x,plot_y):

        for single_line in single_set.split("\n"):
            if len(single_line) > 1:
                x,y = single_line.split(",")
                single_set_x.append(float(x))
                single_set_y.append(float(y))

    ax.clear()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    i = 0
    for single_set_x, single_set_y in zip(plot_x, plot_y):
        ax.plot(\
          single_set_x, single_set_y, \
          marker="o", markersize=7, \
          linewidth=1.5, \
          label=str(i) \
        )

        i+=1

    ax.legend(loc="best")

if __name__ == "__main__":

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






