#!/usr/bin/env python

# Based on
#   https://pythonprogramming.net/live-graphs-matplotlib-tutorial/

def animate(i):

    file_path = "../../results/"
    file_name = "live_plot_data"

    with open(file_path + file_name,'r') as f:
        graph_data = f.read()

    plot_x, plot_y = \
        [ [] for i in range(0,len(graph_data.split("\n\n"))) ], \
        [ [] for i in range(0,len(graph_data.split("\n\n"))) ]

    for single_set, single_set_x, single_set_y \
      in zip(graph_data.split("\n\n"),plot_x,plot_y):

        for i in single_set.split("\n"):
            if len(i) > 1:
                x,y = i.split(",")
                single_set_x.append(x)
                single_set_y.append(y)

    ax1.clear()
    for i, j in zip(plot_x, plot_y):
        ax1.plot(i,j,marker="o")

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    from matplotlib import style
    from matplotlib.ticker import FormatStrFormatter

    #style.use("seaborn-poster")
    style.use("fivethirtyeight")

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.3e'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.3e'))

    ani = animation.FuncAnimation(fig,animate,interval=1000)
    plt.show()






