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
  label_x,
  label_y,
  all_plots,
  plots_R,
  plots_phiScal_inf
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

    colors_iter = iter(colors)
    c = "w"

    try:
        for single_plot, R, phiScal_inf in zip(all_plots, plots_R, plots_phiScal_inf):

            if single_plot[0][0] < single_plot[0][-1]:
                c = next(colors_iter)

            ax.plot(
                single_plot[0],
                single_plot[1],
                marker="o",
                markersize = 5,
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
  label_x,
  label_y,
  all_plots
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
        for single_plot in all_plots:

            ax.plot(
                single_plot[0][:len(single_plot[-1])],
                single_plot[-1],
                marker="o",
                markersize = 5,
                linewidth = 1.5
            )
    except ValueError:
        pass

def my_plotting_p(
  ax,
  label_x,
  label_y,
  all_plots
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
        for single_plot in all_plots:

            ax.plot(
                single_plot[0][:len(single_plot[3])],
                single_plot[3],
                marker="o",
                markersize = 5,
                linewidth = 1.5
            )

            ax.axhline(y=1e-14, linestyle="--", alpha=0.5, linewidth=2)
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
    graph_data = graph_data[-len(colors):]

    graph_data = [ k.split("\n") for k in graph_data ]

    plot_titles = [ k.pop(0) for k in graph_data ]

    plot_R = [ float(k.split(",")[1].split(" = ")[1]) for k in plot_titles ]
    plot_phiScal_inf = [ float(k.split(",")[2].split(" = ")[1]) for k in plot_titles ]

    plot_titles = [ k.split(",")[0] for k in plot_titles ]

    all_plots = [
      [ [] for l in range(0, len(system_names)) ] for k in range(0,len(graph_data))
    ]

    try:
        for single_plot, single_data in zip(all_plots, graph_data):
            for single_line in single_data:
                if len(single_line) > 1:
                    tmp = [ k for k in single_line.split(" ") if len(k) > 1 ]

                for m,n in zip(single_plot, tmp):
                    m.append(float(n))

                single_plot[-1] = [ k for k in single_plot[-1] if k ]
                single_plot[3] = [ k for k in single_plot[3] if k > 1e-12 ]

    except ValueError:
        pass

    my_plotting_phiScal(
      ax_phiScal,
      system_names[0],
      system_names[1],
      all_plots,
      plot_R,
      plot_phiScal_inf
    );

    my_plotting_rho(
      ax_rho,
      system_names[0],
      system_names[-1],
      all_plots
    );

    my_plotting_p(
      ax_p,
      system_names[0],
      system_names[3],
      all_plots
    );

    plt.suptitle(file_to_use.split("_")[6:10], fontsize=16, y=1.001)

if __name__ == "__main__":

    from time import sleep
    from sys import argv
    from matplotlib import pyplot as plt
    from matplotlib import animation as animation
    from matplotlib import style
    from matplotlib.ticker import FormatStrFormatter
    from matplotlib import gridspec

    #sleep(1)

    style.use("seaborn-poster")
    #style.use("fivethirtyeight")

    gs = gridspec.GridSpec(2,2)

    fig = plt.figure()

    fig.set_tight_layout(True)

    ax_phiScal = fig.add_subplot(gs[:,0])
    ax_p = fig.add_subplot(gs[0,1])
    ax_rho= fig.add_subplot(gs[1,1])

    ani_live = animation.FuncAnimation(fig,animate,interval=100)

    plt.show()
