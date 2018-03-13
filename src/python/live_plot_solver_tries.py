#!/usr/bin/env python

file_path = "/home/dimitar/projects/STT_theories/results"
file_results = "STT_phiScal_I_"
file_live = "live_plot_solver"

from itertools import cycle
colors = [ "b", "g", "r", "c", "m", "y" ]
colors_vals = [ i for i in range(0,len(colors)+1)]
colors_cycle = cycle(colors_vals)

min_pressure = 1e-16

def get_recent_file(starts_with):
    import glob
    import os

    return max(glob.glob("{}*".format(starts_with)), key=os.path.getctime)

def clear_plots_params(ax, label_x, label_y, fontsize, ticksize):

    from matplotlib.ticker import FormatStrFormatter

    ax.clear()

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=fontsize)
    ax.xaxis.set_tick_params(labelsize=ticksize)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=fontsize)
    ax.yaxis.set_tick_params(labelsize=ticksize)

    return

def data_loader():

    data_chunk = []
    pos_sharp = 0
    while True:

        file_model = file_path + "/" + file_live

        filename = get_recent_file(file_model)

        while get_recent_file(file_model) == filename:

            data_chunk.clear()
            with open(filename, "r") as f:
                for cnt, line in enumerate(f.readlines()[pos_sharp:]):
                    data_chunk.append(line.strip())
                    if "#" in line and line.strip() not in data_chunk[0]:
                        pos_sharp += cnt
                        break

            if len(data_chunk):
                if "#" in data_chunk[-1]:
                    data_chunk.pop(-1)
                yield data_chunk
            else:
                pos_sharp = 0

        clear_plots_params(ax_phiScal, "r", "phiScal", 12, 10 )
        clear_plots_params(ax_p, "r", "p", 12, 10 )
        clear_plots_params(ax_rho, "r", "rho", 12, 10 )
        pos_sharp = 0

def get_figure():

    from matplotlib import style
    from matplotlib.gridspec import GridSpec

    style.use("seaborn-poster")
    gs = GridSpec(nrows=2,ncols=2)

    fig = plt.figure()
    fig.set_tight_layout(True)

    ax_phiScal = fig.add_subplot(gs[:,0])
    ax_p = fig.add_subplot(gs[0,1])
    ax_rho= fig.add_subplot(gs[1,1])

    clear_plots_params(ax_phiScal, "r", "phiScal", 12, 10 )
    clear_plots_params(ax_p, "r", "p", 12, 10 )
    clear_plots_params(ax_rho, "r", "rho", 12, 10 )

    return fig, ax_phiScal, ax_p, ax_rho

def my_plotting_phiScal(ax, x, y, R, phiScal_inf, color):

    try:
        ax.plot(
          x,y,
          marker="o",markersize=5,
          linewidth=1.5,
          color=color
        )
        if R and R < 50:
            ax.axvline(x=R, alpha=0.4, linestyle = "--", linewidth=1.5)

        ax.axvline(x=phiScal_inf, alpha=0.4, linestyle="-", linewidth=1.5)

    except ValueError:
        pass

    return

def my_plotting_rho(ax, x, y, color):

    try:
        ax.plot(
          x,y,
          marker="o",markersize=5,
          linewidth=1.5,
          color=color
        )

    except ValueError:
        pass

    return

def my_plotting_p(ax, x, y, color):

    try:
        ax.plot(
          x,y,
          marker="o",markersize=5,
          linewidth=1.5,
          color=color
        )

        ax.axhline(y=1e-14, linestyle="--", alpha=0.5, linewidth=2)
    except ValueError:
        pass

    return

def update(frame, ax_phiScal, ax_p, ax_rho ):

    if len(frame) == 0:
        return

    titles = frame.pop(0)
    try:
        plot_R = float(titles.split(",")[1].split(" = ")[1])
        plot_phiScal_inf = float(titles.split(",")[2].split(" = ")[1])
    except IndexError:
        return

    to_be_ploted = [ [] for i in frame[0].split(" ") if len(i) ]

    for line in frame:
        if len(to_be_ploted) == len(line.split(" ")):
            for m, n in zip(to_be_ploted, line.split(" ")):
                try:
                    m.append(float(n))
                except ValueError:
                    return
        else:
            return

    to_be_ploted[-1] = [ i for i in to_be_ploted[-1] if i ]
    to_be_ploted[3] = [ i for i in to_be_ploted[3] if i > min_pressure ]

    c = next(colors_cycle)

    if to_be_ploted[0][0] > to_be_ploted[0][-1]:
        c -= 1

    if c == len(colors):
        c = next(colors_cycle)

        clear_plots_params(ax_phiScal, "r", "phiScal", 12, 10 )
        clear_plots_params(ax_p, "r", "p", 12, 10 )
        clear_plots_params(ax_rho, "r", "rho", 12, 10 )

    my_plotting_phiScal(
      ax_phiScal, to_be_ploted[0], to_be_ploted[1], plot_R, plot_phiScal_inf, colors[c]
    )

    if to_be_ploted[0][0] < to_be_ploted[0][-1]:

        my_plotting_p(
          ax_p, to_be_ploted[0][:len(to_be_ploted[3])], to_be_ploted[3], colors[c]
        )

        my_plotting_rho(
          ax_rho, to_be_ploted[0][:len(to_be_ploted[3])], to_be_ploted[-1][:len(to_be_ploted[3])], colors[c]
        )

    return

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    global ax_phiScal, ax_p, ax_rho

    fig, ax_phiScal, ax_p, ax_rho = get_figure()

    clear_plots_params(ax_phiScal, "r", "phiScal", 12, 10 )
    clear_plots_params(ax_p, "r", "p", 12, 10 )
    clear_plots_params(ax_rho, "r", "rho", 12, 10 )

    animation_live = FuncAnimation(
      fig = fig, func = update, fargs = (ax_phiScal, ax_p, ax_rho),
      frames = data_loader, interval = 10, save_count = 0
    )

    plt.show()
