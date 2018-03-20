#!/usr/bin/env python

file_path = "/home/dimitar/projects/STT_theories/results"
file_results = "STT_phiScal_I_"
file_live = "live_plot_solver"

from itertools import cycle
colors = [ "b", "g", "r", "c", "m", "y" ]
colors_index = [ i for i in range(0,len(colors)+1)]
colors_cycle = cycle(colors_index)

min_pressure = 1e-16

def get_recent_file(starts_with):
    import glob
    import os

    return max(glob.glob("{}*".format(starts_with)), key=os.path.getctime)

def data_loader():

    clear_plots_params(ax_Omega, "r", "Omega", 12, 10 )
    clear_plots_params(ax_Z, "r", "Z", 12, 10 )
    clear_plots_params(ax_J, "r", "J", 12, 10 )

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

        clear_plots_params(ax_Omega, "r", "Omega", 12, 10 )
        clear_plots_params(ax_Z, "r", "Z", 12, 10 )
        clear_plots_params(ax_J, "r", "J", 12, 10 )

        pos_sharp = 0

def clear_plots_params(ax, label_x, label_y, fontsize, ticksize):

    from matplotlib.ticker import FormatStrFormatter

    ax.clear()

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=fontsize)
    ax.xaxis.set_tick_params(labelsize=ticksize)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=fontsize)
    ax.yaxis.set_tick_params(labelsize=ticksize)

def my_plotting(ax,x,y,R,phiScal_inf, color):

    ax.plot(
      x,y,
      marker="o",markersize=5,
      linewidth=1.5,
      color=color
    )

    if x[0] < x[1] and R and R < 50:
        ax.axvline(x=R, alpha=0.4, linestyle = "--", linewidth=1.5)

    ax.axvline(x=phiScal_inf, alpha=0.4, linestyle="-", linewidth=1.5)

def update_J(frame, ax_Omega, ax_Z, ax_J ):

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

    c_i = next(colors_cycle)

    if to_be_ploted[0][0] > to_be_ploted[0][-1]:
        c_i -= 1

    if c_i == len(colors):
        c_i = next(colors_cycle)

        clear_plots_params(ax_Omega, "r", "Omega", 12, 10 )
        clear_plots_params(ax_Z, "r", "Z", 12, 10 )
        clear_plots_params(ax_J, "r", "J", 12, 10 )

    my_plotting(ax_Omega,to_be_ploted[0], to_be_ploted[8], plot_R, plot_phiScal_inf, colors[c_i])
    my_plotting(ax_Z,to_be_ploted[0], to_be_ploted[7], plot_R, plot_phiScal_inf, colors[c_i])
    my_plotting(ax_J,to_be_ploted[0], to_be_ploted[9], plot_R, plot_phiScal_inf, colors[c_i])

    return

def get_figure(nrows):

    from matplotlib import style
    from matplotlib.gridspec import GridSpec

    style.use("seaborn-poster")
    gs = GridSpec(nrows=nrows,ncols=1)

    fig = plt.figure()
    fig.set_tight_layout(True)

    axes = []
    for i in range(nrows):
        axes.append(fig.add_subplot(gs[i,0]))

    return fig, axes

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    global ax_Omega, ax_Z, ax_J

    fig_J, axs = get_figure(3)
    ax_Omega = axs[0]
    ax_Z = axs[1]
    ax_J = axs[2]

    ani_J = FuncAnimation(
      fig = fig_J, func = update_J,
      fargs = (ax_Omega, ax_Z, ax_J),
      frames = data_loader,
      interval = 500
    )

    plt.show()

