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

    from time import sleep

    clear_plots_params(ax_phiScal_left, "r", "phiScal", 12, 10 )
    clear_plots_params(ax_Q_left, "r", "Q", 12, 10 )
    clear_plots_params(ax_phiScal_right, "r", "phiScal", 12, 10 )
    clear_plots_params(ax_Q_right, "r", "Q", 12, 10 )

    data_chunk = []
    pos_sharp = 0
    while True:

        file_model = file_path + "/" + file_live

        filename = get_recent_file(file_model)

        #last_modified = 0
        while get_recent_file(file_model) == filename:

            data_chunk.clear()
            with open(filename, "r") as f:
                #last_modified = Path(filename).stat().st_mtime
                for cnt, line in enumerate(f.readlines()[pos_sharp:]):
                    data_chunk.append(line.strip())
                    if "#" in line and line.strip() not in data_chunk[0]:
                        pos_sharp += cnt
                        break

            if len(data_chunk):
                if "#" in data_chunk[-1]:
                    data_chunk.pop(-1)
                yield data_chunk
                sleep(0.1)
            else:
                pos_sharp = 0

        #clear_plots_params(ax_phiScal, "r", "phiScal", 12, 10 )
        #clear_plots_params(ax_Q, "r", "Q", 12, 10 )

        pos_sharp = 0

def my_plotting(ax,x,y,R,phiScal_inf, color):

    ax.plot(
      x,y,
      marker="o",markersize=5,
      linewidth=1.5,
      color=color,
      label=color
    )

    if x[0] < x[1]:
        ax.set_title(
          "({:.3e}, {:.3e}) -> ({:.3e}, {:.3e})".format(
            x[0], y[0], x[-1], y[-1]
          ),
          fontsize=10
        )

        if R and R < 50:
            ax.axvline(x=R, alpha=0.4, linestyle = "--", linewidth=1.5, label="radii")

    else:
        ax.set_title(
          "({:.3e}, {:.3e}) -> ({:.3e}, {:.3e})".format(
            x[0], y[0], x[-1], y[-1]
          ),
          fontsize=10
        )

    ax.axvline(x=phiScal_inf, alpha=0.4, linestyle="-", linewidth=1.5, label = "fitting")

def update_phiScal(frame, ax_phiScal, ax_Q ):

    if len(frame) == 0:
        print("\n\t\t\t\t frame empty")
        return

    titles = frame.pop(0)
    try:
        plot_R = float(titles.split(",")[1].split(" = ")[1])
    except IndexError:
        print("\n\t\t\t\t IndexError python shit")
        return

    to_be_ploted = [ [] for i in frame[0].split(" ") if len(i) ]

    for line in frame:
        if len(to_be_ploted) == len(line.split(" ")):
            for m, n in zip(to_be_ploted, line.split(" ")):
                try:
                    m.append(float(n))
                except ValueError:
                    print("\n\t\t\t\t ValueError python shit")
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

        clear_plots_params(ax_phiScal_left, "r", "phiScal", 12, 10 )
        clear_plots_params(ax_Q_left, "r", "Q", 12, 10 )
        clear_plots_params(ax_phiScal_right, "r", "phiScal", 12, 10 )
        clear_plots_params(ax_Q_right, "r", "Q", 12, 10 )

    if to_be_ploted[0][0] < to_be_ploted[0][-1]:
        my_plotting(ax_phiScal_left,to_be_ploted[0], to_be_ploted[1], plot_R, plot_phiScal_inf, colors[c_i])
        my_plotting(ax_Q_left,to_be_ploted[0], to_be_ploted[2], plot_R, plot_phiScal_inf, colors[c_i])
    else:
        my_plotting(ax_phiScal_right,to_be_ploted[0], to_be_ploted[1], plot_R, plot_phiScal_inf, colors[c_i])
        my_plotting(ax_Q_right,to_be_ploted[0], to_be_ploted[2], plot_R, plot_phiScal_inf, colors[c_i])

    return

def clear_plots_params(ax, label_x, label_y, fontsize, ticksize):

    from matplotlib.ticker import FormatStrFormatter

    ax.clear()

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=fontsize)
    ax.xaxis.set_tick_params(labelsize=ticksize)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=fontsize)
    ax.yaxis.set_tick_params(labelsize=ticksize)

def get_figure(nrows):

    from matplotlib import style
    from matplotlib.gridspec import GridSpec

    style.use("seaborn-poster")
    gs = GridSpec(nrows=nrows,ncols=1)

    fig = plt.figure()
    fig.set_tight_layout(True)

    return fig, fig.add_subplot(gs[0,0]), fig.add_subplot(gs[1,0])

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    global ax_phiScal, ax_Q

    fig_phiScal, ax_phiScal, ax_Q = get_figure(2)

    ani_phiscal = FuncAnimation(
      fig = fig_phiScal, func = update_phiScal,
      fargs = (ax_phiScal, ax_Q),
      frames = data_loader,
      interval = 500
    )

    plt.show()

