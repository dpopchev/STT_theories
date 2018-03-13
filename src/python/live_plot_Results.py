#!/usr/bin/env python

file_path = "/home/dimitar/projects/STT_theories/results"
file_results = "STT_phiScal_I_"

def get_recent_file(starts_with):
    import glob
    import os

    return max(glob.glob("{}*".format(starts_with)), key=os.path.getctime)

def data_loader():

    data_chunk = []
    while True:

        file_model = file_path + "/" + file_results

        filename = get_recent_file(file_model)

        while get_recent_file(file_model) == filename:

            data_chunk.clear()
            with open(filename, "r") as f:
                data_chunk = f.read().split("#")[1].split("\n")

            clear_plots_params(ax_phiScal_c, "p_c", "phiScal_c", 12, 10 )
            clear_plots_params(ax_M_AR, "AR", "M", 12, 10 )
            clear_plots_params(ax_J_M, "M", "J", 12, 10 )

            yield data_chunk, filename

def my_plotting(ax, x, y):

    try:
        ax.plot(
          x, y,
          marker="o",markersize=5,
          linewidth=1.5
        )
    except ValueError:
        pass

def update(
  frame, ax_phiScal_c, ax_M_AR, ax_J_M
):
    if len(frame[0]) == 0:
        return

    titles = [ i for i in frame[0].pop(0).split(" ") if len(i) ]

    to_be_ploted = [ [] for i in titles if len(i) ]

    for line in frame[0]:
        for m, n in zip(to_be_ploted, line.split(" ")):
            try:
                m.append(float(n))
            except ValueError:
                pass

    my_plotting(ax_phiScal_c, to_be_ploted[0], to_be_ploted[1])
    my_plotting(ax_M_AR, to_be_ploted[5], to_be_ploted[4])
    my_plotting(ax_J_M, to_be_ploted[4], to_be_ploted[7])

    plt.suptitle(frame[1].split("_")[4:], y=0.998)

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

def get_figure():

    from matplotlib import style
    from matplotlib.gridspec import GridSpec

    style.use("seaborn-poster")
    gs = GridSpec(nrows=2,ncols=2)

    fig = plt.figure()
    fig.set_tight_layout(True)

    ax_phiScal_c = fig.add_subplot(gs[0,:])
    ax_M_AR = fig.add_subplot(gs[1,0])
    ax_J_M = fig.add_subplot(gs[1,1])

    clear_plots_params(ax_phiScal_c, "p_c", "phiScal_c", 12, 10 )
    clear_plots_params(ax_M_AR, "AR", "M", 12, 10 )
    clear_plots_params(ax_J_M, "M", "J", 12, 10 )

    return fig, ax_phiScal_c, ax_M_AR, ax_J_M

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    global \
      ax_phiScal_c, ax_M_AR, ax_J_M

    fig, ax_phiScal_c, ax_M_AR, ax_J_M = get_figure()

    animation_live = FuncAnimation(
      fig = fig, func = update,
      fargs = (ax_phiScal_c, ax_M_AR, ax_J_M),
      frames = data_loader, interval = 1000, save_count = 0
    )

    plt.show()
