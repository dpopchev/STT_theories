#!/usr/bin/env python

# generator to yield chunks of file whos content is numbers separated by the character "#"
# it achieves it by loading only the unread part of the file
# using f.readlines()[ <line of last encountered special character:]
def read_line_by_line():
    data_chunk = []
    pos_sharp = 0
    while True:
        data_chunk.clear()
        with open("live_plot_solver_STT_phiScal_I_AkmalPR_beta-1.000e+01_m0.000e+00_lambda0.000e+00","r") as f:
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

# similar to the above generator this does it by keeping track of the offset
def read_by_offset():
    data_chunk = []
    pos = 0
    pos_sharp = 0
    while True:
        data_chunk.clear()
        with open("live_plot_solver_STT_phiScal_I_AkmalPR_beta-1.000e+01_m0.000e+00_lambda0.000e+00","r") as f:
            if pos >= f.seek(0,2):
                pos = 0
                if pos_sharp >= f.seek(0,2):
                    pos_sharp = 0
            f.seek(pos_sharp)
            line = f.readline()
            pos += len(line)
            data_chunk.append(line.strip())
            for line in f.readlines():
                if "#" not in line:
                    data_chunk.append(line.strip())
                    pos += len(line)
                else:
                    pos_sharp = pos
                    break

        yield data_chunk

# simple function to plot calculate the units coefficients
# it an dictionary for the result
def units_coef_clac():
    # mas of sun in kg
    const_msun = 1.9891e30
    # gravitational const in m^3kg^-1s^-2
    const_g = 6.67384e-11
    # speed of light in ms^-1
    const_c = 299792458

    units = {}

    # units of density in g cm^-3
    units["density"] = 1e-3 * const_c**6/(const_g**3*const_msun**2)

    # units of pressure in dyns cm^-3
    units["pressure"] = const_c**8/(const_g**3*const_msun**2)*10

    # units of rad coordinate in km
    units["rad"] = 1e-3 * const_g * const_msun/const_c**2

    # units of moment of inertia
    units["j"] = 1e7 * const_g**2 * const_msun**3/const_c**4

    return units

# load data from my results file and return dictionary with "labels" and "data"
def load_myResFile(filename):

    with open(filename, "r") as f:
        all_data = f.read().split("#")[1]

    all_data = [ i.strip() for i in all_data.split("\n") if len(i) ]

    labels = [ i for i in all_data.pop(0).split(" ") if len(i) ]

    data = [ [] for i in labels ]

    for line in all_data:
        for m, n in zip(all_data, line.split(" ")):
            m.append(float(n))

    return { "labels" : labels, "data" : data }

# load data from Daniela results file and return dictionary with "labels" and "data"
# they also undergo conversion as they are not dimensionless
def load_Daniela_ResFile(filename):

    with open(filename, "r") as f:
        all_data = f.readlines()

    for i, line in enumerate(all_data):
        if "rho" in line:
            all_data.pop(i)

    all_data = [ i.strip() for i in all_data if len(i) ]

    labels = [
        "rho" , "R", "M", "J", "phiScal_c", "p_c"
    ]

    # i am not interested in all data, so only the indexes of the corresponding
    # labels listed above will come here
    indexes =[
        0, 1, 2, 3, 4, 10
    ]

    if len(labels) != len(indexes):
        print("LOL 1")
        exit()

    units = units_coef_clac()

    all_units = [
        units["density"], units["rad"], 1, units["j"], 1, 1
    ]

    if len(all_units) != len(indexes):
        print("LOL 2")
        exit()

    data = [ [] for i in labels ]

    for cnt, line in enumerate(all_data):
        for i, c, d in zip(indexes, all_units, data):
            try:
                d.append(float(line.split(" ")[i])/c)
            except ValueError:
                print("\n ValErr: line {} : {} \n".format(cnt, line))
                break

    return { "labels" : labels, "data" : data }

# fast plotting a line on axes
def fast_plot_line(ax, x, y):

    ax.plot(
        x, y,
        linewidth = 1.5
    )

    return

# fast plotting markers on axes
def fast_plot_markers(ax, x, y):

    ax.plot(
        x, y,
        marker = "o",
        markersize = 3.5,
        alpha = 0.4
    )

    return

# set labels, fontsize and ticksize on the axis
def set_params(ax, label_x, label_y, fontsize, ticksize):

    from matplotlib.ticker import FormatStrFormatter

    ax.clear()

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=fontsize)
    ax.xaxis.set_tick_params(labelsize=ticksize)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=fontsize)
    ax.yaxis.set_tick_params(labelsize=ticksize)

    return

# create figure, attach axes, return them
def fast_get_figure():

    from matplotlib import style
    from matplotlib.gridspec import GridSpec

    style.use("seaborn-poster")
    gs = GridSpec(nrows=1,ncols=1)

    fig = plt.figure()
    fig.set_tight_layout(True)

    ax = fig.add_subplot(gs[:])

    return fig, ax

class compare_graphs:

    def __init__(self):
        self.mine_data = []
        self.mine_labels = []

        self.Daniela_data = []
        self.Daniela_labels = []

    def print_info_mapping(self):

        mapping = {
            0 : "p_c",
            1 : "varPhi_c",
            4 : "M",
            5 : "AR",
            6 : "rho_c",
            7 : "J"
        }

        for i in mapping.keys():
            print("{} -> {}".format(i, mapping[i]))

        return

    @staticmethod
    def units_coef_clac():
        # mas of sun in kg
        const_msun = 1.9891e30
        # gravitational const in m^3kg^-1s^-2
        const_g = 6.67384e-11
        # speed of light in ms^-1
        const_c = 299792458

        units = {}

        # units of density in g cm^-3
        units["density"] = 1e-3 * const_c**6/(const_g**3*const_msun**2)

        # units of pressure in dyns cm^-3
        units["pressure"] = const_c**8/(const_g**3*const_msun**2)*10

        # units of rad coordinate in km
        units["rad"] = 1e-3 * const_g * const_msun/const_c**2

        # units of moment of inertia
        units["j"] = 1e7 * const_g**2 * const_msun**3/const_c**4

        return units

    def load_Daniela_ResFile(self,filename):

        print("GO TO THE FOLDER or will not now where is the central pressure !!!!!")

        with open(filename, "r") as f:
            all_data = f.readlines()

        for i, line in enumerate(all_data):
            if "rho" in line:
                all_data.pop(i)

        all_data = [ i.strip() for i in all_data if len(i.strip()) ]
        labels = [
            "rho" , "R", "M", "J", "phiScal_c", "p_c"
        ]

        # i am not interested in all data, so only the indexes of the corresponding
        # labels listed above will come here
        # in some the central pressure is in different column
        if filename in [ "models_APR_beta-4.5_mphi5e-3.dat", "models_APR_beta-6.0_mphi0.dat", "models_APR_beta-6.0_mphi1e-3.dat"]:
            indexes =[
              0, 1, 2, 3, 4, 10
            ]
        else:
            indexes =[
              0, 1, 2, 3, 4, 11
            ]

        if len(labels) != len(indexes):
            print("LOL 1")
            exit()

        units = self.units_coef_clac()

        all_units = [
            units["density"], units["rad"], 1, units["j"], 1, 1
        ]

        if len(all_units) != len(indexes):
            print("LOL 2")
            exit()

        data = [ [] for i in labels ]

        for cnt, line in enumerate(all_data):
            for i, c, d in zip(indexes, all_units, data):
                try:
                    d.append(float(line.split(" ")[i])/c)
                except ValueError:
                    print("\n ValErr: line {} : {} \n".format(cnt, line))
                    break

        self.Daniela_data.append(data)
        self.Daniela_labels.append(labels)

    def load_myResFile(self, filename):

        with open(filename, "r") as f:
            all_data = f.read().split("#")[1]

        all_data = [ i.strip() for i in all_data.split("\n") if len(i.strip()) ]

        labels = [ i.strip() for i in all_data.pop(0).split(" ") if len(i.strip()) ]

        data = [ [] for i in labels ]

        for line in all_data:
            for m, n in zip(data, line.split(" ")):
                m.append(float(n))

        for index, val in enumerate(data[1]):
            if val < 0:
                data[1][index] *= (-1)

        self.mine_data.append(data)
        self.mine_labels.append(labels)

    @staticmethod
    def fast_get_figure(self):

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=1,ncols=1)

        self.fig = plt.figure()
        self.fig.set_tight_layout(True)

        self.ax = self.fig.add_subplot(gs[:])

    @staticmethod
    def fast_plot(self, index_x = 0, index_y = 0):

        if len(self.mine_data):
            for x, y in zip(self.mine_data, self.mine_data):
                self.ax.plot(
                    x[index_x], y[index_y],
                    linewidth = 1.5
                )

        # mapping mine indexes over hers
        Daniela_indexes = [ 5, 4, -10, -10, 2, 1, 0, 3 ]

        if len(self.Daniela_data):
            for x, y in zip(self.Daniela_data, self.Daniela_data):
                self.ax.plot(
                    x[Daniela_indexes[index_x]],
                    y[Daniela_indexes[index_y]],
                    marker = "o",
                    markersize = 5,
                    alpha = 0.4,
                    linewidth = 0
                )

    @staticmethod
    def set_params(self):

        from matplotlib.ticker import FormatStrFormatter

        fontsize = 12
        ticksize = 10

        self.ax.clear()

        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        self.ax.set_xlabel("x",fontsize=fontsize)
        self.ax.xaxis.set_tick_params(labelsize=ticksize)

        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        self.ax.set_ylabel("y",fontsize=fontsize)
        self.ax.yaxis.set_tick_params(labelsize=ticksize)

    def clear_mine(self):
        self.mine_data.clear()
        self.mine_labels.clear()

    def clear_Daniela(self):
        self.Daniela_data.clear()
        self.Daniela_labels.clear()

    def plot_it(self, index_x = 0, index_y = 0):

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        self.fast_get_figure(self)
        self.set_params(self)
        self.fast_plot(self,index_x, index_y)

        plt.show()
