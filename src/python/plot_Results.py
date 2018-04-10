#!/usr/bin/env python


class ClassToFitMyNeeds:

    def __init__(self):

        # path to my folder with results
        self.path = ""

        # the model of a result filename and place to save the current
        self.model_res = ""
        self.fname_res = ""

        # place to append the results data, headline and labels
        self.data_res = []
        self.headline_res = []
        self.label_res = []

        self.units = self.units_coef_clac()

        # set the defaults values for the vars here
        self.set_defaults(self)

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
        units["density"] = 1e-3 * const_c**6 / (const_g**3 * const_msun**2)

        # units of pressure in dyns cm^-3
        units["pressure"] = const_c**8 / (const_g**3 * const_msun**2) * 10

        # units of rad coordinate in km
        units["rad"] = 1e-3 * const_g * const_msun / const_c**2

        # units of moment of inertia
        units["j"] = 1e7 * const_g**2 * const_msun**3 / const_c**4

        return units

    @staticmethod
    def set_defaults(self):

        self.path = "/home/dimitar/projects/STT_theories/results"

        #self.path_daniela =
            #"/home/dimitar/Documents/Teaching_Materials/University/"
            #+ "Uvod.Fizika.Cherni.Dupki/Doktorant/"
            #+ "Daniela_Static_NS_Massive_SlowRot"

        self.model_res = "STT_phiScal_"
        #self.model_live = "live_plot_phiScal_"

    @staticmethod
    def check_if_NotBlankString(_string):
        return bool(_string and _string.strip())

    # get the latest file from the _path with provided filename model in _model
    @staticmethod
    def get_latest_file(_path, _model):

        import glob
        import os

        if self.check_if_NotBlank(_path) and self.check_if_NotBlank(_model):

            filename = \
                max(
                    glob.glob(_path + "/" + _model + "*"),
                    key=os.path.getctime
                )
            filename = filename.split("/")[-1]

            print("\n latest file in \n\t: \n\t is {} \n".format(_path, filename))
        else:

            filename = ""

            print("\n path and name not set, some filename var!!! \n")

        return filename

    @staticmethod
    def load_file_result(_path, _filename):

        if not (self.check_if_NotBlank(_path) and
                self.check_if_NotBlank(_filename)):

            print("\n path or filename blank !! \n")

        with open(_path + "/" + _filename, "r") as f:
            all_data = f.readlines()

        self.headline_res.append(
            [
                i.strip() for i in
                all_data.pop(0).strip().split(" ")
                if
                "#" not in i and
                len(i.strip())
            ]
        )

        self.label_res.append(
            _filename.split("/")[-1][-41:]
        )

        self.data_res.append(
            [
                [] for i in all_data[0].strip().split(" ") if len(i.strip())
            ]
        )

        for cnt, line in enumerate(all_data):
            for apnd, num in zip(
                self.data_res[-1], [
                    float(i) for i in line.strip().split(" ") if len(i.strip())
                ]
            ):
                apnd.append(num)

        return

    @staticmethod
    def single_grid_placement(self, fig, gs):

        all_ax = []
        all_ax.append(fig.add_subplot(gs[:]))

        return all_ax

    @staticmethod
    def get_figure(self, nrows, ncols, grid_placement):

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=nrows, ncols=ncols)

        fig = plt.figure()
        fig.set_tight_layout(True)

        return fig, grid_placement(self, fig, gs)

    @staticmethod
    def get_index(max_index):

        while True:
            try:
                index = int(input("\t index =... "))
            except ValueError:
                print("\n Not a number is {} \n".format(index))
                continue

            if index < 0 or index >= max_index:
                print("\n {} not in allowed boundaries \n".format(index))
            else:
                break

        return index

    @staticmethod
    def get_xy_current_file(_headline):

        print("\n choose columns to plot \n")
        for cnt, value in enumerate(_headline[-1]):
            print("\t {}: {}".format(cnt, value))

        print("\n x.. \n")
        index_x = self.get_index(len(_headline[-1]))

        print("\n y.. \n")
        index_y = self.get_index(len(_headline[-1]))

        return index_x, index_y

    def plot_latest_ResultFile(self):

        self.fname_res = self.get_latest_file(self.path, self.model_res)

        print("\n will plot {}/{} \n".format(self.path, self.fname_res))

        self.load_file_result(self.path, self.fname_res)

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        index_x, index_y = self.get_xy_current_file(self.headline_res[-1])

        fig, all_ax = self.get_figure(self, 1, 1, self.single_grid_placement)

        ax = all_ax[-1]

        self.set_parms(
            self, ax,
            self.headline_res[-1][index_x],
            self.headline_res[-1][index_y]
        )


class PlotResults:

    def __init__(self):

        # absolute file path and name of the current file to be plotted
        self.path = ""
        self.model_live = ""
        self.model_result = ""
        self.filename_result = ""
        self.filename_live = ""

        # absolute file path and name for daniela results
        self.path_daniela = ""
        self.model_fname_daniela = ""
        self.filename_daniela = ""

        # which index corresponds to what value
        # depends on what type of file I want to plot
        self.mapping = []

        # the data itself
        self.data = []
        self.headline = []
        self.label = []

        self.data_daniela = []

        self.units = self.units_coef_clac(self)

        self.set_defaults(self)

        return

    # reset all data
    def clear(self):

        self.path = ""
        self.model_live = ""
        self.model_result = ""
        self.filename_result = ""
        self.filename_live

        self.path_daniela = ""
        self.model_fname_daniela = ""
        self.filename_daniela = ""

        self.mapping.clear()

        self.data.clear()
        self.headline.clear()
        self.label.clear()

        self.data_daniela.clear()

        self.units = self.units_coef_clac(self)

        self.set_defaults(self)

        return

    def set_defaults(self):

        self.path = "/home/dimitar/projects/STT_theories/results"

        self.path_daniela = \
            "/home/dimitar/Documents/Teaching_Materials/University/" \
            + "Uvod.Fizika.Cherni.Dupki/Doktorant/" \
            + "Daniela_Static_NS_Massive_SlowRot"

        self.model_result = "STT_phiScal_"
        self.model_live = "live_plot_phiScal_"

    @staticmethod
    def check_if_NotBlank(_string):
        return bool(_string and _string.strip())

    # get the latest file from the folder
    # by giving some path and model of a filename to search
    @staticmethod
    def get_latest_file(_path, _model):

        import glob
        import os

        # if the path and model name are not empty get the latest file
        # else complain and not do anything
        if self.check_if_NotBlank(_path) and
            self.check_if_NotBlank(_model):

            filename = \
                max(
                    glob.glob(_path + "/" + _model + "*"),
                    key=os.path.getctime
                )
            filename = filename.split("/")[-1]
            print("\n latest file: \n\t {} \n".format(filename))
        else:

            filename = ""

            print("\n path and name not set, some filename var!!! \n")

        return filename

    @staticmethod
    def get_index(max_index):

        while True:
            try:
                index = int(input("\t index =... "))
            except ValueError:
                print("\n Not a number is {} \n".format(index))
                continue

            if index < 0 or index >= max_index:
                print("\n {} not in allowed boundaries \n".format(index))
            else:
                break

        return index

    # get file chosen by the user
    @staticmethod
    def get_file(_path):

        import os

        list_dir = os.listdir(_path)

        for ind, val in enumerate(list_dir):
            print("{}: {}".format(ind, val))

        return list_dir[self.get_index(len(list_dir))]

    # load my results
    @staticmethod
    def load_file(self):

        with open(self.path + "/" + self.filename, "r") as f:
            all_data = f.readlines()

        self.headline.append(
            [
                i.strip() for i in
                all_data.pop(0).strip().split(" ")
                if
                "#" not in i and
                len(i.strip())
            ]
        )

        self.label.append(
            self.filename.split("/")[-1][-41:]
        )

        self.data.append(
            [
                [] for i in all_data[0].strip().split(" ") if len(i.strip())
            ]
        )

        for cnt, line in enumerate(all_data):
            for apnd, num in zip(
                self.data[-1], [
                    float(i) for i in line.strip().split(" ") if len(i.strip())
                    ]
            ):
                apnd.append(num)

        return

    def plot_latest_ResultFile(self):

        self.filename_result = self.get_latest_file(self.path, self.model_result)




    @staticmethod
    def get_xy_current_file(self):

        print("\n choose columns to plot \n")
        for cnt, value in enumerate(self.headline[-1]):
            print("\t {}: {}".format(cnt, value))

        print("\n x.. \n")
        index_x = self.get_index(self, len(self.headline[-1]))

        print("\n y.. \n")
        index_y = self.get_index(self, len(self.headline[-1]))

        return index_x, index_y

    @staticmethod
    def set_parms(self, ax, label_x, label_y):

        from matplotlib.ticker import FormatStrFormatter

        fontsize = 12
        ticksize = 10

        ax.clear()

        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_xlabel(label_x, fontsize=fontsize)
        ax.xaxis.set_tick_params(labelsize=ticksize)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_ylabel(label_y, fontsize=fontsize)
        ax.yaxis.set_tick_params(labelsize=ticksize)

        return

    @staticmethod
    def single_grid_placement(self, fig, gs):

        all_ax = []
        all_ax.append(fig.add_subplot(gs[:]))

        return all_ax

    @staticmethod
    def get_figure(self, nrows, ncols, grid_placement):

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=nrows, ncols=ncols)

        fig = plt.figure()
        fig.set_tight_layout(True)

        return fig, grid_placement(self, fig, gs)

    def plot_latest_file(self):

        self.clear()
        self.set_default_paths()
        self.model_fname = "STT_phiScal_J_"

        self.load_file_latest()

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        index_x, index_y = self.get_xy_current_file(self)

        fig, all_ax = self.get_figure(self, 1, 1, self.single_grid_placement)

        ax = all_ax[-1]

        self.set_parms(
            self, ax, self.headline[-1][index_x], self.headline[-1][index_y]
        )

        ax.plot(
            self.data[-1][index_x], self.data[-1][index_y],
            linewidth=1.5,
            label=self.label[-1]
        )

        ax.legend(loc="best", fontsize=8)
        plt.show()

        return

    @staticmethod
    def set_filenames(self, _path):

        import os

        list_dir = os.listdir(_path)

        for ind, val in enumerate(list_dir):
            print("{}: {}".format(ind, val))

        return list_dir[self.get_index(self, len(list_dir))]

    @staticmethod
    def units_coef_clac(self):
        # mas of sun in kg
        const_msun = 1.9891e30
        # gravitational const in m^3kg^-1s^-2
        const_g = 6.67384e-11
        # speed of light in ms^-1
        const_c = 299792458

        units = {}

        # units of density in g cm^-3
        units["density"] = 1e-3 * const_c**6 / (const_g**3 * const_msun**2)

        # units of pressure in dyns cm^-3
        units["pressure"] = const_c**8 / (const_g**3 * const_msun**2) * 10

        # units of rad coordinate in km
        units["rad"] = 1e-3 * const_g * const_msun / const_c**2

        # units of moment of inertia
        units["j"] = 1e7 * const_g**2 * const_msun**3 / const_c**4

        return units

    @staticmethod
    def load_file_daniela(self):

        with open(self.path_daniela + "/" + self.filename_daniela, "r") as f:
            all_data = f.readlines()

        # get rid of inline text
        for i, line in enumerate(all_data):
            if "rho" in line:
                del all_data[i]

        # i am not interested in all data, so only the indexes here
        # will be saved
        # since in some files the central pressure is in different column
        # we check it and take into account
        if self.filename_daniela in [
            "models_APR_beta-4.5_mphi5e-3.dat",
            "models_APR_beta-6.0_mphi0.dat",
            "models_APR_beta-6.0_mphi1e-3.dat"
        ]:
            indexes = [ 0, 1, 2, 3, 4, 10 ]
        else:
            indexes = [ 0, 1, 2, 3, 4, 11 ]

        _units = [
            self.units["density"], self.units["rad"],
            1, self.units["j"], 1, 1
        ]

        self.data_daniela.append( [ [] for i in indexes ] )

        for cnt, line in enumerate(all_data):
            for i, u, d in zip(indexes, _units, self.data_daniela[-1]):
                try:
                    d.append(float(line.split(" ")[i]) / u)
                except ValueError:
                    print(
                        "\n ValueError: line {}: {} \n".format(
                            cnt, line
                        )
                    )
                    break

        self.data_daniela[-1][-2] = \
            [ (-1)*i for i in self.data_daniela[-1][-2] ]

    def plot_compare_OneSet_mine_daniela(self):

        self.clear()

        self.set_default_paths()

        if self.check_for_path(self.path) and self.check_for_path(self.path_daniela):
            self.filename = self.set_filenames(self, self.path)
            self.filename_daniela = self.set_filenames(self, self.path_daniela)
        else:
            print("\n paths missing {} \t {},pass \n".format(
                self.path, self.path_daniela
                )
            )
            pass

        self.load_file(self)
        self.load_file_daniela(self)

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        index_x, index_y = self.get_xy_current_file(self)

        fig, all_ax = self.get_figure(self, 1, 1, self.single_grid_placement)

        ax = all_ax[-1]

        self.set_parms(
            self, ax, self.headline[-1][index_x], self.headline[-1][index_y]
        )

        daniela_mapping = {
            "rho_c": 0,
            "AR": 1,
            "M": 2,
            "J": 3,
            "phiScal_c": 4,
            "p_c": 5
        }

        ax.plot(
            self.data[-1][index_x], self.data[-1][index_y],
            linewidth=1.5,
            label=self.label[-1]
        )

        ax.plot(
            self.data_daniela[-1][daniela_mapping[self.headline[-1][index_x]]],
            self.data_daniela[-1][daniela_mapping[self.headline[-1][index_y]]],
            marker = "o",
            markersize = 5,
            alpha = 0.4,
            linewidth=0,
            label=self.filename_daniela
        )

        ax.legend(loc="best", fontsize=8)
        plt.show()

        return

class PlotResults_phiScal:

    def __init__(self):

        self.path = "/home/dimitar/projects/STT_theories/results"
        self.name = "STT_phiScal_"

        self.mapping = {
            0: "p_c",
            1: "phiScal_c",
            2: "M",
            3: "AR",
            4: "rho_c",
            5: "delta_phiScal"
        }

        self.data = []

        self.filename = self.get_latest(self)

        self.load_full_file(self)

        return

    @staticmethod
    def get_latest(self):

        import glob
        import os

        filename = \
            max(
                glob.glob(self.path + "/" + self.name + "*"),
                key=os.path.getctime
            )

        print("\n latest file: \n\t {} \n".format(filename))

        return filename

    def clear(self):

        self.data.clear()
        self.filename = ""

        return

    def reload(self):

        self.clear()
        self.filename = self.get_latest(self)

        self.load_full_file(self)

        return

    @staticmethod
    def load_full_file(self):

        import itertools

        with open(self.filename, "r") as f:
            all_data = f.read().split("#")

        headline = []
        for chunk in [m for m in all_data if len(m)]:

            headline.append(chunk.split("\n")[0])

        self.data = [[] for i in range(len(headline))]

        for my_chunk, block_data, cnt in zip(self.data, [m for m in all_data if len(m)], itertools.count(0)):

            block_data = block_data.split("\n")[1:]

            self.data[cnt] = [[] for k in block_data[0].split(" ") if len(k)]

            for line in block_data:
                for i, j in zip(self.data[cnt], [float(k) for k in line.split(" ") if len(k)]):

                    i.append(j)

        return

    @staticmethod
    def get_figure(self):

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=1, ncols=1)

        fig = plt.figure()
        fig.set_tight_layout(True)

        return fig, fig.add_subplot(gs[:])

    def help_mapping(self):

        for i in self.mapping.keys():
            print("{} -> {}".format(i, self.mapping[i]))

        return

    @staticmethod
    def set_param(self, ax, col_x, col_y):

        from matplotlib.ticker import FormatStrFormatter

        fontsize = 12
        ticksize = 10

        ax.clear()

        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_xlabel(self.mapping[col_x], fontsize=fontsize)
        ax.xaxis.set_tick_params(labelsize=ticksize)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_ylabel(self.mapping[col_y], fontsize=fontsize)
        ax.yaxis.set_tick_params(labelsize=ticksize)

        return

    def plot_it(self, col_x, col_y):

        import itertools

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        fig, ax = self.get_figure(self)

        self.set_param(self, ax, col_x, col_y)

        for cnt, chunk in enumerate(self.data):

            ax.plot(
                chunk[col_x], chunk[col_y],
                linewidth=1.5
            )

        plt.show()

        return


class PlotResults_phiScal_J:

    def __init__(self):

        self.path = "/home/dimitar/projects/STT_theories/results"
        self.name = "STT_phiScal_J_"

        self.mapping = {
            0: "p_c",
            1: "phiScal_c",
            2: "M",
            3: "AR",
            4: "rho_c",
            5: "J",
            6: "delta_phiScal",
            7: "delta_PhiMetr",
            8: "delta_Omega"
        }

        self.data = []

        self.filename = self.get_latest(self)

        self.load_full_file(self)

        return

    @staticmethod
    def get_latest(self):

        import glob
        import os

        filename = \
            max(
                glob.glob(self.path + "/" + self.name + "*"),
                key=os.path.getctime
            )

        print("\n latest file: \n\t {} \n".format(filename))

        return filename

    def clear(self):

        self.data.clear()
        self.filename = ""

        return

    def reload(self):

        self.clear()
        self.filename = self.get_latest(self)

        self.load_full_file(self)

        return

    @staticmethod
    def load_full_file(self):

        import itertools

        with open(self.filename, "r") as f:
            all_data = f.readlines()

        headline = [i for i in all_data.pop(
            0).strip().split(" ") if len(i) and "#" not in i]

        if len(headline) != len(self.mapping.keys()):
            print(
                "\n len headline diff from len self.mapping.keys() \n",
                headline,
                list(self.mapping.keys()),
                "\n PASSING \n"
            )
            pass

        self.data = [[]
                     for i in all_data[0].strip().split(" ") if len(i.strip())]

        for line in [i.strip() for i in all_data if len(i.strip())]:
            for apnd, num in zip(self.data, [float(i.strip()) for i in line.split(" ")]):
                apnd.append(num)

        return

    @staticmethod
    def get_figure(self):

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=1, ncols=1)

        fig = plt.figure()
        fig.set_tight_layout(True)

        return fig, fig.add_subplot(gs[:])

    def help_mapping(self):

        for i in self.mapping.keys():
            print("{} -> {}".format(i, self.mapping[i]))

        return

    @staticmethod
    def set_param(self, ax, col_x, col_y):

        from matplotlib.ticker import FormatStrFormatter

        fontsize = 12
        ticksize = 10

        ax.clear()

        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_xlabel(self.mapping[col_x], fontsize=fontsize)
        ax.xaxis.set_tick_params(labelsize=ticksize)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_ylabel(self.mapping[col_y], fontsize=fontsize)
        ax.yaxis.set_tick_params(labelsize=ticksize)

        return

    def plot_single_myRes(self, col_x, col_y):

        import itertools

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        fig, ax = self.get_figure(self)

        self.set_param(self, ax, col_x, col_y)

        for cnt, chunk in enumerate(self.data):

            ax.plot(
                chunk[col_x], chunk[col_y],
                linewidth=1.5
            )

        plt.show()

        return


class PlotResults_phiScal_J_compare:

    def __init__(self):

        self.data_mine = []
        self.headline_mine = []

        self.data_daniela = []
        self.headline_daniela = []

        self.mapping = {
            0: "p_c",
            1: "varPhi_c",
            4: "M",
            5: "AR",
            6: "rho_c",
            7: "J"
        }

        self.units = self.units_coef_clac(self)

    def help_mapping(self):

        for i in self.mapping.keys():
            print("{} -> {}".format(i, self.mapping[i]))

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
        units["density"] = 1e-3 * const_c**6 / (const_g**3 * const_msun**2)

        # units of pressure in dyns cm^-3
        units["pressure"] = const_c**8 / (const_g**3 * const_msun**2) * 10

        # units of rad coordinate in km
        units["rad"] = 1e-3 * const_g * const_msun / const_c**2

        # units of moment of inertia
        units["j"] = 1e7 * const_g**2 * const_msun**3 / const_c**4

        return units

    def load_ResFile_daniela(self, filename):

        print("\n GO TO THE FOLDER or will not now where is the central pressure \n")

        with open(filename, "r") as f:
            all_data = f.readlines()

        # get rid of inline text lines
        for i, line in enumerate(all_data):
            if "rho" in line:
                all_data.pop(i)

        all_data = [i.strip() for i in all_data if len]
        headline = [
            "rho", "AR", "M", "J", "phiScal_c", "p_c"
        ]

        # i am not interested in all data, so only the indexes of the corresponding
        # labels listed above will come here
        # in some the central pressure is in different column
        if filename in [
            "models_APR_beta-4.5_mphi5e-3.dat",
            "models_APR_beta-6.0_mphi0.dat",
            "models_APR_beta-6.0_mphi1e-3.dat"
        ]:
            indexes = [0, 1, 2, 3, 4, 10]
        else:
            indexes = [0, 1, 2, 3, 4, 11]

        if len(headline) != len(indexes):
            print(
                "\n headline and indexes len is different... \n",
                headline,
                indexes,
                "\n PASSING \n"
            )
            pass

        _units = [
            self.units["density"], self.units["rad"],
            1, self.units["j"], 1, 1
        ]

        if len(_units) != len(indexes):
            print(
                "\n amount of units and indexes is different \n",
                _units,
                indexes,
                "\n PASSING \n"
            )
            pass

        data = [[] for i in headline]

        for cnt, line in enumerate(all_data):
            for i, u, d in zip(indexes, _units, data):
                try:
                    d.append(float(line.split(" ")[i]) / u)
                except ValueError:
                    print("\n ValueError: line {}: {}\n".format(cnt, line))
                    break

        self.data_daniela.append(data)
        self.headline_daniela.append(headline)


if __name__ == "__main__":

    import sys

    print("\n Hello from {} \n", sys.argv[0])
