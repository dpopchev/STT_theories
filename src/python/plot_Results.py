#!/usr/bin/env python

class PlotResults:

    def __init__(self):

        self.path = "/home/dimitar/projects/STT_theories/results"
        self.name = "STT_"

        self.mapping = {
            0 : "p_c",
            1 : "phiScal_c",
            2 : "M",
            3 : "AR",
            4 : "rho_c",
            5 : "delta_phiScal"
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
            key = os.path.getctime
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

        with open( self.filename, "r" ) as f:
            all_data = f.read().split("#")

        headline = []
        for chunk in [ m for m in all_data if len(m) ]:

            headline.append(chunk.split("\n")[0])

        self.data = [ [] for i in range(len(headline)) ]

        for my_chunk, block_data, cnt in zip( self.data, [ m for m in all_data if len(m) ], itertools.count(0) ):

            block_data = block_data.split("\n")[1:]

            self.data[cnt] = [ [] for k in block_data[0].split(" ") if len(k) ]

            for line in block_data:
                for i,j in zip(self.data[cnt], [ float(k) for k in line.split(" ") if len(k) ]):

                    i.append(j)

        return

    @staticmethod
    def get_figure(self):

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=1,ncols=1)

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
        ax.set_xlabel(self.mapping[col_x],fontsize=fontsize)
        ax.xaxis.set_tick_params(labelsize=ticksize)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_ylabel(self.mapping[col_y],fontsize=fontsize)
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
              linewidth = 1.5
            )

        plt.show()

        return

if __name__ == "__main__":

    print("\n Hello world \n")
