#!/usr/bin/env python3

import argparse
import configparser
import os
import itertools
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def units_coef_calc():
        """
        Calculates the unit coefficients

        Parameters
        ----------

        Returns
        -------
        : dictionary
            "density" in (double) g cm^-3
            "pressure" in (double) dyns cm^-3
            "r" in (double) km
            "j" in (double) m^2 kg
        """

        # mas of sun in kg
        const_msun = 1.9891e30
        # gravitational const in m^3kg^-1s^-2
        const_g = 6.67384e-11
        # speed of light in ms^-1
        const_c = 299792458

        global units
        units = {}

        # units of density in g cm^-3
        units["density"] = 1e-3 * const_c**6 / (const_g**3 * const_msun**2)

        # units of pressure in dyns cm^-3
        units["pressure"] = const_c**8 / (const_g**3 * const_msun**2) * 10

        # units of rad coordinate in km
        units["R"] = 1e-3 * const_g * const_msun / const_c**2

        # units of moment of inertia
        units["J"] = 1e7 * const_g**2 * const_msun**3 / const_c**4

        return

def map_ms_ls_c():

    global map_ms, map_ls, map_c

    #~ all marker styles which are nice to know
    all_marker_styles = [
        "s", "8", ">", "<", "^", "v", "o",
        "X", "P", "d", "D", "H", "h", "*", "p",
        "$\\bowtie$", "$\\clubsuit$", "$\\diamondsuit$", "$\\heartsuit$",
        "$\\spadesuit$",
        "$\\bigotimes$", "$\\bigoplus$",
    ]

    #~ all line styles which are nice to know
    #~ plus how to create a line style on your own
    all_line_styles = [
        "-.", "--",
        #~ (0, (5, 1, 1, 1, 1, 1)),
        (0, (5, 1, 1, 1, 1, 1, 1, 1)),
        #~ (0, (5, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
        (0, (8, 1, 1, 1, 3, 1, 1, 1))
    ]

    #~ all colors which are nice to know
    #~ plus some strange from matplotlib xkcdb
    all_colors = [
        "#e6194B", "#3cb44b", "#4363d8",
        "b", "g", "r", "c", "m", "y"
    ]

    #~ we are expecting to plot only those EOSs
    EOS_max18 = [
        "SLy", "APR2", "APR3", "APR4", "FPS", "WFF1", "WFF2", "WFF3", "BBB2",
        "ENG", "MPA1", "MS1", "MS2",  "MS1b", "GNH3", "H3", "H4", "ALF2", "ALF4"
    ]

    #~ map each EOS to different marker style
    map_ms = {
        "SLy": "s", "APR2": "8", "APR3": "$\\bigotimes$", "APR4": "<", "FPS":
        "^", "WFF1": "v", "WFF2": "o", "WFF3": "X", "BBB2": "P",
        "ENG":"$\\bigoplus$", "MPA1": "D", "MS1": "H", "MS2": "h",
        "MS1b": "*", "GNH3": "$\\spadesuit$", "H3": "$\\bowtie$",
        "H4": "$\\clubsuit$", "ALF2": "$\\diamondsuit$", "ALF4": "$\\heartsuit$"
    }

    #~ map each lambda value to different linestyle
    map_ls = {
        0: "-.", 1e-1: "--", 1e0: (0, (5, 1, 1, 1, 1, 1, 1, 1)),
        1e1: (0, (8, 1, 1, 1, 3, 1, 1, 1)),
        "GR": "-"
    }

    #~ map each m value to different color
    map_c = {
        0: "#e6194B", 5e-3: "#3cb44b", 5e-2: "#4363d8", "GR": "#f58231"
    }

    return

def load_YvsX_config(config):
    """
    the input is a list of tuples, each of them having as first element as
    entry key and the second a string, which should be converted to a list

    basically this takes the confg as a list of tuples and converts it to
    dictionary
    """

    #~ keys which values should be converted to numbers using eval
    eval_keys = [ "beta", "m", "lambda", "x_col", "y_col" ]

    #~ keys which values should be lists
    list_keys = [ "eos_name", "beta", "m", "lambda" ]

    config_dict = {}

    for entry in config:

        if entry[0] in eval_keys and entry[0] in list_keys:
            config_dict.update(
                { entry[0]: [
                    eval(_) for _ in entry[1].split(",") if _.strip()
                ] }
            )
        elif entry[0] not in eval_keys and entry[0] in list_keys:
            config_dict.update(
                { entry[0]: [
                    _.strip() for _ in entry[1].split(",") if _.strip()
                ] }
            )
        elif entry[0] in eval_keys:
            config_dict.update(
                { entry[0]: eval(entry[1]) }
            )
        else:
            config_dict.update(
                { entry[0]: entry[1].strip() }
            )

    return config_dict

def get_YvsX_ax():
    """
    set plot variables and return axes to created figure
    """

    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)

    plt.rc('font', size=15)

    plt.rc('figure', autolayout=True)
    plt.rc('figure', figsize=(7.2,4.45))


    plt.rc('axes', titlesize=16)
    plt.rc('axes', labelsize=17)

    plt.rc('lines', linewidth=2)
    plt.rc('lines', markersize=10)

    plt.rc('legend', fontsize=13)

    plt.rc('text', usetex=True)
    plt.rc('mathtext', fontset="stix")

    #~ plt.rc('font', family='serif')
    plt.rc('font', family='STIXGeneral')

    fig, ax = plt.subplots()

    fig.set_rasterized(True)

    return ax

def set_axYvsX_parms(ax, x_label = "", y_label = ""):

    #~ fontsize=16
    #~ x_ticksize = 16
    #~ y_ticksize = 16

    #~ direction of the tickks
    ax.tick_params(direction="in")

    #~ set the labels and their size
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    #~ how big the tick labels should be
    #~ ax.xaxis.set_tick_params(labelsize=x_ticksize)
    #~ ax.yaxis.set_tick_params(labelsize=y_ticksize)

    #~ if diffrent format should be applied to the tick numbers
    #~ ax.xaxis.set_major_formatter(FormatStrFormatter(x_format_str))
    #~ ax.yaxis.set_major_formatter(FormatStrFormatter(y_format_str))

    return

def get_YvsX_data( fpath, fname, EOSname, EOSbeta, EOSm, EOSlambda, col_x, col_y):
    """
    retrieve the data of the model
    """

    EOSmodel = "_".join( [
        fname,
        EOSname,
        "beta{:.3e}".format(EOSbeta),
        "m{:.3e}".format(EOSm),
        "lambda{:.3e}".format(EOSlambda)
    ] )

    fullPathModel = os.path.join( fpath, EOSname, EOSmodel)

    print("\n Will load data for: \n\t {} \n".format(fullPathModel))

    # TODO if file not exists what do
    data = np.loadtxt(fullPathModel, comments="#", delimiter=" ",
    usecols=(col_x, col_y))

    return data

def convert_NumScientific(num):

    snum = "{:.2e}".format(num)
    mantisa, power = snum.split("e")

    mantisa = "{}".format(mantisa)
    power = "{}".format(power)

    return " ${{{:.0f}}} \\times 10^{{{:.0f}}}$".format(
        float(mantisa), float(power)
    ) if float(mantisa) \
    else " ${{{:.0f}}}$".format(
        float(mantisa)
    )

def plot_MvsR_GR(config):

    config = load_YvsX_config(config)

    ax = get_YvsX_ax()
    set_axYvsX_parms(ax, x_label = config["x_label"], y_label = config["y_label"])

    for model in itertools.product( config["eos_name"], config["beta"],
    config["m"], config["lambda"] ):

        data = get_YvsX_data(
            config["path"], config["base_name"], model[0],
            model[1], model[2], model[3], config["x_col"], config["y_col"]
        )

        #~ remove entries who have mass less than min_mass
        min_mass = 0.5
        data = data[~(data[:,1]<min_mass), :]

        #~ include only "stable" masses, those before the max mass
        data = np.delete(data, np.s_[np.argmax(data[:,1]):], axis=0)

        ax.plot(
            data[:,0]*units.get(config["x_unit"], 1),
            data[:,1]*units.get(config["y_unit"], 1),
            label = None,
            markevery = 0.1,
            marker = map_ms.get(model[0], None),
            color = map_c.get(model[2], None),
            linestyle = map_ls.get(model[3], None)
        )

        #~ NOW SAME PROCEDURE BUT FOR GR CASE
        data = get_YvsX_data(
            config["path"], config["base_name"], model[0],
            0, 0, 0, config["x_col"], config["y_col"]
        )

        data = data[~(data[:,1]<min_mass), :]
        data = np.delete(data, np.s_[np.argmax(data[:,1]):], axis=0)

        ax.plot(
            data[:,0]*units.get(config["x_unit"], 1),
            data[:,1]*units.get(config["y_unit"], 1),
            label = None,
            markevery = 0.1,
            marker = map_ms.get(model[0], None),
            color = map_c.get("GR", None),
            linestyle = map_ls.get("GR", None)
        )

    # lets add a legend for marker styles
    ax.add_artist( ax.legend(
        handles = [
            Line2D(
                [], [], color = "k",
                marker = map_ms.get(_, None),
                linewidth = 0, linestyle = None,
                label = _
            ) for _ in config["eos_name"]
        ],
        loc = "upper left",
        ncol = 4
    ) )

    # lets add a legend for line style
    ax.add_artist( ax.legend(
        handles = [
            Line2D(
                [], [], color = "k",
                marker = None,
                linestyle = map_ls.get(_, None),
                label = (
                    "$\lambda = {}$".format(_) if _ != "GR" else "GR"
                )
            ) for _ in [ *config["lambda"], "GR" ]
        ],
        loc = "lower right",
        ncol = 1,
        handlelength = 3
    ) )

    # lets add a legend for color
    ax.add_artist( ax.legend(
        handles = [
            mpatches.Patch(
                color = map_c.get(_, None),
                label = (
                    "$m = $ {}".format(convert_NumScientific(_)) if _ != "GR" else "GR"
                )
            ) for _ in [ *config["m"], "GR" ]
        ],
        loc = "center right",
        ncol = 1,
        #~ handlelength = 3
    ) )

    plt.savefig(
        'MvsR_STT_GR.eps', format="eps",
        #~ bbox_inches='tight',
        dpi=200,
        pad_inches=0
    )

    plt.show()

    return

def create_uniI_data(config):

    def _load_config(config):
        """
        the input is a list of tuples, each of them having as first element as
        entry key and the second a string, which should be converted to a list

        basically this takes the confg as a list of tuples and converts it to
        dictionary
        """

        #~ keys which values should be converted to numbers using eval
        eval_keys = [ "beta", "m", "lambda", "col_r", "col_m", "col_j" ]

        #~ keys which values should be lists
        list_keys = [ "eos_name", "beta", "m", "lambda" ]

        config_dict = {}

        for entry in config:

            if entry[0] in eval_keys and entry[0] in list_keys:
                config_dict.update(
                    { entry[0]: [
                        eval(_) for _ in entry[1].split(",") if _.strip()
                    ] }
                )
            elif entry[0] not in eval_keys and entry[0] in list_keys:
                config_dict.update(
                    { entry[0]: [
                        _.strip() for _ in entry[1].split(",") if _.strip()
                    ] }
                )
            elif entry[0] in eval_keys:
                config_dict.update(
                    { entry[0]: eval(entry[1]) }
                )
            else:
                config_dict.update(
                    { entry[0]: entry[1].strip() }
                )

        return config_dict

    def _get_data( fpath, fname, EOSname, EOSbeta, EOSm, EOSlambda, col_r,
    col_m, col_j ):
        """
        retrieve the data of the model
        """

        EOSmodel = "_".join( [
            fname,
            EOSname,
            "beta{:.3e}".format(EOSbeta),
            "m{:.3e}".format(EOSm),
            "lambda{:.3e}".format(EOSlambda)
        ] )

        fullPathModel = os.path.join( fpath, EOSname, EOSmodel)

        print("\n Will load data for: \n\t {} \n".format(fullPathModel))

        # TODO if file not exists what do
        data = np.loadtxt(fullPathModel, comments="#", delimiter=" ",
        usecols=(col_r, col_m, col_j))

        return data

    config = _load_config(config)

    for model in itertools.product( config["eos_name"], config["beta"],
    config["m"], config["lambda"] ):

        data = _get_data(
            config["path"], config["base_name"], model[0],
            model[1], model[2], model[3],
            config["col_r"], config["col_m"], config["col_j"]
        )

        #remove entries who have mass less than min_mass
        min_mass = 0.5
        data = data[~(data[:,1]<min_mass), :]

        #include only "stable" masses, those before the max mass
        data = np.delete(data, np.s_[np.argmax(data[:,1]):], axis=0)

        EOSmodel = "_".join( [
            config["base_name"],
            model[0],
            "beta{:.3e}".format(model[1]),
            "m{:.3e}".format(model[2]),
            "lambda{:.3e}".format(model[3]),
        ] )

        fullModelPath = os.path.join(
            config["path"], model[0], "Fitting"
        )

        pathlib.Path( fullModelPath ).mkdir(parents=True, exist_ok=True)

        target = os.path.join(fullModelPath, "_".join([EOSmodel, "tildeI"]))
        print("\n\t will convert universal I \n\t\t {} \n".format( target ) )
        np.savetxt(
            target,
            np.column_stack( (
                data[:,1]/data[:,0], data[:,2]/(data[:,1]*data[:,0]**2)
            ) ),
            delimiter = " ",
            newline = "\n",
            header = "M/R I/(MR**2)",
            comments="# "
        )

        target = os.path.join(fullModelPath, "_".join([EOSmodel, "barI2"]))
        print("\n\t will convert universal I \n\t\t {} \n".format( target ) )
        np.savetxt(
            target,
            np.column_stack( (
                data[:,1]/data[:,0], data[:,2]/(data[:,1]**3)
            ) ),
            delimiter = " ",
            header = "M/R I/(M**3)",
            comments="# "
        )

    return

if __name__ == "__main__":

    def _parser_init():

        parser = argparse.ArgumentParser(
            prog = "My personal result plotter",
            description="""
            Quick plotting tool to suit my needs. Basically it has limited regimes of
            work, each corresponding to different keyword and thus different
            configuration settings and output.
            """
        )

        #argument for the config file
        parser.add_argument(
            "--config",
            action = "store",
            nargs = "?",
            type = str,
            default = "result_plotter.config",
            const = "result_plotter.config",
            metavar = "path to config file",
            dest = "ConfigFile",
            help = "path to configuration file for each option, presented as args"
        )

        #argument for MvsR with GR included
        parser.add_argument(
            "--MvsR_GR",
            action = "store_const",
            #~ nargs = "?",
            #~ type = str,
            #~ default = "quick_plotter.config",
            const = "MvsR_GR",
            #~ metavar = "",
            dest = "ConfigSection",
            required = False,
            help = """
                if called, it will use MvsR_GR section in config file to produce
                a plot
            """
        )

        #argument for creating all uni I in separated directories
        parser.add_argument(
            "--create_uniI_data",
            action = "store_const",
            #~ nargs = "?",
            #~ type = str,
            #~ default = "quick_plotter.config",
            const = "create_uniI_data",
            #~ metavar = "",
            dest = "ConfigSection",
            required = False,
            help = """
                if called, it will use just create all uni I files needed to
                create the universality plots and read path to results from
                create_uniI_data section of the config file
            """
        )

        #argument for uni I
        parser.add_argument(
            "--tildeI_GR",
            action = "store_const",
            #~ nargs = "?",
            #~ type = str,
            #~ default = "quick_plotter.config",
            const = "tildeI_GR",
            #~ metavar = "",
            dest = "ConfigSection",
            required = False,
            help = """
                if called, it will use tildeI_GR section in config file to produce
                a plot
            """
        )

        return parser

    def _get_config(ConfigFile):

        # TODO no input control !!!

        config = configparser.ConfigParser()
        config.read(ConfigFile)

        return config

    units_coef_calc()
    map_ms_ls_c()

    args = _parser_init().parse_args()

    config = _get_config(args.ConfigFile)

    if args.ConfigSection == "MvsR_GR":

        plot_MvsR_GR(config.items(args.ConfigSection))

    elif args.ConfigSection == "create_uniI_data":

        create_uniI_data(config.items(args.ConfigSection))

    else:
        print("\n {} unknown, terminating... \n", args.ConfigSection)
