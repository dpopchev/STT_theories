#!/usr/bin/env python

import argparse
import configparser
import os
import itertools
import pathlib
import random

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

from matplotlib.lines import Line2D
from matplotlib import gridspec
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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
        (0, (5, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
        #~ (0, (8, 1, 1, 1, 3, 1, 1, 1))
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
        0: "--",
        1e-1: "-.",
        1e0: (0, (6, 1, 1, 1, 1, 1)),
        1e1: (0, (6, 1, 1, 1, 1, 1, 1, 1)),
        "GR": "-"
    }

    #~ map each m value to different color
    map_c = {
        0: "#e6194B", 5e-3: "#3cb44b", 5e-2: "#4363d8", "GR": "#f58231"
    }

    return

def luminosity_color(color, amount=1.2):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color

    c = colorsys.rgb_to_hls(*mc.to_rgb(c))

    return colorsys.hls_to_rgb(abs(c[0]), abs(1 - amount * (1 - c[1])), abs(c[2]))

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

    plt.rc('figure', figsize=(7.2,4.8))

    plt.rc('axes', titlesize=16)
    plt.rc('axes', labelsize=16)

    plt.rc('lines', linewidth=2)
    plt.rc('lines', markersize=10)

    plt.rc('legend', fontsize=8)

    plt.rc('text', usetex=True)
    plt.rc('mathtext', fontset="stix")

    #~ plt.rc('font', family='serif')
    plt.rc('font', family='STIXGeneral')

    fig, ax = plt.subplots()

    fig.set_rasterized(True)
    fig.set_tight_layout(False)

    return ax

def get_uniI_ax():
    """
    set plot variables and return axes to created figure
    """

    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)

    plt.rc('font', size=15)

    plt.rc('figure', figsize=(7.2,4.8))

    plt.rc('axes', titlesize=16)
    plt.rc('axes', labelsize=16)

    plt.rc('lines', linewidth=2)
    plt.rc('lines', markersize=10)

    plt.rc('legend', fontsize=8)

    plt.rc('text', usetex=True)
    plt.rc('mathtext', fontset="stix")

    #~ plt.rc('font', family='serif')
    plt.rc('font', family='STIXGeneral')

    fig, (ax_up, ax_down) = plt.subplots(
        2,1, sharex=True,
        gridspec_kw=dict(height_ratios=[2, 1])
    )

    fig.subplots_adjust(hspace=0)
    fig.set_rasterized(True)
    fig.set_tight_layout(False)

    return ax_up, ax_down

def get_uniI_data( fpath_main, fpath_fitting, fname, EOSname, EOSbeta, EOSm,
EOSlambda, col_x, col_y, tilde ):
    """
    retrieve the data of the model
    """

    EOSmodel = "_".join( [
        fname,
        EOSname,
        "beta{:.3e}".format(EOSbeta),
        "m{:.3e}".format(EOSm),
        "lambda{:.3e}".format(EOSlambda),
        tilde
    ] )

    fullPathModel = os.path.join( fpath_main, EOSname, fpath_fitting, EOSmodel)

    print("\n Will load data for: \n\t {} \n".format(fullPathModel))

    # TODO if file not exists what do
    data = np.loadtxt(fullPathModel, comments="#", delimiter=" ",
    usecols=(col_x, col_y))

    min_compactness = 0.09
    data = data[~(data[:,0]<min_compactness), :]

    return data

def plotMarkers_getFits_uniI(
    config, pars, tilde, ax_up, ax_down, ax_in = None, fit_results = {}
):

    def _get_polyfit_res_tildeI(xp, yp):

        coef, rest = polyfit(
            x = xp, y = yp,
            deg = [ 0, 1, 4 ],
            w = np.sqrt(yp),
            full = True
        )

        #~ calcualte the chi reduce
        chi_red = rest[0][0]/(len(xp) - 3)
        p = lambda x: coef[0] + coef[1]*x + coef[4]*x**4

        return coef, chi_red, p

    def _get_polyfit_res_barI(xp, yp):

        coef, rest = polyfit(
            x = xp, y = yp,
            deg = [ 1, 2, 3, 4 ],
            #~ w = np.sqrt(yp),
            full = True
        )

        #~ calcualte the chi reduce
        chi_red = rest[0][0]/(len(xp) - 4)
        p = lambda x: coef[1]*x**-1 + coef[2]*x**-2 + coef[3]*x**-3 + coef[4]*x**-4

        return coef, chi_red, p

    # TODO this choice is arbitrary, it should not be hardcoded
    # MS2 is soft; SLy is moderate; APR3 is stiff
    #~ PlaceMarkersForEOS = ["MS2", "SLy", "APR3"]

    plot_alpha = 0.5 # give some transperancy
    #~ fit_results = {} # to print what we have achieved at the end for the model

    # accumulate data for the model here, and use them to make the fit
    data_to_fit = np.empty(shape=(0,2))

    # fill the up plot with markers while gathering data for the fit
    for EOSname in config["eos_name"]:

        data = get_uniI_data(
            config["path_main"], config["path_fitting"],
            config["base_name"], EOSname, pars[0], pars[1], pars[2],
            config["x_col"], config["y_col"], tilde
        )

        data_to_fit = np.append(data_to_fit, data, axis=0)
        if pars[2] > 1:
            continue

        ax_up.plot(
            data[:,0], data[:,1],
            label = None,
            markevery = random.uniform(0.14,0.18),
            markersize = 6,
            linewidth = 0,
            marker = map_ms.get(EOSname, None),
            color = map_c.get(
                pars[1] if pars[0] else "GR",None
            ),
            linestyle = None,
            alpha = plot_alpha
        )

        # do not want any markers in the zoomed box
        if ax_in and False:
            ax_in.plot(
                data[:,0], data[:,1],
                label = None,
                markevery = None,
                markersize = 6,
                linewidth = 0,
                marker = map_ms.get(EOSname, None),
                color = map_c.get(pars[1],None),
                linestyle = None,
                alpha = plot_alpha
            )

    # get the coef in list, the Chi reduced score, and the polynom itself
    coef, chi_red, p = (
        _get_polyfit_res_tildeI(
            data_to_fit[:,0], data_to_fit[:,1]
        )  if tilde == "tildeI"
        else  _get_polyfit_res_barI(
            data_to_fit[:,0]**-1, data_to_fit[:,1]
        )
    )

    # average over all EOS of all the residuals
    delta_all = 0
    n_all = 0

    # average over all EOS of largest residual per EOS
    n_all_max = 0
    delta_all_max = 0

    # the largest residual over all EOS
    delta_max = 0

    # fill the down plot with residiums while gathering data for avreages
    for EOSname in config["eos_name"]:

        # change the markevery a little bit random

        data = get_uniI_data(
            config["path_main"], config["path_fitting"],
            config["base_name"], EOSname, pars[0], pars[1], pars[2],
            config["x_col"], config["y_col"], tilde
        )

        _data = np.abs( 1 - data[:,1]/p(data[:,0]) )

        delta_all += np.sum(_data)
        n_all += _data.size

        delta_all_max += np.amax(_data)
        n_all_max += 1

        delta_max = delta_max if delta_max > np.amax(_data) else np.amax(_data)

        if pars[2] >= 1:
            continue

        ax_down.plot(
            data[:,0],
            _data,
            label = None,
            linewidth = 0,
            markersize = 6,
            markevery = random.uniform(0.14,0.18),
            marker = map_ms.get(EOSname, None),
            color = map_c.get(pars[1] if pars[0] else "GR", None),
            linestyle = None,
            alpha = plot_alpha,
            zorder = 90 if pars[0] == 0 else 10
        )

    avg_L_1 = delta_all/n_all
    avg_L_inf = delta_all_max/n_all_max
    L_inf_worst = delta_max

    fit_results.update(
        {
            "beta {}, m = {}, lambda = {}".format(pars[0], pars[1], pars[2]):
            "a0 = {:.3e}, a1 = {:.3e}, a4 = {:.3e}, chi = {:.3e}, L1 = {:.3e}, Linf = {:.3e}, LinfW = {:.3e}".format(
                coef[0], coef[1], coef[4], chi_red, avg_L_1, avg_L_inf, L_inf_worst
        ) } if tilde == "tildeI"
        else {
            "beta {}, m = {}, lambda = {}".format(pars[0], pars[1], pars[2]):\
            "a1 = {:.3e}, a2 = {:.3e}, a3 = {:.3e}, a4 = {:.3e}, chi = {:.3e}, L1 = {:.3e}, Linf = {:.3e}, LinfW = {:.3e}".format(
            coef[1], coef[2], coef[3], coef[4], chi_red, avg_L_1, avg_L_inf, L_inf_worst
        ) }
    )

    #give min and max compactnesses by hand
    #~ p_x = np.linspace(np.amin(data_to_fit[:,0]), np.amax(data_to_fit[:,0]), 100)
    # TODO maybe not hardcode the compactneses like that? now 0.09 to 0.35
    p_x = np.linspace(0.09, 0.35, 100)

    ax_up.plot(
        p_x,
        p(p_x),
        label = None,
        color = luminosity_color(map_c.get(pars[1] if pars[0] else "GR", None)),
        linestyle = map_ls.get(pars[2] if pars[0] else "GR", None),
        zorder = 90 if pars[0] == 0 else 85,
        linewidth = 2.5
    )

    if ax_in:
        ax_in.plot(
            p_x,
            p(p_x),
            label = None,
            color = luminosity_color(map_c.get(pars[1] if pars[0] else "GR", None)),
            linestyle = map_ls.get(pars[2] if pars[0] else "GR", None),
            zorder = 90 if pars[0] == 0 else 85,
            linewidth = 2.5
        )

    # beta = 0 means GR, then make the worst and not so worst fit
    if pars[0] == 0:

        ax_up.fill_between(
            p_x,
            p(p_x)*(1 + avg_L_inf),
            p(p_x)*(1 - avg_L_inf),
            facecolor=map_c.get("GR", None),
            alpha= 0.4,
            zorder = 0
        )

        if ax_in:
            ax_in.fill_between(
                p_x,
                p(p_x)*(1 + avg_L_inf),
                p(p_x)*(1 - avg_L_inf),
                facecolor=map_c.get("GR", None),
                alpha= 0.4,
                zorder = 0
            )

        ax_up.fill_between(
            p_x,
            p(p_x)*( 1 + L_inf_worst ),
            p(p_x)*( 1 - L_inf_worst ),
            facecolor=map_c.get("GR", None),
            alpha= 0.3,
            zorder = 0
        )

        if ax_in:
            ax_in.fill_between(
                p_x,
                p(p_x)*( 1 + L_inf_worst ),
                p(p_x)*( 1 - L_inf_worst ),
                facecolor=map_c.get("GR", None),
                alpha= 0.3,
                zorder = 0
            )

    return fit_results

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

def plot_MvsR_GR(config, GRonly = False):

    config = load_YvsX_config(config)

    ax = get_YvsX_ax()
    set_axYvsX_parms(ax, x_label = config["x_label"], y_label = config["y_label"])

    for model in itertools.product( config["eos_name"], config["beta"],
    config["m"], config["lambda"] ):

        min_mass = 0.5

        # TODO GRonly variable should go to the config file
        if not GRonly:

            data = get_YvsX_data(
                config["path"], config["base_name"], model[0],
                model[1], model[2], model[3], config["x_col"], config["y_col"]
            )

            #~ remove entries who have mass less than min_mass
            data = data[~(data[:,1]<min_mass), :]

            #~ include only "stable" masses, those before the max mass
            #~ data = np.delete(data, np.s_[np.argmax(data[:,1]):], axis=0)

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
        #~ data = np.delete(data, np.s_[np.argmax(data[:,1]):], axis=0)

        ax.plot(
            data[:,0]*units.get(config["x_unit"], 1),
            data[:,1]*units.get(config["y_unit"], 1),
            label = None,
            markevery = 0.1,
            marker = map_ms.get(model[0], None),
            color = map_c.get("GR", None),
            linestyle = map_ls.get("GR", None),
            zorder = 100
        )

    handle_EOSnames = [
        Line2D(
            [], [],
            color = "k",
            marker = map_ms.get( _, None),
            linewidth = 0,
            linestyle = None,
            label = _
        ) for _ in config["eos_name"]
    ]

    handle_linestyles = [
        Line2D(
            [], [],
            color = "k",
            marker = None,
            linestyle = map_ls.get(_, None),
            label = (
                "$\lambda = {}$".format(_)
                if _ != "GR" else "GR" )
        ) for _ in [ *config["lambda"], "GR" ]
    ] if not GRonly else []

    handle_markers = [
        mpatches.Patch(
            color = map_c.get(_, None),
            label = (
                "$m= $ {}".format(convert_NumScientific(_))
                if _ != "GR" else "GR"
            )
        ) for _ in [ *config["m"], "GR" ]
    ] if not GRonly else []

    # lets add a legend for marker styles
    ax.add_artist( ax.legend(
        handles = [
            *handle_EOSnames, *handle_linestyles, *handle_markers
        ],
        loc = "lower right",
        ncol = 1,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
    ) )

    ax.set_xlim(7.75, 15.25)
    ax.set_ylim(0.25, 3)

    plt.savefig(
        'MvsR_GR.eps' if GRonly else "MvsR_STT_GR.eps",
        format="eps",
        dpi=600,
        pad_inches=0,
        bbox_inches='tight',
        papertype = "a4"
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

        target = os.path.join(fullModelPath, "_".join([EOSmodel, "barI"]))
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

def plot_tildeI_GR_zoomBox(config):

    config = load_YvsX_config(config)

    ax_up, ax_down = get_uniI_ax()

    ax_in = zoomed_inset_axes(ax_up, 3, loc=8)

    set_axYvsX_parms(ax_up, x_label = "", y_label = config["y_up_label"])
    set_axYvsX_parms(ax_down, x_label = config["x_label"], y_label = config["y_down_label"])
    ax_down.set_yscale("log")

    FitResults = {}

    for pars in itertools.product(config["beta"], config["m"], config["lambda"]):

        FitResults = plotMarkers_getFits_uniI(config, pars, "tildeI", ax_up, ax_down, ax_in)

    FitResults.update(plotMarkers_getFits_uniI(config, [0,0,0], "tildeI", ax_up, ax_down, ax_in))

    for k, v in FitResults.items():
        print("\n\t {} \n\t\t {}".format(k, v))

    handle_EOSnames = [
        Line2D(
            [], [],
            color = "k",
            marker = map_ms.get( _, None),
            linewidth = 0,
            linestyle = None,
            label = _
        ) for _ in config["eos_name"]
    ]

    handle_linestyles = [
        Line2D(
            [], [],
            color = "k",
            marker = None,
            linestyle = map_ls.get(_, None),
            label = (
                "$\lambda = {}$".format(_)
                if _ != "GR" else "GR" )
        ) for _ in [ *config["lambda"], "GR" ]
    ]

    # colors will be presented as patches
    handle_colors = [
        mpatches.Patch(
            color = map_c.get(_, None),
            label = (
                "$m= $ {}".format(convert_NumScientific(_))
                if _ != "GR" else "GR"
            )
        ) for _ in [ *config["m"], "GR" ]
    ]

    # colors will be presented as solid lines
    #~ handle_colors = [
        #~ Line2D(
            #~ [], [],
            #~ color = map_c.get(_, None),
            #~ marker = None,
            #~ label = (
                #~ "$\lambda = {}$".format(_)
                #~ if _ != "GR" else "GR" )
        #~ ) for _ in [ *config["lambda"], "GR" ]
    #~ ]

    legend_markers = ax_up.legend(
        handles = [ *handle_EOSnames ],
        loc = 2,
        ncol = 1,
        frameon = False,
        markerscale = 0.6,
        fancybox=True,
        #~ framealpha = 0.5,
        handlelength = 4,
        bbox_to_anchor=(1, 1),
        #~ borderaxespad=0.1,
    )

    legend_linestyles = ax_up.legend(
        handles = [ *handle_linestyles ],
        loc = 2,
        ncol = 1,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
        title = "Line patterns"
        #~ bbox_to_anchor=(0.79, 0.5),
        #~ borderaxespad=0.1
    )

    #~ depricated way of seting font size of tittle of legend
    #~ seem not the way any more in matplotlib 3
    #~ legend_linestyles.set_title('location', prop={"size":10})
    legend_linestyles.get_title().set_fontsize(10)

    legend_colors = ax_up.legend(
        handles = [ *handle_colors ],
        loc = 4,
        ncol = 1,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
        title = "Line colours"
        #~ bbox_to_anchor=(0.99, 0.22),
        #~ borderaxespad=0.1
    )

    legend_colors.get_title().set_fontsize(10)

    ax_up.add_artist(legend_markers)

    ax_up.add_artist( legend_linestyles )

    ax_up.add_artist( legend_colors )

    ax_up.set_xlim(0.09, 0.325)
    ax_up.set_ylim(0.28, 0.55)

    ax_in.set_xlim(0.255, 0.275)
    ax_in.set_ylim(0.43, 0.45)
    ax_in.xaxis.set_visible(False)
    ax_in.yaxis.set_visible(False)
    #~ ax_in.xaxis.set_tick_params(labelsize=8)
    #~ ax_in.yaxis.set_tick_params(labelsize=8)
    #~ ax_in.tick_params(axis="both",direction="in", pad=-10)

    mark_inset(
        ax_up, ax_in,
        loc1=2, loc2=4,
        fc="black", ec="black",
        zorder=110
    )

    ax_down.set_ylim(1e-3, 1.5e0)
    #~ ax_down.axhline(y=0.1, linewidth=2, color='r', alpha = 0.5)

    plt.savefig(
        'TildeI_all.eps', format="eps",
        dpi=600,
        pad_inches=0,
        bbox_inches='tight',
        papertype = "a4",
        bbox_extra_artists=(legend_markers,),
    )

    plt.show()

    return

def plot_tildeI_GR(config):

    config = load_YvsX_config(config)

    ax_up, ax_down = get_uniI_ax()

    set_axYvsX_parms(ax_up, x_label = "", y_label = config["y_up_label"])
    set_axYvsX_parms(ax_down, x_label = config["x_label"], y_label = config["y_down_label"])
    ax_down.set_yscale("log")

    FitResults = {}

    for pars in itertools.product(config["beta"], config["m"], config["lambda"]):

        FitResults = plotMarkers_getFits_uniI(config, pars, "tildeI", ax_up, ax_down)

    FitResults.update(plotMarkers_getFits_uniI(config, [0,0,0], "tildeI", ax_up, ax_down))

    for k, v in FitResults.items():
        print("\n\t {} \n\t\t {}".format(k, v))

    handle_EOSnames = [
        Line2D(
            [], [],
            color = "k",
            marker = map_ms.get( _, None),
            linewidth = 0,
            linestyle = None,
            label = _
        ) for _ in config["eos_name"]
    ]

    handle_linestyles = [
        Line2D(
            [], [],
            color = "k",
            marker = None,
            linestyle = map_ls.get(_, None),
            label = (
                "$\lambda = {}$".format(_)
                if _ != "GR" else "GR" )
        ) for _ in [ *config["lambda"], "GR" ]
    ]

    handle_colors = [
        mpatches.Patch(
            color = map_c.get(_, None),
            label = (
                "$m= $ {}".format(convert_NumScientific(_))
                if _ != "GR" else "GR"
            )
        ) for _ in [ *config["m"], "GR" ]
    ]

    # colors will be presented as solid lines
    #~ handle_colors = [
        #~ Line2D(
            #~ [], [],
            #~ color = map_c.get(_, None),
            #~ marker = None,
            #~ label = (
                #~ "$\lambda = {}$".format(_)
                #~ if _ != "GR" else "GR" )
        #~ ) for _ in [ *config["lambda"], "GR" ]
    #~ ]

    legend_markers = ax_up.legend(
        handles = [ *handle_EOSnames ],
        loc = 2,
        ncol = 2,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
        #~ bbox_to_anchor=(1, 1),
        #~ borderaxespad=0.1,
    )

    ax_up.add_artist(legend_markers)

    if len(handle_linestyles) > 2:
        legend_linestyles = ax_up.legend(
            handles = [ *handle_linestyles ],
            loc = 4,
            ncol = 1,
            frameon = True,
            markerscale = 0.6,
            fancybox=True,
            framealpha = 0.5,
            handlelength = 4,
            title = r"$m_\varphi = 5\times10^{-3}$"
            #~ bbox_to_anchor=(0.79, 0.5),
            #~ borderaxespad=0.1
        )

        legend_linestyles.get_title().set_fontsize(10)

        ax_up.add_artist( legend_linestyles )

    if len(handle_colors) > 2:
        legend_colors = ax_up.legend(
            handles = [ *handle_colors ],
            loc = 4,
            ncol = 1,
            frameon = True,
            markerscale = 0.6,
            fancybox=True,
            framealpha = 0.5,
            handlelength = 4,
            title = r"$\lambda = 0.1$"
            #~ bbox_to_anchor=(0.99, 0.22),
            #~ borderaxespad=0.1
        )

        ax_up.add_artist( legend_colors )

        legend_colors.get_title().set_fontsize(10)

    ax_up.set_xlim(0.09, 0.325)
    ax_up.set_ylim(0.28, 0.55)

    ax_down.set_ylim(1e-3, 1.5e0)
    #~ ax_down.axhline(y=0.1, linewidth=2, color='r', alpha = 0.5)

    plt.savefig(
        'TildeI_lambda1e-1.eps', format="eps",
        dpi=600,
        pad_inches=0,
        bbox_inches='tight',
        papertype = "a4",
        bbox_extra_artists=(legend_markers,)
    )

    plt.show()

    return

def plot_barI_GR(config):

    config = load_YvsX_config(config)

    ax_up, ax_down = get_uniI_ax()

    set_axYvsX_parms(ax_up, x_label = "", y_label = config["y_up_label"])
    set_axYvsX_parms(ax_down, x_label = config["x_label"], y_label = config["y_down_label"])
    ax_down.set_yscale("log")

    FitResults = {}

    for pars in itertools.product(config["beta"], config["m"], config["lambda"]):

        FitResults = plotMarkers_getFits_uniI(config, pars, "barI", ax_up, ax_down)

    FitResults.update(plotMarkers_getFits_uniI(config, [0,0,0], "barI", ax_up, ax_down))

    for k, v in FitResults.items():
        print("\n\t {} \n\t\t {}".format(k, v))

    handle_EOSnames = [
        Line2D(
            [], [],
            color = "k",
            marker = map_ms.get( _, None),
            linewidth = 0,
            linestyle = None,
            label = _
        ) for _ in config["eos_name"]
    ]

    handle_linestyles = [
        Line2D(
            [], [],
            color = "k",
            marker = None,
            linestyle = map_ls.get(_, None),
            label = (
                "$\lambda = {}$".format(_)
                if _ != "GR" else "GR" )
        ) for _ in [ *config["lambda"], "GR" ]
    ]

    handle_colors = [
        mpatches.Patch(
            color = map_c.get(_, None),
            label = (
                "$m= $ {}".format(convert_NumScientific(_))
                if _ != "GR" else "GR"
            )
        ) for _ in [ *config["m"], "GR" ]
    ]

    #~ handle_colors = [

        #~ Line2D(
            #~ [], [],
            #~ color = map_c.get(_, None),
            #~ marker = None,
            #~ label = (
                #~ "$\lambda = {}$".format(_)
                #~ if _ != "GR" else "GR" )
        #~ ) for _ in [ *config["lambda"], "GR" ]

    #~ ]

    ax_up.add_artist( ax_up.legend(
        handles = [ *handle_EOSnames ],
        loc = 2,
        ncol = 2,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
        #~ bbox_to_anchor=(1, 1),
        #~ borderaxespad=0.1
    ) )

    ax_up.add_artist( ax_up.legend(
        handles = [ *handle_linestyles ],
        loc = 7,
        ncol = 1,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
        #~ bbox_to_anchor=(0.79, 0.5),
        #~ borderaxespad=0.1
    ) )

    ax_up.add_artist( ax_up.legend(
        handles = [ *handle_colors ],
        loc = 4,
        ncol = 1,
        frameon = True,
        markerscale = 0.6,
        fancybox=True,
        framealpha = 0.5,
        handlelength = 4,
        #~ bbox_to_anchor=(0.99, 0.22),
        #~ borderaxespad=0.1
    ) )

    # lets add a legend for marker styles

    ax_up.set_xlim(0.09, 0.325)
    ax_up.set_ylim(5, 35)

    ax_down.set_ylim(1e-3, 1.5e0)
    #~ ax_down.axhline(y=0.1, linewidth=2, color='r', alpha = 0.5)

    plt.savefig(
        'barI_GR.eps', format="eps",
        dpi=600,
        pad_inches=0,
        bbox_inches = "tight",
        papertype = "a4"
    )

    plt.show()

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

        # the config file
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

        # MvsR with GR included
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

        # creating all uni I in separated directories
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

        # create tildeI_GR
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

        # create tildeI_GR
        parser.add_argument(
            "--tildeI_GR_zoomBox",
            action = "store_const",
            #~ nargs = "?",
            #~ type = str,
            #~ default = "quick_plotter.config",
            const = "tildeI_GR_zoomBox",
            #~ metavar = "",
            dest = "ConfigSection",
            required = False,
            help = """
                it will use the tilde I section
            """
        )

        # create barI_GR
        parser.add_argument(
            "--barI_GR",
            action = "store_const",
            #~ nargs = "?",
            #~ type = str,
            #~ default = "quick_plotter.config",
            const = "barI_GR",
            #~ metavar = "",
            dest = "ConfigSection",
            required = False,
            help = """
                if called, it will use barI_GR section in config file to produce
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

    elif args.ConfigSection == "tildeI_GR":

        plot_tildeI_GR(config.items(args.ConfigSection))

    elif args.ConfigSection == "tildeI_GR_zoomBox":
        #TODO Documentation for this needed

        plot_tildeI_GR_zoomBox(config.items("tildeI_GR"))

    elif args.ConfigSection == "barI_GR":

        plot_barI_GR(config.items(args.ConfigSection))

    else:
        print("\n {} unknown, terminating... \n", args.ConfigSection)
