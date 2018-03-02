#!/usr/bin/env python

file_path = "/home/dimitar/projects/STT_theories/results"
file_name_results = "STT_phiScal_I"

system_names = [ "r", "phiScal", "Q", "p", "LambdaMetr", "m" ]

parameters = {
    "names" : [ "beta", "m", "lambda" ],
    "values" : {
        "beta" : [ 0, -6, -10 ],
        "m" : [ 0, 1e-3, 5e-2 ],
        "lambda" : [ 0, 1e-1, 1, 10, 100 ]
    }
}

eos_names = [ "EOSII", "AkmalPR" ]

def my_ploting(ax, label_x, x, label_y, y):

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
        ax.plot(
          x, y,
          linewidth=1.5
        )
    except ValueError:
        pass

def animate():

    import glob
    import os

    all_files = glob.glob(file_path + file_name_results + "*")
    while not len(all_files):
        all_files = glob.glob(file_path + file_name_results + "*")

    file_to_use = max(all_files, key=os.path.getctime)

    with open(file_to_use,'r') as f:
        graph_data = f.read()

    graph_data = [ k.split("\n") for k in graph_data.split("# ") if len(k) > 1 ][1]

    plot_labels = graph_data.pop(0).split(" ")

    plots_all = [ [] for k in range(len(plot_labels)) ]

    try:
        for single_line in graph_data:
            if len(single_line) > 1:
                tmp = [ j for j in single_line.split(" ") if len(j) > 1 ]
                for m,n in zip(plots_all, tmp):
                    m.append(float(n))
    except ValueError:
        graph_data = []
        plot_labels = []
        plots_all = []
        pass

    my_ploting(
      ax_phiScal,
      plot_labels[0], plots_all[0],
      plot_labels[1], plots_all[1],
    )

    my_ploting(
      ax_mR,
      plot_labels[3], plots_all[3],
      plot_labels[2], plots_all[2],
    )

    #ax_phiScal.set_title(file_to_use.split("_")[3:])
    plt.suptitle(file_to_use.split("_")[3:], fontsize=16, y=1.001)

    graph_data = []
    plot_labels = []
    plots_all = []

def get_parm_order():

    all_choices = {
        1 : {
            "change" : parameters["names"][0],
            "fixed" : [ parameters["names"][1], parameters["names"][2] ],
            "file" : "{f_name:}_{eos_name:}_beta{change:.3e}_m{fix1:.3e}_lambda{fix2:.3e}"
        },

        2 : {
            "change" : parameters["names"][1],
            "fixed" : [ parameters["names"][0], parameters["names"][2] ],
            "file" : "{f_name:}_{eos_name:}_beta{fix1:.3e}_m{change:.3e}_lambda{fix2:.3e}"
        },

        3 : {
            "change" : parameters["names"][2],
            "fixed" : [ parameters["names"][0], parameters["names"][1] ],
            "file" : "{f_name:}_{eos_name:}_beta{fix1:.3e}_m{fix2.3e}_lambda{change:.3e}"
        }
    }

    print("\n Choose graph type: \n")
    for i in all_choices.keys():
        print(
          "{}. change {}; \t fixed {} and {}".format(
            i, all_choices[i]["change"],
            all_choices[i]["fixed"][0], all_choices[i]["fixed"][1]
          )
        )

    while True:
        try:
            choice = int(input("\n Type index of parameter order:... "))

        except ValueError:
            print("{} not valid, try again... \n".format(choice))
            continue

        else:
            if choice > 3 or choice < 1:
                print("{} not valid, try again... \n".format(choice))
                continue
            else:
                return all_choices[choice]

def set_fix_parm(fix_parm_name):

    print("\n Set parameter {}: \n".format(fix_parm_name))
    for i, value in enumerate(parameters["values"][fix_parm_name]):
        print(
          "\t {}. {:.3e}".format(i, value)
        )

    while True:
        try:
            choice = int(input("\n Type index of value:... "))

        except ValueError:
            print("{} not valid, try again... \n".format(choice))
            continue

        else:
            if choice > len(parameters["values"][fix_parm_name]) or choice < 0:
                print("{} not valid, try again... \n".format(choice))
                continue
            else:
                return {
                  "name" : fix_parm_name,
                  "value" : parameters["values"][fix_parm_name][choice]
                }

def set_change_parm(change_parm_name):

    print("\n Choose parameter range {}: \n".format(change_parm_name))
    for i, value in enumerate(parameters["values"][change_parm_name]):
        print(
          "\t {}. {:.3e}".format(i, value)
        )

    while True:
        try:
            choice1 = int(input("\n Start value index:... "))

        except ValueError:
            print("{} not valid, try again... \n".format(choice1))
            continue

        else:
            if choice1 > len(parameters["values"][change_parm_name]) or choice1 < 0:
                print("{} not valid, try again... \n".format(choice1))
                continue
            else:
                break

    while True:
        try:
            choice2 = int(input("\n End value index:... "))

        except ValueError:
            print("{} not valid, try again... \n".format(choice2))
            continue

        else:
            if choice2 > len(parameters["values"][change_parm_name]) or choice2 < 0:
                print("{} not valid, try again... \n".format(choice1))
            else:
                break

    return {
      "name" : change_parm_name,
       "value" : parameters["values"][change_parm_name][choice1 : choice2+1]
    }

def set_eos_name():

    print("\n Choose EOS name: \n")
    for i, value in enumerate(eos_names):
        print(
          "\t {}. {}".format(i, value)
        )

    while True:
        try:
            choice = int(input("\n EOS name:... "))

        except ValueError:
            print("{} not valid, try again... \n".format(choice))
            continue

        else:
            if choice > len(eos_names) or choice < 0:
                print("{} not valid, try again... \n".format(choice))
                continue
            else:
                return \
                  { "name" : eos_names[choice] }

def load_data(filename, data = {}):

    data.clear()

    with open(filename, "r") as f:
        tmp_data = f.read()

    tmp_data = tmp_data.split("#")[1].split("\n")

    data["labels"] = [ k for k in tmp_data.pop(0).split(" ") if len(k) ]

    data["data"] = [ [] for k in range(len(data["labels"])) ]

    for i in tmp_data:
        for m, n in zip( i.split(" "), data["data"] ):
            if len(m):
                n.append(float(m))

    return

def get_figure_phiScal_mR():

    from matplotlib import style
    from matplotlib.gridspec import GridSpec

    style.use("seaborn-poster")
    gs = GridSpec(nrows=1,ncols=2)

    fig = plt.figure()
    fig.set_tight_layout(True)

    ax_phiScal = fig.add_subplot(gs[0,0])
    ax_mR = fig.add_subplot(gs[0,1])

    return ax_phiScal, ax_mR

def my_plot_set_params( ax, label_x, label_y, label_fontsize, label_ticksize ):

    from matplotlib.ticker import FormatStrFormatter

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_xlabel(label_x, fontsize=label_fontsize)
    ax.xaxis.set_tick_params(labelsize=label_ticksize)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_ylabel(label_y, fontsize=label_fontsize)
    ax.yaxis.set_tick_params(labelsize=label_ticksize)

    return

def my_plot_set_data( ax, x, y, linewidth):

    ax.plot(x, y, linewidth = linewidth)

if __name__ == "__main__":

    from matplotlib import pyplot as plt

    parm_order = get_parm_order()

    fix_parm_1 = set_fix_parm(parm_order["fixed"][0])

    fix_parm_2 = set_fix_parm(parm_order["fixed"][1])

    change_parm = set_change_parm(parm_order["change"])

    eos_name_val = set_eos_name()["name"]

    f_name_val = file_name_results
    fix1_val = fix_parm_1["value"]
    fix2_val = fix_parm_2["value"]

    ax_phiScal, ax_mR = get_figure_phiScal_mR()

    data = {}
    for change_val in change_parm["value"]:

        current_file = parm_order["file"].format(
          f_name = f_name_val, eos_name = eos_name_val,
          fix1 = fix1_val, fix2 = fix2_val, change = change_val
        )

        print( current_file )
        load_data(file_path + "/" +current_file, data)

        my_plot_set_params(ax_phiScal, data["labels"][0], data["labels"][1], 12, 10)
        my_plot_set_params(ax_mR, data["labels"][5], data["labels"][4], 12, 10)

        my_plot_set_data( ax_phiScal, data["data"][0], data["data"][1], 1.5)
        my_plot_set_data( ax_mR, data["data"][5], data["data"][4], 1.5)


    plt.show()
