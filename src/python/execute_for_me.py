#!/usr/bin/env python

import sys
import shutil
import os
import pathlib
import time
import itertools
import glob

from subprocess import Popen, PIPE, STDOUT

#~ path to make file
MAKEFILE_path = "/home/dimitar/projects/STT_theories/"

#~ path to were the resuts are kept
RESULTS_path = "/home/dimitar/projects/STT_theories/results"
RESULTS_fname = "STT_phiScal_J"

EOS_start_n = 2

EOS_all = [
    "", "EOSII", "PAL6", "SLy", "APR1", "APR2", "APR3", "APR4", "FPS", "WFF1",
    "WFF2", "WFF3", "BBB2", "BPAL12", "ENG", "MPA1", "MS1", "MS2", "MS1b",
    "PS", "GS1", "GS2", "BGN1H1", "GNH3", "H1", "H2", "H3", "H4", "H5",
    "H6", "H7", "PCL2", "ALF1", "ALF2", "ALF3", "ALF4"
]

#~ only those who have max mass in units of the sun of 2
#~ the list i would like to use
#~ EOS_max18 = [
    #~ "SLy", "APR2", "APR3", "APR4", "FPS", "WFF1", "WFF2", "WFF3", "BBB2", "ENG",
    #~ "MPA1", "MS1", "MS2",  "MS1b", "PS", "GNH3", "H3", "H4", "ALF2", "ALF4"
#~ ]
EOS_max18 = [
    "PS", "GNH3", "H3", "H4", "ALF2", "ALF4"
]
#~ EOS_max2 = [
    #~ "SLy", "APR3", "APR4", "WFF1", "WFF2", "ENG",
    #~ "MPA1", "MS1", "MS1", "MS1b", "H4", "ALF2"
#~ ]

#~ to avoid rewriting the code in C and give me opportunity to be more flexible
#~ will only map the ones I wan to use to all of them as nested list of dictionary
EOS_mapping = [
    EOS_all.index(_) for _ in EOS_max18
]

#~ the source file for EOS
EOS_file = "/home/dimitar/projects/STT_theories/src/EOS.c"

#~ the line number, minus one since the list starts at 0
EOS_line = 9 - 1
#~ the substitute of the line parametrized
EOS_sub = "#define DEF_EOS_MODEL_NUM {:.0f}\n"

#~ path to the ODE_phiScal_J source file
ODE_PATH = "/home/dimitar/projects/STT_theories/src/ODE_phiScal_J.c"
#~ the line numbers for the value of the guess and initial infinity
P_C_line = 7 - 1
P_C_sub = "#define P_START {:.2e} // 5e-5; GR 1e-5\n"
P_C_end_line = 8 - 1

PHISCAL_C_line = 9 - 1
PHISCAL_C_sub = "#define GVPHISCAL_FF {:.2e}\n"

R_INF_line = 12 - 1
R_INF_sub = "#define R_INF_PHISCAL {:.2e}\n"

BETA_line = 16 - 1
BETA_sub = "#define BETA {:.2e}\n"

M_line = 17 - 1
M_sub = "#define M {:.2e}\n"

LAMBDA_line = 18 - 1
LAMBDA_sub = "#define LAMBDA {:.2e}\n"

EOS_content = None
with open(EOS_file, "r") as f:
    EOS_content = f.readlines()

ODE_content = None
with open(ODE_PATH, "r") as f:
    ODE_content = f.readlines()

P_C_init = float(ODE_content[P_C_line].strip().split()[2])

P_C_end = float(ODE_content[P_C_end_line].strip().split()[2])

PHISCAL_C_init = float(ODE_content[PHISCAL_C_line].strip().split()[2])

R_INF_init = float(ODE_content[R_INF_line].strip().split()[2])

BETA_init = float(ODE_content[BETA_line].strip().split()[2])

M_init = float(ODE_content[M_line].strip().split()[2])

LAMBDA_init = float(ODE_content[LAMBDA_line].strip().split()[2])

#~ for EOS_cur in range(EOS_start_n, len(EOS_all)):
for EOS_cur in EOS_mapping:

    EOS_content[EOS_line] = EOS_sub.format(EOS_cur)

    with open(EOS_file, "w") as f:
        f.writelines(EOS_content)

    #~ ODE_content[P_C_line] = P_C_sub.format(P_C_init)
    #~ ODE_content[PHISCAL_C_line] = PHISCAL_C_sub.format(PHISCAL_C_init)
    #~ ODE_content[R_INF_line] = R_INF_sub.format(R_INF_init)
    #~ ODE_content[BETA_line] = BETA_sub.format(BETA_init)
    #~ ODE_content[M_line] = M_sub.format(M_init)
    #~ ODE_content[LAMBDA_line] = LAMBDA_sub.format(LAMBDA_init)

    #~ with open(ODE_PATH, "w") as f:
        #~ f.writelines(ODE_content)

    all_beta = [ -6 ]
    all_m = [ 0, 5e-3, 1e-2, 5e-2 ]
    all_lambda = [ 0, 1e-1, 1e0, 1e1 ]

    for p_beta, p_m, p_lambda in itertools.product(
        all_beta,all_m, all_lambda, repeat=1
    ):

        print(
            "\n Starting "
            "\n eos N. name: {}. {}"
            "\n beta = {:.3e}; m = {:.3e}; lambda = {:.3e}; \n".format(
                EOS_cur, EOS_all[EOS_cur], p_beta, p_m, p_lambda
            )
        )

        ODE_content[BETA_line] = BETA_sub.format(p_beta)
        ODE_content[M_line] = M_sub.format(p_m)
        ODE_content[LAMBDA_line] = LAMBDA_sub.format(p_lambda)

        with open(ODE_PATH, "w") as f:
            f.writelines(ODE_content)

        r_boom = None
        last_p = None
        last_phiScal = None
        amount_files = 0

        filename_base = "_".join( [
            RESULTS_fname,
            EOS_all[EOS_cur],
            "beta{:.3e}".format(p_beta),
            "m{:.3e}".format(p_m),
            "lambda{:.3e}".format(p_lambda)
        ] )

        while not r_boom:

            with Popen(
                ["make", "run", "-C", MAKEFILE_path],
                stdout=PIPE,
                stderr=STDOUT,
                bufsize=1
            ) as p:

                for line in p.stdout:

                    if "integrate_phiscal says system went boom at" in line.decode():

                        r_boom = float(line.decode().strip().split()[-1])

                    elif "p_c" in line.decode():

                        last_p = float(line.decode().strip().split()[2])
                        #~ print(
                            #~ "\n\t p_c = {:.2e}, v_init = {:.2e}, r_inf = {:.2e}".format(

                            #~ last_p,
                            #~ float(last_phiScal if last_phiScal else ODE_content[PHISCAL_C_line].strip().split()[-1]),
                            #~ float(ODE_content[R_INF_line].strip().split()[-1])
                        #~ ) )

                    elif "setting v_phiScal" in line.decode():

                        last_phiScal = float(line.decode().strip().split()[3])

                        #~ print("\n\t\t last phiScal_c = {:.2e}".format(last_phiScal))

                        if abs(last_phiScal) <= 1e-2:

                            last_phiScal = None

            #~ if there was a boom and we are not withing the 95% range of final central pressure
            if r_boom and last_p < P_C_end*0.92:

                amount_files += 1

                result_cur = os.path.join(
                    RESULTS_path,
                    filename_base
                )

                result_save = result_cur + "_{:.0f}".format(amount_files)

                #~ print(
                    #~ "\n\n\t BOOM HAPPEND \t at {:.6e} \t try {:.6e}\n"
                    #~ "\n\t save current progress \n\t\t at {}\n".format(
                        #~ r_boom, r_boom*0.95, result_save
                    #~ )
                #~ )

                shutil.move(result_cur, result_save)

                if last_p:
                    ODE_content[P_C_line] = P_C_sub.format(last_p)

                if last_phiScal:
                    ODE_content[PHISCAL_C_line] = PHISCAL_C_sub.format(last_phiScal)

                ODE_content[R_INF_line] = R_INF_sub.format(r_boom*0.95)

                with open(ODE_PATH, "w") as f:
                    f.writelines(ODE_content)

                r_boom = None
                last_p = None
                last_phiScal = None

            else:
                r_boom = True

        result_full_cur = os.path.join( RESULTS_path, EOS_all[EOS_cur] )

        pathlib.Path( result_full_cur ).mkdir(parents=True, exist_ok=True)

        #~ print("\n will write out {} \n".format(
            #~ os.path.join( result_full_cur, filename_base )
        #~ ) )

        if amount_files:

            with open( os.path.join( result_full_cur, filename_base ), "w") as f:

                for _ in [
                    *[
                        os.path.join( RESULTS_path, filename_base + "_{:.0f}".format(_) )
                        for _ in range(1, amount_files+1)
                    ],
                    os.path.join(RESULTS_path, filename_base)
                ]:

                    #~ print("\n\t adding {} \n".format(
                        #~ _
                    #~ ) )

                    with open(_, "r") as finside:
                        f.write(finside.read())
        else:

            shutil.move(
                os.path.join(RESULTS_path, filename_base),
                os.path.join(RESULTS_path, EOS_all[EOS_cur], filename_base)
            )

        _file_data = None
        with open(os.path.join(RESULTS_path, EOS_all[EOS_cur], filename_base), "r") as f:
            _file_data = f.readlines()

        _file_data = [
            __ for _, __ in enumerate(_file_data) if _ == 0 or __[0] != "#"
        ]

        with open(os.path.join(RESULTS_path, EOS_all[EOS_cur], filename_base), "w") as f:
            f.writelines(_file_data)

        #~ RESET THEM ALL BECAUSE I WANT SO
        ODE_content[P_C_line] = P_C_sub.format(P_C_init)
        ODE_content[PHISCAL_C_line] = PHISCAL_C_sub.format(PHISCAL_C_init)
        ODE_content[R_INF_line] = R_INF_sub.format(R_INF_init)
        ODE_content[BETA_line] = BETA_sub.format(BETA_init)
        ODE_content[M_line] = M_sub.format(M_init)
        ODE_content[LAMBDA_line] = LAMBDA_sub.format(LAMBDA_init)

        with open(ODE_PATH, "w") as f:
            f.writelines(ODE_content)

        for _ in glob.glob(os.path.join(RESULTS_path, "live_*")):
            os.remove(_)

        for _ in glob.glob(os.path.join(RESULTS_path, "STT_phiScal_J*")):
            os.remove(_)
