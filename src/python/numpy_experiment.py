#!/usr/bin/env python

from math import exp
from math import log10
from math import fabs
from math import floor
from math import
tiny = 1e-30
min_p_power = -16
units = units_coef_clac
EOSdata = {
    "K": [],
    "Gamma": [],
    "a": [],
    "density": [],
    "pressure": [],
    "N": []
}

def EOS_init():

    log_p1 = 34.269;
    Gamma1 = 2.830;
    Gamma2 = 3.445;
    Gamma3 = 3.348;

    EOSdata["a"].append(0)

    EOSdata["K"].append(6.80110e-9)
    EOSdata["Gamma"].append(1.58425)
    EOSdata["density"].append(2.44034e7)

    EOSdata["K"].append(1.06186e-6)
    EOSdata["Gamma"].append(1.28733)
    EOSdata["density"].append(3.78358e11)

    EOSdata["K"].append(5.32697e+1)
    EOSdata["Gamma"].append(0.62223)
    EOSdata["density"].append(2.62780e12)

    EOSdata["K"].append(3.99874e-8)
    EOSdata["Gamma"].append(1.35692)
    EOSdata["density"].append(0.00000e00)

    EOSdata["K"] = [ i*units["c"]**2 for i in EOSdata["K"] ]

    EOSdata["K"].append(0.00000e00)
    EOSdata["Gamma"].append(Gamma1)
    EOSdata["density"].append(10**14.7)

    EOSdata["K"].append(0.00000e00)
    EOSdata["Gamma"].append(Gamma2)
    EOSdata["density"].append(1.00000e15)

    EOSdata["K"].append(0.00000e00)
    EOSdata["Gamma"].append(Gamma3)
    EOSdata["density"].append(9.99999e99)

    EOSdata["K"][4] = 10**log_p1 / EOSdata["density"][4]**EOSdata["Gamma"][4]

    if EOSdata["Gamma"][4] == EOSdata["Gamma"][3]:
        EOSdata["density"][3] = 1e14
    else:
        EOSdata["density"][3] = \
          (EOSdata["K"][3]/EOSdata["K"][4])**(1/(EOSdata["Gamma"][4]-EOSdata["Gamma"][3]))

    for i in range(5,7):
        EOSdata["K"][i] = \
            EOSdata["K"][i-1]
            * EOSdata["density"][i-1]**(
                EOSdata["Gamma"][i-1] - EOSdata["Gamma"][i]
            )

    for i in range(1,7):
        EOSdata["a"][i] = \
            EOSdata["a"][i-1]
            + EOSdata["K"][i-1] / units["c"]**2 / (EOSdata["Gamma"][i-1] - 1)
            * EOSdata["density"][i-1]**(EOSdata["Gamma"][i-1]-1)
            - EOSdata["K"][i] / units["c"]**2 / (EOSdata["Gamma"][i] - 1)
            * EOSdata["density"][i-1]**(EOSdata["Gamma"][i] - 1)

    EOS_data["p"] = [
        EOSdata["K"][i]*EOSdata["density"][i]**EOSdata["Gamma"][i]
        for i in range(0,7)
    ]

    return

def get_power(num):

    return floor(log10(fabs(num + 1e-30)))

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

        units["c"] = const_c

        return units

def EOS_eq(p):

    p *= units["pressure"]
    index = 0

    while p > EOS_data["pressure"][i]:
        i++

    eps = \
        ( 1 + EOS_data["a"][i])
        * (p/EOS_data["K"][i])**(1/EOS_data["Gamma"][i])
        + p/EOS_data["c"]**2 / (EOS_data["Gamma"][i] - 1)

    return eps/units["density"]

def foo(y, t, beta_phiScal, beta_phiScal, lambda_phiScal):

    r = t

    phiScal, Q, p, LambdaMetr, m = y

    Vhat = 2 * m_phiScal**2 * phiScal**2 + lambda_phiScal * phiScal**4
    Vhat_dphiScal = 4 * m_phiScal**2 * phiScal
        + 4 * lambda_phiScal * phiScal**3

    alpha = beta_phiScal*phiScal
    A = exp(1/2 * beta_phiScal * phiScal**2)
    exp_2LambdaMetr = exp(2*LambdaMetr)

    step4_A = A**4
    rho = 0

    if not R and p and get_power(p) <= min_p_power:
        R = r
        AR = A*r
        p = 0
        rho = 0
    elif get_power(p) <= min_p_power:
        p = 0
        rho = 0
    else:
        rho = EOS_eq(p)

    PhiMetr_dr = \
      r*1/2*(
        8*pi*step4_A*rho*exp_2LambdaMetr
        + Q**2
        + 1/2*exp_2LambdaMetr*Vhat
        + (1/r)**2 * (exp_2LambdaMetr-1)
      )

    phiScal_dr = Q

    LambdaMetr_dr = \
      r*1/2*(
        8*pi*step4_A*rho*exp_2LambdaMetr
        +Q**2
        +1/2*exp_2LambdaMetr*Vhat
        - (1/r)**2*(exp_2LambdaMetr-1)
      )

    Q_dr = \
      4*pi*alpha*step4_A*(rho - 3*p)*exp_2LambdaMetr
      + 1/4*Vhat_dphiScal*exp_2LambdaMetr
      - (PhiMetr_dr - LambdaMetr_dr + 2/r )*Q

    p_dr = -(rho + p)*(PhiMetr_dr + alpha*Q)

    m_dr = \
      r**2*(
        4*pi*step4_A*rho
        + 1/2*1/exp_2LambdaMetr*Q**2
        + 1/4*Vhat
      )

    dydx = []
    dydx.append(phiScal_dr)
    dydx.append(Q_dr)
    dydx.append(p_dr)
    dydx.append(LambdaMetr_dr)
    dydx.append(m_dr)

    return dydx

def descriptancy(v, beta_phiScal, beta_phiScal, lambda_phiScal):

    R = 0
    AR = 0

    y0 = [ v, 0, 3e-4, 0, 0 ]

    result =

EOS_init()



