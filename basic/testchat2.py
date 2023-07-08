import math

def INIT(fq00, rr00, th00, ph00, es00, dl00):
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD

    RE = 6370.0
    RREF = 7370.0
    RREFP = RREF
    RBTM = 9370.0
    CVLCTY = 3.00e5

    XME = 9.109534e-31
    XMP = 1.6726485e-27
    QE = 1.6021892e-19
    EP = 8.854185e-12
    XK = 1.380662e-23
    GM = 3.98603e14

    PF0 = [
        QE**2 / (4.0 * PAI**2 * EP * XME),
        QE**2 / (4.0 * PAI**2 * EP * XMP),
        QE**2 / (4.0 * PAI**2 * EP * (XMP * 4.0)),
        QE**2 / (4.0 * PAI**2 * EP * (XMP * 16.0))
    ]

    GF0 = [
        -QE / (2.0 * PAI * XME) * 1.0e-3,
        QE / (2.0 * PAI * XMP) * 1.0e-3,
        QE / (2.0 * PAI * (XMP * 4.0)) * 1.0e-3,
        QE / (2.0 * PAI * (XMP * 16.0)) * 1.0e-3
    ]

    XN0 = 50000.0
    THERM = 1000.0
    SHINV = [
        XMP * (GM / RREF**2) / XK / THERM * 1.0e-3,
        XMP * 4.0 * (GM / RREF**2) / XK / THERM * 1.0e-3,
        XMP * 16.0 * (GM / RREF**2) / XK / THERM * 1.0e-3
    ]
    SHINVP = SHINV[0]
    ETA = [0.0, 0.152, 0.82, 0.025]

    RRLI = 6370.0 + 90.0
    HWLI = 140.0
    XLPP = 4.0
    ERPP = 0.2
    DRR1 = 5.0
    DTH1 = 0.05
    DPH1 = 0.2
    DRR2 = 2.5
    DTH2 = 0.01
    DPH2 = 0.1

    RRMAX = 4.0 * RE
    RRMIN = 6370.0
    ERRR = 1.0e-5
    ERRA = 1.0e-5
    HMIN = 1.0e-10
    DISOUT = 200
    NLPMAX = 30000

    INDOUT = 21
    INDMSH = 51
    IEOF = -1

    F0 = fq00
    R0 = rr00
    T0 = th00
    P0 = ph00
    E0 = es00
    D0 = dl00

    return (F0, R0, T0, P0, E0, D0)

import math

def DIPOLE(PATH):
    # Constants
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD
    RE = 6370.0
    CVLCTY = 3.00e5
    GMDM = 8.07e6

    # Common block variables
    RR = PATH[0]
    TH = PATH[1]
    PH = PATH[2]
    COSTH = math.cos(TH * RAD)
    SINTH = math.sin(TH * RAD)
    COSTHS = 1.0 + 3.0 * COSTH ** 2
    YR = -2.0 * COSTH / math.sqrt(COSTHS)
    YT = -SINTH / math.sqrt(COSTHS)
    YP = 0.0

    DYRDRR = 0.0
    DYTDRR = 0.0
    DYPDRR = 0.0
    DYRDTH = -2.0 / COSTHS * YT
    DYTDTH = 2.0 / COSTHS * YR
    DYPDTH = 0.0
    DYRDPH = 0.0
    DYTDPH = 0.0
    DYPDPH = 0.0

    BB = GMDM * math.sqrt(COSTHS) / RR ** 3

    DBBDRR = -3.0 / RR
    DBBDTH = -3.0 * SINTH * COSTH / COSTHS
    DBBDPH = 0.0

    return

import math

def DENS(PATH):
    # Constants
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD
    RE = 6370.0
    CVLCTY = 3.00e5

    # Common block variables
    NSPEC = 4
    XNS = [0.0] * NSPEC
    DNSDRR = [0.0] * NSPEC
    DNSDTH = [0.0] * NSPEC
    DNSDPH = [0.0] * NSPEC
    XNDE = [0.0] * NSPEC
    DNDDRR = [0.0] * NSPEC
    XNLI = 0.0
    DNLDRR = 0.0
    XNAO = 0.0
    DNADRR = 0.0
    DNADTH = 0.0
    DNADPH = 0.0
    XNOP = 0.0
    DNODRR = 0.0
    DNODTH = 0.0
    DNODPH = 0.0

    FLPRM(PATH)
    DENSDE(PATH)
    DENSLI(PATH)
    DENSAO(PATH)
    DENSOP(PATH)

    for IS in range(1, NSPEC + 1):
        XNS[IS - 1] = XNDE[IS - 1] * XNLI * XNAO * XNOP

        DNSDRR[IS - 1] = DNDDRR[IS - 1] + DNLDRR + DNADRR + DNODRR
        DNSDTH[IS - 1] = DNADTH + DNODTH
        DNSDPH[IS - 1] = DNADPH + DNODPH

    return


def DENSAO(PATH):
    # Constants
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD
    RE = 6370.0
    CVLCTY = 3.00e5

    # Common block variables
    XNAO = 1.0
    DNADRR = 0.0
    DNADTH = 0.0
    DNADPH = 0.0

    return


def DENSDE(PATH):
    # Constants
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD
    RE = 6370.0
    CVLCTY = 3.00e5

    # Common block variables
    NSPEC = 4
    XNDE = [0.0] * NSPEC
    DNDDRR = [0.0] * NSPEC

    GPH = RREF * (PATH[0] - RREF) / PATH[0]
    DGPHDR = (RREF / PATH[0]) ** 2
    ETEXP2 = ETA[1] * math.exp(-GPH * SHINV[1])
    ETEXP3 = ETA[2] * math.exp(-GPH * SHINV[2])
    ETEXP4 = ETA[3] * math.exp(-GPH * SHINV[3])

    SUM1 = ETEXP2 + ETEXP3 + ETEXP4
    SUM2 = ETEXP2 * SHINV[1] + ETEXP3 * SHINV[2] + ETEXP4 * SHINV[3]

    XNDE[0] = math.sqrt(SUM1)
    XNDE[1] = ETEXP2 / XNDE[0]
    XNDE[2] = ETEXP3 / XNDE[0]
    XNDE[3] = ETEXP4 / XNDE[0]

    DNDDRR[0] = -DGPHDR * SUM2 / SUM1 / 2.0
    DNDDRR[1] = -DGPHDR * SHINV[1] - DNDDRR[0]
    DNDDRR[2] = -DGPHDR * SHINV[2] - DNDDRR[0]
    DNDDRR[3] = -DGPHDR * SHINV[3] - DNDDRR[0]

    return


def DENSLI(PATH):
    # Constants
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD
    RE = 6370.0
    CVLCTY = 3.00e5

    # Common block variables
    XNLI = 0.0
    DNLDRR = 0.0

    if PATH[0] - RRLI >= 0:
        ALTNRM = (PATH[0] - RRLI) / HWLI
        XNLI = 1.0 - math.exp(-ALTNRM ** 2)
        DNLDRR = math.exp(-ALTNRM ** 2) * ALTNRM * 2.0 / HWLI
    else:
        XNLI = 0.0
        DNLDRR = 0.0

    return


def DENSOP(PATH):
    # Constants
    PAI = math.asin(1.0) * 2.0
    RAD = PAI / 180.0
    DEG = 1.0 / RAD
    RE = 6370.0
    CVLCTY = 3.00e5

    # Common block variables
    XNOP = XN0
    DNODRR = 0.0
    DNODTH = 0.0
    DNODPH = 0.0

    return
