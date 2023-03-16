#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 16 4:06 PM 2023
Created in PyCharm
Created as QGP_Scripts/gaus_dsigma_fits

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit as cf

from Measure import Measure


def main():
    sigma = 0.3
    x, y = get_data(sigma)
    # plot_data(x, y)

    func = cos_1_sin
    p0 = [1, 0.2]

    fit_data(x, y, func, p0)
    plt.show()
    print('donzo')


def fit_data(x, y, func, p0):
    popt, pcov = cf(func, x, y, p0=p0)

    plot_data(x, np.array(y) - func(x, *popt), title='Fit Residuals', y_label=r'$\Delta\sigma$ Residuals')

    xs_fit = np.linspace(min(x), max(x), 1000)
    plot_data(x, y)
    plt.plot(xs_fit, func(xs_fit, *popt), color='red')

    print(', '.join([str(Measure(val, err)) for val, err in zip(popt, np.sqrt(np.diag(pcov)))]))


def plot_data(x, y, title='Observable', y_label=r'$\Delta\sigma$'):
    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='gray', zorder=0)
    ax.scatter(x, y)
    ax.set_xlabel('Azimuthal Partition Width w')
    ax.set_ylabel(y_label)
    ax.set_title(title)
    fig.tight_layout()


def get_data(sigma=0.8):
    if sigma == 0.3:
        x = np.array([0.0, 0.06346651825433926, 0.12693303650867852, 0.1903995547630178, 0.25386607301735703, 0.3173325912716963, 0.3807991095260356, 0.4442656277803748, 0.5077321460347141, 0.5711986642890533, 0.6346651825433925, 0.6981317007977318, 0.7615982190520711, 0.8250647373064104, 0.8885312555607496, 0.9519977738150889, 1.0154642920694281, 1.0789308103237674, 1.1423973285781066, 1.2058638468324459, 1.269330365086785, 1.3327968833411243, 1.3962634015954636, 1.4597299198498028, 1.5231964381041423, 1.5866629563584815, 1.6501294746128208, 1.71359599286716, 1.7770625111214993, 1.8405290293758385, 1.9039955476301778, 1.967462065884517, 2.0309285841388562, 2.0943951023931957, 2.1578616206475347, 2.221328138901874, 2.284794657156213, 2.3482611754105527, 2.4117276936648917, 2.475194211919231, 2.53866073017357, 2.6021272484279097, 2.6655937666822487, 2.729060284936588, 2.792526803190927, 2.8559933214452666, 2.9194598396996057, 2.982926357953945, 3.0463928762082846, 3.1098593944626236, 3.173325912716963, 3.236792430971302, 3.3002589492256416, 3.3637254674799806, 3.42719198573432, 3.490658503988659, 3.5541250222429985, 3.6175915404973376, 3.681058058751677, 3.744524577006016, 3.8079910952603555, 3.8714576135146945, 3.934924131769034, 3.998390650023373, 4.0618571682777125, 4.1253236865320515, 4.188790204786391, 4.25225672304073, 4.3157232412950695, 4.3791897595494085, 4.442656277803748, 4.506122796058087, 4.569589314312426, 4.6330558325667655, 4.696522350821105, 4.759988869075444, 4.823455387329783, 4.886921905584122, 4.950388423838462, 5.013854942092801, 5.07732146034714, 5.14078797860148, 5.204254496855819, 5.267721015110158, 5.331187533364497, 5.394654051618837, 5.458120569873176, 5.521587088127515, 5.585053606381854, 5.648520124636194, 5.711986642890533, 5.775453161144872, 5.838919679399211, 5.902386197653551, 5.96585271590789, 6.029319234162229, 6.092785752416569, 6.156252270670908, 6.219718788925247, 6.283185307179586])
        y = np.array([0.0, 1.5929089091604675e-06, 6.329115969260358e-06, 1.4083868907372506e-05, 2.465829651101094e-05, 3.779116142696177e-05, 5.317719611834396e-05, 7.048097845524081e-05, 8.935608610689469e-05, 0.0001094604571465485, 0.00013046963832680745, 0.00015208681516775194, 0.00017404937386514968, 0.00019613213358683074, 0.00021814768653861255, 0.00023993422443875992, 0.0002613919562950394, 0.00028242046832687925, 0.00030295170050879255, 0.0003229361982540896, 0.0003423389304931662, 0.00036113569870452333, 0.0003793101967879964, 0.00039685170309117224, 0.00041375333676330756, 0.00042999250327076244, 0.0004456026105742611, 0.0004605643558998479, 0.0004748768139444387, 0.0004885394544702631, 0.0005015519793727941, 0.0005139142236086081, 0.0005256260965356557, 0.0005366875479740629, 0.0005470985488850999, 0.0005568369172518212, 0.0005659466939478275, 0.0005744059844728089, 0.0005822147679829526, 0.0005893730179355183, 0.0005958807013571044, 0.0006017377805129231, 0.0006069442181523832, 0.0006114999869339666, 0.0006154050820015855, 0.0006186352101209269, 0.0006212389242369609, 0.0006231921794586093, 0.0006244950982901842, 0.0006251477891459079, 0.0006251503152583382, 0.0006245026749305826, 0.0006232048015198388, 0.0006212565831574257, 0.0006186578938096043, 0.0006153829277958378, 0.0006114828915682602, 0.0006069321706489439, 0.0006017307661810856, 0.0005958787016256406, 0.0005893760101444712, 0.0005822227251343626, 0.000574418874931093, 0.0005659644810860898, 0.0005568595590320635, 0.000547077472318247, 0.0005366714475495415, 0.0005256149319825543, 0.000513907948425818, 0.0005015505401388398, 0.0004885427900250261, 0.00047488485441737893, 0.0004605770215585503, 0.00044561981043789167, 0.0004300141334711016, 0.0004137342216928852, 0.0003968373185483376, 0.0003793004431927871, 0.00036113046028418694, 0.0003423380724835168, 0.0003229395636047405, 0.0003029591060597081, 0.00028243170046760024, 0.0002614067661617092, 0.00023995232278140044, 0.00021813280516091993, 0.00019612136580504647, 0.00017404244198504504, 0.0001520833932385779, 0.0001304693480695951, 0.00010946286529245342, 8.93607037255606e-05, 7.048726280411977e-05, 5.318455544167655e-05, 3.7798962930879476e-05, 2.4653157492449118e-05, 1.4080964268026186e-05, 6.327824115470726e-06, 1.5925870928690244e-06, -2.220446049250313e-16])
    elif sigma == 0.8:
        x = np.array([0.0, 0.06346651825433926, 0.12693303650867852, 0.1903995547630178, 0.25386607301735703, 0.3173325912716963, 0.3807991095260356, 0.4442656277803748, 0.5077321460347141, 0.5711986642890533, 0.6346651825433925, 0.6981317007977318, 0.7615982190520711, 0.8250647373064104, 0.8885312555607496, 0.9519977738150889, 1.0154642920694281, 1.0789308103237674, 1.1423973285781066, 1.2058638468324459, 1.269330365086785, 1.3327968833411243, 1.3962634015954636, 1.4597299198498028, 1.5231964381041423, 1.5866629563584815, 1.6501294746128208, 1.71359599286716, 1.7770625111214993, 1.8405290293758385, 1.9039955476301778, 1.967462065884517, 2.0309285841388562, 2.0943951023931957, 2.1578616206475347, 2.221328138901874, 2.284794657156213, 2.3482611754105527, 2.4117276936648917, 2.475194211919231, 2.53866073017357, 2.6021272484279097, 2.6655937666822487, 2.729060284936588, 2.792526803190927, 2.8559933214452666, 2.9194598396996057, 2.982926357953945, 3.0463928762082846, 3.1098593944626236, 3.173325912716963, 3.236792430971302, 3.3002589492256416, 3.3637254674799806, 3.42719198573432, 3.490658503988659, 3.5541250222429985, 3.6175915404973376, 3.681058058751677, 3.744524577006016, 3.8079910952603555, 3.8714576135146945, 3.934924131769034, 3.998390650023373, 4.0618571682777125, 4.1253236865320515, 4.188790204786391, 4.25225672304073, 4.3157232412950695, 4.3791897595494085, 4.442656277803748, 4.506122796058087, 4.569589314312426, 4.6330558325667655, 4.696522350821105, 4.759988869075444, 4.823455387329783, 4.886921905584122, 4.950388423838462, 5.013854942092801, 5.07732146034714, 5.14078797860148, 5.204254496855819, 5.267721015110158, 5.331187533364497, 5.394654051618837, 5.458120569873176, 5.521587088127515, 5.585053606381854, 5.648520124636194, 5.711986642890533, 5.775453161144872, 5.838919679399211, 5.902386197653551, 5.96585271590789, 6.029319234162229, 6.092785752416569, 6.156252270670908, 6.219718788925247, 6.283185307179586])
        y = np.array([0.0, 2.348083129277752e-06, 9.378883902982698e-06, 2.1052183178039822e-05, 3.730136925561285e-05, 5.8033787689462137e-05, 8.313249860242679e-05, 0.00011245623045003962, 0.0001458413119497506, 0.00018310300486487456, 0.00022403714184932554, 0.00026842190438886374, 0.00031601971516038996, 0.0003665792180866448, 0.00041983731884767325, 0.000475514724901982, 0.0005333429104006082, 0.0005930306042269719, 0.0006542884363805873, 0.0007168254837477198, 0.0007803511351955408, 0.0008445768524095246, 0.0009092178111921936, 0.0009739944107089449, 0.0010386336410093122, 0.0011028462617825863, 0.001166421771275114, 0.0012290917899603904, 0.0012906203402795746, 0.0013507826044575483, 0.001409365392020967, 0.001466167471628324, 0.001520999776424392, 0.001573685493466051, 0.0016240600488132695, 0.0016719225751744976, 0.0017172269194679068, 0.0017597983804151784, 0.0017995195526082641, 0.0018362840754734544, 0.0018699962617597932, 0.0019005707110354708, 0.0019279319196542077, 0.0019520138978756096, 0.0019727598039535887, 0.001990049531353205, 0.002003985590472074, 0.002014466777872592, 0.002021469816218424, 0.002024979275372041, 0.0020249874579578164, 0.0020214943415505227, 0.002014507578645963, 0.0020040425544208262, 0.0019901225011133006, 0.0019726887858084186, 0.001951959326505881, 0.0019278936512818756, 0.0019005485564775326, 0.0018699899861890357, 0.001836293398168043, 0.0017995441467084095, 0.0017598378727036734, 0.0017172808901895764, 0.0016719905579205552, 0.0016239974831004922, 0.001573638283491774, 0.0015209674696502473, 0.0014661495688945148, 0.0014093613477600409, 0.001350791827067488, 0.0012906421926788259, 0.001229125590343938, 0.0011664667940978202, 0.001102901738980533, 0.0010385856569584684, 0.0009739591066666442, 0.0009091944315231526, 0.0008445646014683117, 0.000780349179380635, 0.0007168329538507212, 0.0006543044302316225, 0.0005930541896329578, 0.0005333731283773124, 0.00047555059317994086, 0.0004198086788074207, 0.000366559092646046, 0.00031600712901413974, 0.0002684158649373103, 0.00022403664343007001, 0.00018310703315738142, 0.00014584884853274804, 0.00011246625747052263, 8.314400351361328e-05, 5.8045768026038225e-05, 3.729359530280085e-05, 2.1047841402666023e-05, 9.376969552166159e-06, 2.3476087451790306e-06, -2.220446049250313e-16])
    else:
        print(f'Don\'t have data for sigma={sigma}!')
        x = None
        y = None

    return x, y


def cos_1(x, a):
    return a * (1 - np.cos(x))


def cos_1_sin(x, a, b):
    return a * (1 - np.cos(x) + b * np.sin(1.5 * x))


if __name__ == '__main__':
    main()