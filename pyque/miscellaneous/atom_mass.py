#!/usr/bin/env python3
# created by Michel, changed on Aug 9 at 2017 12:46 PM by Qi Zhang
"""
:mod:`atom_mass` -- Atomic mass for elements
============================================

.. module:: atom_mass
   :platform: Unix, Windows, Mac, Linux
   :synopsis: 8-digits atomic mass for each element up to Bi.
.. moduleauthor:: Michel Marcondes <mld2189@columbia.edu>
.. moduleauthor:: Qi Zhang <qz2280@columbia.edu>
"""

# Dictionary of atomic mass
atom_mass = {'H': 1.00800000,
             'He': 4.00260200,
             'Li': 6.94000000,
             'Be': 9.01218310,
             'B': 10.8100000,
             'C': 12.0110000,
             'N': 14.0070000,
             'O': 15.9990000,
             'F': 18.9984032,
             'Ne': 20.1797000,
             'Na': 22.9897693,
             'Mg': 24.3050000,
             'Al': 26.9815385,
             'Si': 28.0850000,
             'P': 30.9737620,
             'S': 32.0600000,
             'Cl': 35.4500000,
             'Ar': 39.9480000,
             'K': 39.0983000,
             'Ca': 40.0780000,
             'Sc': 44.9559080,
             'Ti': 47.8670000,
             'V': 50.9415000,
             'Cr': 51.9961000,
             'Mn': 54.9380440,
             'Fe': 55.8450000,
             'Co': 58.9331940,
             'Ni': 58.6934000,
             'Cu': 63.5460000,
             'Zn': 65.3800000,
             'Ga': 69.7230000,
             'Ge': 72.6300000,
             'As': 74.9215950,
             'Se': 78.9710000,
             'Br': 79.9040000,
             'Kr': 83.7980000,
             'Rb': 85.4678000,
             'Sr': 87.6200000,
             'Y': 88.9058400,
             'Zr': 91.2240000,
             'Nb': 92.9063700,
             'Mo': 95.9500000,
             'Tc': 97.9100000,
             'Ru': 101.0700000,
             'Rh': 102.9055000,
             'Pd': 106.4200000,
             'Ag': 107.8682000,
             'Cd': 112.4140000,
             'In': 114.8180000,
             'Sn': 118.7100000,
             'Sb': 121.7600000,
             'Te': 127.6000000,
             'I': 126.9044700,
             'Xe': 131.2930000,
             'Cs': 132.9054520,
             'Ba': 137.3270000,
             'La': 138.9054700,
             'Ce': 140.1160000,
             'Pr': 140.9076600,
             'Nd': 144.2420000,
             'Pm': 144.9100000,
             'Sm': 150.3600000,
             'Eu': 151.9640000,
             'Gd': 157.2500000,
             'Tb': 158.9253500,
             'Dy': 162.5000000,
             'Ho': 164.9303300,
             'Er': 167.2590000,
             'Tm': 168.9342200,
             'Yb': 173.0450000,
             'Lu': 174.9668000,
             'Hf': 178.4900000,
             'Ta': 180.9478800,
             'W': 183.8400000,
             'Re': 186.2070000,
             'Os': 190.2300000,
             'Ir': 192.2170000,
             'Pt': 195.0840000,
             'Au': 196.9665690,
             'Hg': 200.5920000,
             'Tl': 204.3800000,
             'Pb': 207.2000000,
             'Bi': 208.9804000}


def get_atom_mass(name: str) -> float:
    """
    Give the name of an element, and its atomic mass is returned.

    :param name: The name of an element.
    :return: The atomic mass of an element.
    """
    return atom_mass[name]
