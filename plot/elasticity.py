#!/usr/bin/env python3
# created at Oct 20, 2017 12:49 AM by Qi Zhang

from compute.elasticity import *
from plot.plot_basic import *


class ElasticityOutputPlotter(SingleAxes):
    def __init__(self, file: str):
        super().__init__()
        self.reo = ElasticityOutputReader(file)
        self.ec = ElasticityCalculator(file)

    def plot_cij_vs_pressures(self):
        data = self.reo.read_cij_vs_pressures()
        for i, p in enumerate(zip(data.keys(), data.values())):
            for pp in zip([float(p[0])] * len(p[1]), p[1]):
                print(*pp)
                plt.scatter(*pp, label="value" + str(i))

    def _plot_x_vs_pressure(self, func: Callable, items):
        """
        This is a generic method that can be used to implement the following methods.

        :param func: a function takes one or more arguments
        :param items: can be a list or zip object
        :return: a list of lines that were added to plot
        """
        try:
            for item in items:
                if isinstance(item, tuple):
                    # If the `func` function needs more than one arguments,
                    # an item will be a tuple and need to be distinguished.
                    return plt.plot(self.ec.pressures, [func(*item) for item in items])
                else:
                    # If the `func` function only needs one argument.
                    return plt.plot(self.ec.pressures, [func(item) for item in items])
        except ValueError:
            # Rewrite error information
            print('Pressure and ' + what_to_be_plot + ' do not have the same length! Python will exit!')
            exit()

    def plot_bulk_modulus_voigt_average_vs_pressure(self):
        """
        This method plots the Voigt average calculated by method `derive_bulk_modulus_voigt_average` versus
        corresponding pressure.

        :return:
        """
        what_to_be_plot = 'Voigt average of bulk modulus'
        return self._plot_x_vs_pressure(self.ec.derive_bulk_modulus_voigt_average, self.ec.elastic_tensors)

    def plot_bulk_modulus_reuss_average_vs_pressure(self):
        """

        :return:
        """
        return self._plot_x_vs_pressure(self.ec.derive_bulk_modulus_reuss_average, self.ec.compliance_tensors)

    def plot_shear_modulus_voigt_average_vs_pressure(self):
        return self._plot_x_vs_pressure(self.ec.derive_shear_modulus_voigt_average, self.ec.elastic_tensors)

    def plot_shear_modulus_reuss_average_vs_pressure(self):
        return self._plot_x_vs_pressure(self.ec.derive_shear_modulus_reuss_average, self.ec.compliance_tensors)

    def plot_bulk_modulus_vrh_average_vs_pressure(self):
        print(len(self.ec.pressures),
              len(list(zip(self.ec.elastic_tensors, self.ec.compliance_tensors))))
        return self._plot_x_vs_pressure(self.ec.derive_bulk_modulus_vrh_average,
                                        list(zip(self.ec.elastic_tensors, self.ec.compliance_tensors)))

    def plot_shear_modulus_vrh_average_vs_pressure(self):
        return self._plot_x_vs_pressure(self.ec.derive_shear_modulus_vrh_average,
                                        list(zip(self.ec.elastic_tensors, self.ec.compliance_tensors)))

    def plot_isotropic_poisson_ratio_vs_pressure(self):
        return self._plot_x_vs_pressure(self.ec.derive_isotropic_poisson_ratio,
                                        zip(self.ec.elastic_tensors, self.ec.compliance_tensors))

    def plot_universal_elastic_anisotropy_vs_pressure(self):
        return self._plot_x_vs_pressure(self.ec.derive_universal_elastic_anisotropy,
                                        zip(self.ec.elastic_tensors, self.ec.compliance_tensors))
