#!/usr/bin/env python3
# created at Oct 20, 2017 12:49 AM by Qi Zhang

from matplotlib.lines import Line2D
from pyespresso.calculator.elasticity import *
from pyespresso.plotting.plot_basic import *


class ElasticityOutputPlotter(SingleAxes):
    def __init__(self, file: str):
        super().__init__()
        self._reo = ElasticityOutputLexer(file)
        self._ec = ElasticityCalculator(file)

    def get_cijs(self, crystal_class):
        cijs = np.empty(5)
        for i, p in enumerate(self._ec.pressures):
            cijs = np.vstack(
                (cijs, self._ec.grab_independent_tensor_values(self._ec.elastic_tensors[i], crystal_class)))
        cijs = cijs.transpose()
        return cijs

    def plot_cij_vs_pressures(self, crystal_class):
        cijs = self.get_cijs(crystal_class)
        for i in range(5):
            plt.plot(self._ec.pressures, cijs[i][1:], label=crystal_classes[crystal_class]['names'][i])

    def plot_cij_vs_volumes(self, crystal_class, volumes):
        cijs = self.get_cijs(crystal_class)
        for i in range(5):
            plt.plot(volumes, cijs[i][1:], label=crystal_classes[crystal_class]['names'][i])

    def plot_bulk_modulus_voigt_average_vs_pressure(self) -> Line2D:
        """
        This method plots the Voigt average calculated by method `derive_bulk_modulus_voigt_average` versus
        corresponding pressure.

        :return: a line that were added to plotters
        """
        return self._plot_f_vs_pressure(self._ec.derive_bulk_modulus_voigt_average, self._ec.elastic_tensors)

    def plot_bulk_modulus_reuss_average_vs_pressure(self) -> Line2D:
        """
        This method plots the Reuss average calculated by method `derive_bulk_modulus_reuss_average` versus
        corresponding pressure.

        :return: a line that were added to plotters
        """
        return self._plot_f_vs_pressure(self._ec.derive_bulk_modulus_reuss_average, self._ec.compliance_tensors)

    def plot_shear_modulus_voigt_average_vs_pressure(self) -> Line2D:
        return self._plot_f_vs_pressure(self._ec.derive_shear_modulus_voigt_average, self._ec.elastic_tensors)

    def plot_shear_modulus_reuss_average_vs_pressure(self) -> Line2D:
        return self._plot_f_vs_pressure(self._ec.derive_shear_modulus_reuss_average, self._ec.compliance_tensors)

    def plot_bulk_modulus_vrh_average_vs_pressure(self) -> Line2D:
        return self._plot_f_vs_pressure(self._ec.derive_bulk_modulus_vrh_average,
                                        zip(self._ec.elastic_tensors, self._ec.compliance_tensors))

    def plot_shear_modulus_vrh_average_vs_pressure(self) -> Line2D:
        return self._plot_f_vs_pressure(self._ec.derive_shear_modulus_vrh_average,
                                        zip(self._ec.elastic_tensors, self._ec.compliance_tensors))

    def plot_isotropic_poisson_ratio_vs_pressure(self) -> Line2D:
        return self._plot_f_vs_pressure(self._ec.derive_isotropic_poisson_ratio,
                                        zip(self._ec.elastic_tensors, self._ec.compliance_tensors))

    def plot_universal_elastic_anisotropy_vs_pressure(self) -> Line2D:
        return self._plot_f_vs_pressure(self._ec.derive_universal_elastic_anisotropy,
                                        zip(self._ec.elastic_tensors, self._ec.compliance_tensors))

    def _plot_f_vs_pressure(self, f: Callable,
                            items: Union[List[np.ndarray], Iterator[Tuple[np.ndarray, np.ndarray]]]) -> Line2D:
        """
        This is a generic method that can be used to implement the following methods.

        :param f: a function takes one or more arguments.
        :param items: can be a list or zip object
        :return: a line that were added to plotters
        """
        # If `f` takes more than one arguments, an item will be a tuple, and need to be distinguished from the case
        # where it only takes one argument, where an item is just a numpy array.
        # This is because "tuple parameter unpacking" is removed from Python 3.x.
        line, = plt.plot(self._ec.pressures, [f(*item) if isinstance(item, tuple) else f(item) for item in items])
        return line
