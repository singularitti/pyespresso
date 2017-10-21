#!/usr/bin/env python3
# created at Oct 19, 2017 11:21 PM by Qi Zhang

from numpy.linalg import inv

from read_file.elasticity import *


class ElasticityCalculator:
    def __init__(self, file):
        """
        Usually we know elastic tensors, which are calculated by Quantum ESPRESSO, add we need to calculate
        corresponding compliance tensors as well as other important values manually.

        :param file: a file giving an 6x6 elastic tensor under each pressure.
        """
        eor = ElasticityOutputReader(file)
        self.pressures, self.elastic_tensors = eor.read_elastic_tensor()
        self.compliance_tensors: List[np.ndarray] = [self.create_compliance_tensor(e) for e in self.elastic_tensors]

    @staticmethod
    def create_compliance_tensor(elastic_tensor: np.ndarray) -> np.ndarray:
        """
        This method generates a compliance tensor s_ij, which is the inverse of elastic tensor.

        :return: compliance tensor
        """
        return inv(elastic_tensor)

    @staticmethod
    def derive_bulk_modulus_voigt_average(elastic_tensor: np.ndarray) -> float:
        """
        This method gives the Voigt average of bulk modulus, which is the upper bound on K for polycrystalline material.

        :return: Voigt average of Bulk modulus
        """
        c11, c22, c33, c12, c23, c31 = _get_values_by_indices(elastic_tensor,
                                                              [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)])
        return 1 / 9 * (c11 + c22 + c33 + 2 * (c12 + c23 + c31))

    @staticmethod
    def derive_bulk_modulus_reuss_average(compliance_tensor: np.ndarray) -> float:
        """
        This method gives the Reuss average of bulk modulus, which is the lower bound on K for polycrystalline material,

        :param compliance_tensor:
        :return: Reuss average of bulk modulus
        """
        s11, s22, s33, s12, s23, s31 = _get_values_by_indices(compliance_tensor,
                                                              [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)])
        return 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s31))

    @staticmethod
    def derive_shear_modulus_voigt_average(elastic_tensor: np.ndarray) -> float:
        """
        This method gives the Voigt average of shear modulus, which is the upper bound on G for polycrystalline material.

        :param elastic_tensor:
        :return: Voigt average of shear modulus
        """
        c11, c22, c33, c12, c23, c31, c44, c55, c66 = _get_values_by_indices(elastic_tensor,
                                                                             [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2),
                                                                              (2, 0), (3, 3), (4, 4), (5, 5)])
        return 1 / 15 * (c11 + c22 + c33 - c12 - c23 - c31 + 3 * (c44 + c55 + c66))

    @staticmethod
    def derive_shear_modulus_reuss_average(compliance_tensor: np.ndarray) -> float:
        """
        This method gives the Reuss average of shear modulus, which is the lower bound on G for polycrystalline material.

        :param compliance_tensor:
        :return: Reuss average of shear modulus
        """
        s11, s22, s33, s12, s23, s31, s44, s55, s66 = _get_values_by_indices(compliance_tensor,
                                                                             [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2),
                                                                              (2, 0), (3, 3), (4, 4), (5, 5)])
        return 15 / (4 * (s11 + s22 + s33 - 4 * (s12 + s23 + s31) + 3 * (s44 + s55 + s66)))

    def derive_bulk_modulus_vrh_average(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """
        Voigt--Reuss--Hill average of bulk modulus.

        :param elastic_tensor:
        :param compliance_tensor:
        :return: Voigt--Reuss--Hill average of bulk modulus
        """
        kv = self.derive_bulk_modulus_voigt_average(elastic_tensor)
        kr = self.derive_bulk_modulus_reuss_average(compliance_tensor)
        return (kv + kr) / 2

    def derive_shear_modulus_vrh_average(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """
        Voigt--Reuss--Hill average of shear modulus.

        :param elastic_tensor:
        :param compliance_tensor:
        :return: Voigt--Reuss--Hill average of shear modulus
        """
        gv = self.derive_shear_modulus_voigt_average(elastic_tensor)
        gr = self.derive_shear_modulus_reuss_average(compliance_tensor)
        return (gv + gr) / 2

    def derive_isotropic_poisson_ratio(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """
        A number describing lateral response to loading.

        :param elastic_tensor:
        :param compliance_tensor:
        :return: isotropic Poisson ratio
        """
        kvrh = self.derive_bulk_modulus_vrh_average(elastic_tensor, compliance_tensor)
        gvrh = self.derive_shear_modulus_vrh_average(elastic_tensor, compliance_tensor)
        return (3 * kvrh - 2 * gvrh) / (6 * kvrh + 2 * gvrh)

    def derive_universal_elastic_anisotropy(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """

        :param elastic_tensor:
        :param compliance_tensor:
        :return:
        """
        kv = self.derive_bulk_modulus_voigt_average(elastic_tensor)
        kr = self.derive_bulk_modulus_reuss_average(compliance_tensor)
        gv = self.derive_shear_modulus_voigt_average(elastic_tensor)
        gr = self.derive_shear_modulus_reuss_average(compliance_tensor)
        return 5 * (gv / gr) + (kv / kr) - 6


def _get_values_by_indices(matrix: np.ndarray, indices: List[Tuple[int, int]]) -> List:
    """
    Given a 2-dimensional matrix $m$, and a list of indices $(i, j)$ where $i$, $j$ denote the index for the
    0th and 1st axis, respectively. Then return a list of values $m(i, j)$.

    :param matrix: a numpy array of floats, integers, etc.
    :param indices: a list of 2-tuple-integers
    :return: a list of values corresponding to those indices.
    """
    return [matrix[index] for index in indices]
