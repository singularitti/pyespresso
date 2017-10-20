#!/usr/bin/env python3
# created at Oct 19, 2017 11:21 PM by Qi Zhang

from numpy.linalg import inv

from output.read_file import *


class Elasticity:
    def __init__(self, filename):
        self.reo = ReadElasticityOutput()
        self.pressure_list, self.elastic_tensor_list = self.reo.read_elastic_tensor(filename)

    @staticmethod
    def create_compliance_tensor(elastic_tensor: np.ndarray) -> np.ndarray:
        """
        This method generates a compliance tensor s_ij, which is the inverse of elastic tensor.

        :return: compliance tensor
        """
        return inv(elastic_tensor)

    @staticmethod
    def bulk_modulus_voigt_average(elastic_tensor: np.ndarray) -> float:
        """
        This method gives the Voigt average of bulk modulus, which is the upper bound on K for polycrystalline material.

        :return: Voigt average of Bulk modulus
        """
        c11 = elastic_tensor[0, 0]
        c22 = elastic_tensor[1, 1]
        c33 = elastic_tensor[2, 2]
        c12 = elastic_tensor[0, 1]
        c23 = elastic_tensor[1, 2]
        c31 = elastic_tensor[2, 0]
        return 1 / 9 * (c11 + c22 + c33 + 2 * (c12 + c23 + c31))

    @staticmethod
    def bulk_modulus_reuss_average(compliance_tensor: np.ndarray) -> float:
        """
        This method gives the Reuss average of bulk modulus, which is the lower bound on K for polycrystalline material,

        :param compliance_tensor:
        :return:
        """
        s11 = compliance_tensor[0, 0]
        s22 = compliance_tensor[1, 1]
        s33 = compliance_tensor[2, 2]
        s12 = compliance_tensor[0, 1]
        s23 = compliance_tensor[1, 2]
        s31 = compliance_tensor[2, 0]
        return 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s31))

    @staticmethod
    def shear_modulus_voigt_average(elastic_tensor: np.ndarray) -> float:
        """
        This method gives the Voigt average of shear modulus, which is the upper bound on G for polycrystalline material.

        :param elastic_tensor:
        :return:
        """
        c11 = elastic_tensor[0, 0]
        c22 = elastic_tensor[1, 1]
        c33 = elastic_tensor[2, 2]
        c12 = elastic_tensor[0, 1]
        c23 = elastic_tensor[1, 2]
        c31 = elastic_tensor[2, 0]
        c44 = elastic_tensor[3, 3]
        c55 = elastic_tensor[4, 4]
        c66 = elastic_tensor[5, 5]
        return 1 / 15 * (c11 + c22 + c33 - c12 - c23 - c31 + 3 * (c44 + c55 + c66))

    @staticmethod
    def shear_modulus_reuss_average(compliance_tensor: np.ndarray) -> float:
        """
        This method gives the Reuss average of shear modulus, which is the lower bound on G for polycrystalline material.

        :param compliance_tensor:
        :return:
        """
        s11 = compliance_tensor[0, 0]
        s22 = compliance_tensor[1, 1]
        s33 = compliance_tensor[2, 2]
        s12 = compliance_tensor[0, 1]
        s23 = compliance_tensor[1, 2]
        s31 = compliance_tensor[2, 0]
        s44 = compliance_tensor[3, 3]
        s55 = compliance_tensor[4, 4]
        s66 = compliance_tensor[5, 5]
        return 15 / (4 * (s11 + s22 + s33 - 4 * (s12 + s23 + s31) + 3 * (s44 + s55 + s66)))

    def bulk_modulus_vrh_average(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """
        Voigt--Reuss--Hill average of bulk modulus.

        :param elastic_tensor:
        :param compliance_tensor:
        :return: Voigt--Reuss--Hill average of bulk modulus
        """
        kv = self.bulk_modulus_voigt_average(elastic_tensor)
        kr = self.bulk_modulus_reuss_average(compliance_tensor)
        return (kv + kr) / 2

    def shear_modulus_vrh_average(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """
        Voigt--Reuss--Hill average of shear modulus.

        :param elastic_tensor:
        :param compliance_tensor:
        :return: Voigt--Reuss--Hill average of shear modulus
        """
        gv = self.shear_modulus_voigt_average(elastic_tensor)
        gr = self.shear_modulus_reuss_average(compliance_tensor)
        return (gv + gr) / 2

    def isotropic_poisson_ratio(self, elastic_tensor: np.ndarray, compliance_tensor: np.ndarray) -> float:
        """
        A number describing lateral response to loading.

        :param elastic_tensor:
        :param compliance_tensor:
        :return:
        """
        kvrh = self.bulk_modulus_vrh_average(elastic_tensor, compliance_tensor)
        gvrh = self.shear_modulus_vrh_average(elastic_tensor, compliance_tensor)
        return (3 * kvrh - 2 * gvrh) / (6 * kvrh + 2 * gvrh)
