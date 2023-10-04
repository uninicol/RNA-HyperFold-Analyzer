from abc import ABC, abstractmethod


class RnaStats(ABC):
    @abstractmethod
    def modularity(self):
        """
        Restituisce la modularità dell'ipergrafo
        :return: la modularità dell'ipergrafo
        """
        pass

    @abstractmethod
    def get_subset_conductance(self, subset: set) -> float:
        """
        Restituisce la conduttanza di una partizione
        :param subset: la partizione
        :return: la conduttanza della partizione
        """
        pass

    @abstractmethod
    def get_partitions_conductance(self) -> enumerate[float]:
        """
        Restituisce la conduttanza di tutte le partizioni
        :return: l'enumerazione contenente la conduttanza di tutte le partizioni
        """
        pass

    @abstractmethod
    def get_n_between_centrality(self, n: int = 1) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi
        :param n: connectedness requirement
        :return: la n-between-centrality dei nucleotidi
        """
        pass
