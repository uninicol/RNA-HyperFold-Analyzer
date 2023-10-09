from abc import ABC, abstractmethod
import hypernetx as hnx

import IncidenceProducer


class RnaStats(ABC):
    def __init__(self, producer: IncidenceProducer) -> None:
        incidence_dict = producer.get_incidence_dict()
        self.H = hnx.Hypergraph(incidence_dict)
        del incidence_dict
        pass

    @abstractmethod
    def modularity(self):
        """
        Restituisce la modularità dell'ipergrafo
        :return: la modularità dell'ipergrafo
        """
        pass

    @abstractmethod
    def subset_conductance(self, subset: set) -> float:
        """
        Restituisce la conduttanza di una partizione
        :param subset: la partizione
        :return: la conduttanza della partizione
        """
        pass

    @abstractmethod
    def partitions_conductance(self) -> list[float]:
        """
        Restituisce la conduttanza di tutte le partizioni
        :return: la lista contenente la conduttanza di tutte le partizioni
        """
        pass

    @abstractmethod
    def n_between_centrality(self, n: int = 1) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi
        :param n: connectedness requirement
        :return: la n-between-centrality dei nucleotidi
        """
        pass
