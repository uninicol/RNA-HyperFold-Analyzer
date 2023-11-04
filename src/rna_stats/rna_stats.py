from abc import ABC, abstractmethod

import hypernetx as hnx

from src.incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


class RnaStats(ABC):
    def __init__(self, producer: TemperatureIncidenceProducer, temperature: int) -> None:
        incidence_dict = producer.get_temperature_incidence_dict(temperature)
        self.H = hnx.Hypergraph(incidence_dict)
        del incidence_dict
        pass

    @abstractmethod
    def secondary_structures(self) -> dict:
        """
        Restituisce il dizionario contenente le strutture secondarie rilevate
        """
        pass

    @abstractmethod
    def modularity(self) -> float:
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
