from abc import ABC, abstractmethod


class RnaHypergraphStats(ABC):

    @abstractmethod
    def secondary_structures(self) -> dict:
        """
        Restituisce il dizionario contenente le strutture secondarie rilevate nella temperatura selezionata
        """
        pass

    @abstractmethod
    def modularity(self) -> float:
        """
        Restituisce la modularità dell'ipergrafo nella temperatura selezionata
        :return: la modularità dell'ipergrafo
        """
        pass

    @abstractmethod
    def subset_conductance(self, subset: set) -> float:
        """
        Restituisce la conduttanza di una partizione nella temperatura selezionata
        :param subset: la partizione
        :return: la conduttanza della partizione
        """
        pass

    @abstractmethod
    def partitions_conductance(self) -> list[float]:
        """
        Restituisce la conduttanza di tutte le partizioni nella temperatura selezionata
        :return: la lista contenente la conduttanza di tutte le partizioni
        """
        pass

    @abstractmethod
    def s_between_centrality(self, s: int = 1) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi nella temperatura selezionata
        :param s: connectedness requirement
        :return: la n-between-centrality dei nucleotidi
        """
        pass


class TemporalRnaStats(ABC):
    pass
