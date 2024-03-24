from abc import ABC, abstractmethod

import hypernetx as hnx


class StructuralHypergraphAnalysis(ABC):
    @abstractmethod
    def secondary_structures(self) -> dict:
        """
        Restituisce il dizionario contenente le strutture secondarie rilevate nella temperatura selezionata
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

    @abstractmethod
    def connection_differences(self, hypergraph: hnx.Hypergraph):
        """
        Restituisce le differenze di connessione nucleotide-nucleotide
        :param hypergraph : ipergrafo da mettere a confronto
        """
        pass

    @abstractmethod
    def structure_differences(self, hypergraph: hnx.Hypergraph) -> dict:
        """
        Restituisce un dizionario che indica le strutture aggiunte o rimosse dall'ipergrafo preso in input
        :param hypergraph : ipergrafo da mettere a confronto
        """
        pass

    @abstractmethod
    def get_nucleotides_change_structure(self, hypergraph: hnx.Hypergraph) -> list:
        """
        Restituisce i nucleotidi che hanno subito un cambiamento di struttura
        :param hypergraph : ipergrafo da mettere a confronto
        """

    pass


class CommunityHypergraphAnalysis(ABC):
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


class TemporalRnaStats(ABC):
    @abstractmethod
    def get_nucleotide_sensibility_to_changes(
        self, start_temp: int, end_temp: int, plot=False, plot_size: tuple = (20, 10)
    ) -> dict:
        """
        Restituisce un dizionario che indica, per ogni nucleotide, quante volte ha cambiato struttura in un range di temperature
        :param start_temp: la temperatura iniziale
        :param end_temp: la temperatura finale
        :param plot: indica se fare il grafico della conduttanza
        :param plot_size: se viene richiesto il grafico, definisce la sua grandezza
        """
        pass

    @abstractmethod
    def get_structure_differences(self, start_temp, end_temp):
        """
        Restituisce un dizionario contenente il numero di strutture create o rimosse in un range di temperature dalla
        temperatura di partenza
        :param start_temp: la temperatura iniziale
        :param end_temp: la temperatura finale
        """
        pass
