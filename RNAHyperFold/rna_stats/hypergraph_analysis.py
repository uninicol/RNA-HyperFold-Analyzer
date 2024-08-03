from abc import ABC, abstractmethod

import hypernetx as hnx


class StructuralHypergraphAnalysis(ABC):
    """Classe astratta per l'analisi degli ipergrafi strutturali."""

    @abstractmethod
    def secondary_structures(self) -> dict:
        """Restituisce il dizionario contenente le strutture secondarie rilevate nella temperatura selezionata.

        Returns:
            dict: Il dizionario delle strutture secondarie.
        """
        pass

    @abstractmethod
    def s_between_centrality(self, s: int = 1) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi nella temperatura selezionata.

        Args:
            s (int): Requisito di connessione.

        Returns:
            dict: La n-between-centrality dei nucleotidi.
        """
        pass

    @abstractmethod
    def connection_differences(self, hypergraph: hnx.Hypergraph) -> tuple[dict, dict]:
        """
        Restituisce le differenze di connessione nucleotide-nucleotide tra l'ipergrafo corrente e un altro ipergrafo.

        Args:
            hypergraph (hnx.Hypergraph): L'ipergrafo da mettere a confronto.

        Returns:
            tuple[dict, dict]: Due dizionari che rappresentano le connessioni aggiunte e rimosse.
        """
        pass

    @abstractmethod
    def structure_differences(self, hypergraph: hnx.Hypergraph) -> dict:
        """
        Restituisce un dizionario che indica le strutture aggiunte o rimosse dall'ipergrafo preso in input.

        Args:
            hypergraph (hnx.Hypergraph): L'ipergrafo da mettere a confronto.

        Returns:
            dict: Il dizionario delle strutture aggiunte o rimosse.
        """
        pass

    @abstractmethod
    def get_nucleotides_change_structure(self, hypergraph: hnx.Hypergraph) -> list:
        """
        Restituisce i nucleotidi che hanno subito un cambiamento di struttura.

        Args:
            hypergraph (hnx.Hypergraph): L'ipergrafo da mettere a confronto.

        Returns:
            list: La lista dei nucleotidi che hanno cambiato struttura.
        """
        pass


class CommunityHypergraphAnalysis(ABC):
    """Classe astratta per l'analisi delle comunità negli ipergrafi."""

    @abstractmethod
    def modularity(self) -> float:
        """
        Restituisce la modularità dell'ipergrafo nella temperatura selezionata.

        Returns:
            float: La modularità dell'ipergrafo.
        """
        pass

    @abstractmethod
    def subset_conductance(self, subset: set) -> float:
        """
        Restituisce la conduttanza di una partizione nella temperatura selezionata.

        Args:
            subset (set): La partizione.

        Returns:
            float: La conduttanza della partizione.
        """
        pass

    @abstractmethod
    def partitions_conductance(self) -> list[float]:
        """
        Restituisce la conduttanza di tutte le partizioni nella temperatura selezionata.

        Returns:
            list[float]: La lista contenente la conduttanza di tutte le partizioni.
        """
        pass


class TemporalRnaStats(ABC):
    """Classe astratta per l'analisi statistica temporale dell'RNA."""

    @abstractmethod
    def get_nucleotide_sensibility_to_changes(
        self, start_temp: int, end_temp: int, plot=False, plot_size: tuple = (20, 10)
    ) -> dict:
        """
        Restituisce un dizionario che indica, per ogni nucleotide, quante volte ha cambiato struttura in un range di temperature.

        Args:
            start_temp (int): La temperatura iniziale.
            end_temp (int): La temperatura finale.
            plot (bool): Indica se fare il grafico della conduttanza.
            plot_size (tuple): Se viene richiesto il grafico, definisce la sua grandezza.

        Returns:
            dict: Il dizionario delle sensibilità dei nucleotidi ai cambiamenti di temperatura.
        """
        pass

    @abstractmethod
    def get_structure_differences(self, start_temp: int, end_temp: int) -> dict:
        """
        Restituisce un dizionario contenente il numero di strutture create o rimosse in un range di temperature dalla temperatura di partenza.

        Args:
            start_temp (int): La temperatura iniziale.
            end_temp (int): La temperatura finale.

        Returns:
            dict: Il dizionario delle differenze strutturali.
        """
        pass
