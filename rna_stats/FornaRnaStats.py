import hypernetx as hnx
import hypernetx.algorithms.hypergraph_modularity as hmod
import matplotlib.pyplot as plt

from incidence_producers.FornaIncidenceProducer import FornaIncidenceProducer
from RnaStats import RnaStats


class FornaRnaStats(RnaStats):
    """Classe che raccoglie delle statistiche su una sequenza di RNA utilizzando un ipergrafo"""

    def __init__(self, producer: FornaIncidenceProducer) -> None:
        super().__init__(producer)
        self.__partitions: list = []
        self.__precomputed_H: list[set] = []

    def partitions(self) -> list:
        """Computa delle partizioni dell'ipergrafo"""
        if len(self.__precomputed_H) == 0:
            self.__precomputed_H = hmod.precompute_attributes(self.H)
        if len(self.__partitions) == 0:
            self.__partitions = hmod.kumar(self.__precomputed_H)
        return self.__partitions

    def partition(self, n: int) -> set:
        """
        Restituisce la partizione scelta
        :param n: numero della partizione
        :return: la partizione scelta
        """
        self.partitions()
        return self.__partitions[n]

    def modularity(self) -> float:
        """
        Restituisce la modularità dell'ipergrafo
        :return: la modularità dell'ipergrafo
        """
        self.partitions()
        return hmod.modularity(self.__precomputed_H, self.__partitions)

    def subset_conductance(self, subset: set) -> float:
        """
        Restituisce la conduttanza di una partizione
        :param subset: la partizione
        :return: la conduttanza della partizione
        """
        subset2 = [n for n in self.H.nodes if n not in subset]
        ws = sum((self.H.degree(node) for node in subset))
        was = 0
        for edge in self.H.edges:
            he_vertices = self.H.edges[edge]
            if len([n for n in he_vertices if n in subset]) == 0:
                continue
            if len([n for n in he_vertices if n in subset2]) == 0:
                continue
            was += len(he_vertices)
        return was / ws

    def partitions_conductance(self) -> list[float]:
        """
        Restituisce la conduttanza di tutte le partizioni
        :return: la lista contenente la conduttanza di tutte le partizioni
        """
        return [self.subset_conductance(subset) for subset in self.partitions()]

    def n_between_centrality(self, n: int = 1) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi
        :param n: connectedness requirement
        :return: la n-between-centrality dei nucleotidi
        """
        return hnx.algorithms.s_betweenness_centrality(self.H, n)

    def plot_hypergraph(self, size: tuple = (40, 40)) -> None:
        """Disegna un grafico che rappresenta l'ipergrafo costruito"""
        plt.subplots(figsize=size)
        hnx.draw(self.H)
        plt.show()

    def plot_partitions_conductance(self) -> None:
        """Disegna un grafico che rappresenta la conduttanza delle partizioni"""
        cond = self.partitions_conductance()
        seq = []
        values = []
        for i, c in cond:
            seq.append(i)
            values.append(c)
        plt.bar(seq, values)
        plt.title("Conductance")
        plt.xlabel("Partition")
        plt.ylabel("Conductance")
        plt.show()

    def plot_n_between_centrality(self, n: int = 1) -> None:
        """
        Disegna un grafico che rappresenta la n-between-centrality dei nucleotidi
        :param n: connectedness requirement
        """
        centrality = self.n_between_centrality(n)
        seq = list(centrality.keys())
        centr = list(centrality.values())

        plt.bar(seq, centr)
        plt.title(f"{n}-centrality")
        plt.xlabel("Nucleotides")
        plt.ylabel("Centrality")
        plt.show()
