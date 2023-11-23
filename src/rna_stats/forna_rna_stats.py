import hypernetx as hnx
import hypernetx.algorithms.hypergraph_modularity as hmod
import matplotlib.pyplot as plt

from src.rna_stats.rna_stats import RnaHypergraphStats


class FornaRnaStats(RnaHypergraphStats):
    """Classe che raccoglie delle statistiche su una sequenza di RNA utilizzando un ipergrafo"""

    def __init__(self, HG: hnx.Hypergraph) -> None:
        self.HG = HG
        self.__partitions: list = []

    def secondary_structures(self) -> dict:
        structures = {}
        for key, value in self.HG.incidence_dict.items():
            if not key.startswith("l") and not key.startswith("db"):
                structures[key] = value
        return structures

    def partitions(self) -> list:
        """Computa delle partizioni dell'ipergrafo"""
        if len(self.__partitions) == 0:
            self.__partitions = hmod.kumar(self.HG)
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
        return hmod.modularity(self.HG, self.__partitions)

    def subset_conductance(self, subset: set) -> float:
        """
        Restituisce la conduttanza di una partizione
        :param subset: la partizione
        :return: la conduttanza della partizione
        """
        return hmod.conductance(self.HG, subset)

    def partitions_conductance(self) -> list[float]:
        """
        Restituisce la conduttanza di tutte le partizioni
        :return: la lista contenente la conduttanza di tutte le partizioni
        """
        return [self.subset_conductance(subset) for subset in self.partitions()]

    def s_between_centrality(self, s=1) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi
        :param s: connectedness requirement
        :return: la n-between-centrality dei nucleotidi
        """
        return hnx.algorithms.s_betweenness_centrality(self.HG, s)

    def plot_hypergraph(self, size: tuple = (40, 40)) -> None:
        """Disegna un grafico che rappresenta l'ipergrafo costruito"""
        plt.subplots(figsize=size)
        hnx.draw(self.HG, **{'layout_kwargs': {'seed': 39}})
        plt.show()

    def plot_partitions_conductance(self) -> None:
        """Disegna un grafico che rappresenta la conduttanza delle partizioni"""
        cond = self.partitions_conductance()
        seq = []
        values = []
        for i, c in enumerate(cond):
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
        centrality = self.s_between_centrality()
        seq = list(centrality.keys())
        centr = list(centrality.values())

        plt.bar(seq, centr)
        plt.title(f"{n}-centrality")
        plt.xlabel("Nucleotides")
        plt.ylabel("Centrality")
        plt.show()
