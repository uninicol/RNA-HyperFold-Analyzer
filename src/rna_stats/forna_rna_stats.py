from collections import defaultdict

import hypernetx as hnx
import hypernetx.algorithms.hypergraph_modularity as hmod
import matplotlib.pyplot as plt
import networkx as nx

from hypergraph_folding.temperature_hypergraph import TemperatureFoldingHypergraph
from rna_stats.rna_stats import RnaHypergraphStats, TemporalRnaStats


class RnaStats(RnaHypergraphStats):
    """Classe che raccoglie delle statistiche su una sequenza di RNA utilizzando un ipergrafo"""

    def __init__(self, HG: hnx.Hypergraph) -> None:
        self.HG = HG
        self.__partitions: list = []
        self.__plotter = RnaStatsPlotter()

    def plot_hypergraph(self, size: tuple = (40, 40)) -> None:
        """Disegna un grafico che rappresenta l'ipergrafo costruito"""
        self.__plotter.plot_hypergraph(self.HG, size)

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

    def partitions_conductance(self, plot=False, plot_size=(20, 10)) -> list[float]:
        """
        Restituisce la conduttanza di tutte le partizioni
        :return: la lista contenente la conduttanza di tutte le partizioni
        """
        conductances = [self.subset_conductance(subset) for subset in self.partitions()]
        if plot:
            self.__plotter.plot_partitions_conductance(conductances, size=plot_size)
        return conductances

    def n_between_centrality(self, n=1, plot=False, plot_size=(20, 10)) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi
        :param n: connectedness requirement
        :return: la n-between-centrality dei nucleotidi
        """
        centrality = hnx.algorithms.s_betweenness_centrality(self.HG, n)
        if plot:
            self.__plotter.plot_n_between_centrality(centrality, n=n, size=plot_size)
        return centrality

    def connection_differences(self, hypergraph: hnx.Hypergraph):
        """
        Restituisce
        :param hypergraph : ipergrafo da mettere a confronto
        """
        if len(self.HG.nodes) != len(hypergraph.nodes):
            raise Exception("Ipergrafi hanno un numero diverso di nodi")
        this_connections = [v for k, v in self.HG.incidence_dict.items() if k[0:2] == 'db']
        other_connections = [v for k, v in hypergraph.incidence_dict.items() if k[0:2] == 'db']
        return [item for item in this_connections if item not in other_connections]

    def structure_differences(self, hypergraph: hnx.Hypergraph) -> dict:
        """
        Restituisce un dizionario che indica le strutture aggiunte o rimosse dall'ipergrafo preso in input
        :param hypergraph : ipergrafo da mettere a confronto
        """
        if self.HG is hypergraph:
            return {}
        if len(self.HG.nodes) != len(hypergraph.nodes):
            raise Exception("Ipergrafi hanno un numero diverso di nodi")
        this_structures = self.secondary_structures()
        other_structures = RnaStats(hypergraph).secondary_structures()

        this_count = defaultdict(int)
        for name in this_structures.keys():
            this_count[name[0]] = this_count[name[0]] + 1
        other_count = defaultdict(int)
        for name in other_structures.keys():
            other_count[name[0]] = other_count[name[0]] + 1

        result = {key: this_count[key] - other_count[key] for key in this_count if
                  key in other_count and this_count[key] != other_count[key]}
        return result

    def get_nucleotides_change_structure(self, hypergraph: hnx.Hypergraph) -> list:
        """
        Restituisce i nucleotidi che hanno subito un cambiamento di struttura
        :param hypergraph : ipergrafo da mettere a confronto
        """
        if len(self.HG.nodes) != len(hypergraph.nodes):
            raise Exception("Ipergrafi non hanno lo stesso numero di nucleotidi")
        this_structures = self.secondary_structures()
        other_structures = RnaStats(hypergraph).secondary_structures()
        differences = []
        # Scorro le strutture dei due ipergrafi e individuo i nucleotidi che cambiano struttura
        for this_name, this_structure in this_structures.items():
            for other_name, other_structure in other_structures.items():
                if this_name == other_name:
                    differences.append([item for item in this_structure if item not in other_structure])
                    break

        return [element for row in differences for element in row]


class TemperatureFoldingStats(TemporalRnaStats):

    def __init__(self, THG: TemperatureFoldingHypergraph) -> None:
        self.THG = THG
        self.__plotter = TemperatureFoldingStatsPlotter()

    def get_nucleotide_sensibility_to_changes(self, start_temp: int, end_temp: int, plot=False, plot_size=(20, 10)):
        self.THG.insert_temperature_range(start_temp, end_temp)
        counts = defaultdict(int)
        h1 = self.THG.get_hypergraph(start_temp)
        for temp in range(start_temp + 1, end_temp + 1):
            h2 = self.THG.get_hypergraph(temp)
            if h1 is h2:
                continue
            st = RnaStats(h1)
            elements = st.get_nucleotides_change_structure(h2)
            for elem in elements:
                counts[elem] += 1
        if plot:
            self.__plotter.plot_nucleotide_sensibility_to_changes(counts, start_temp, end_temp, size=plot_size)
        return counts

    def get_structure_differences(self, start_temp, end_temp):
        self.THG.insert_temperature_range(start_temp, end_temp)
        diffs = {}
        h1 = self.THG.get_hypergraph(start_temp)
        for temp in range(start_temp + 1, end_temp + 1):
            h2 = self.THG.get_hypergraph(temp)
            st = RnaStats(h1)
            elements = st.structure_differences(h2)
            diffs[temp] = elements
        return diffs


class RnaStatsPlotter:

    def plot_hypergraph(self, HG: hnx.Hypergraph, size) -> None:
        """Disegna un grafico che rappresenta l'ipergrafo costruito"""
        if len(HG.nodes) > 250:
            plt.subplots(figsize=size)
            G = hmod.two_section(HG).to_networkx()
            nx.draw(G, node_size=50)
            plt.show()
        else:
            plt.subplots(figsize=size)
            hnx.draw(HG, **{'layout_kwargs': {'seed': 39}})
            plt.show()

    def plot_partitions_conductance(self, conductances, size=(20, 10)) -> None:
        """Disegna un grafico che rappresenta la conduttanza delle partizioni"""
        plt.subplots(figsize=size)
        seq = []
        values = []
        for i, c in enumerate(conductances):
            seq.append(i)
            values.append(c)
        plt.bar(seq, values)
        plt.title("Conductance")
        plt.xlabel("Partition")
        plt.ylabel("Conductance")
        plt.show()

    def plot_n_between_centrality(self, centrality, n: int = 1, size=(20, 10)) -> None:
        """
        Disegna un grafico che rappresenta la n-between-centrality dei nucleotidi
        :param n: connectedness requirement
        """
        plt.subplots(figsize=size)
        seq = list(centrality.keys())
        centr = list(centrality.values())

        plt.bar(seq, centr)
        plt.title(f"{n}-centrality")
        plt.xlabel("Nucleotides")
        plt.ylabel("Centrality")
        plt.show()


class TemperatureFoldingStatsPlotter:

    def plot_nucleotide_sensibility_to_changes(self, sensibilities, start_temp, end_temp, size=(20, 10)):
        ordered_keys = sorted(sensibilities.keys())
        ordered_values = [sensibilities[key] for key in ordered_keys]
        seq = list(ordered_keys)
        centr = list(ordered_values)
        plt.subplots(figsize=size)
        plt.bar(seq, centr)
        plt.show()
