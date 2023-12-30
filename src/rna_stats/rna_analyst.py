from collections import defaultdict
from statistics import mean

import hypernetx as hnx
import hypernetx.algorithms.hypergraph_modularity as hmod
import matplotlib.pyplot as plt
import networkx as nx

from hypergraph_folding.temperature_hypergraph import TemperatureFoldingHypergraph
from rna_stats.hypergraph_analysis import StructuralHypergraphAnalysis, CommunityHypergraphAnalysis, TemporalRnaStats


class RnaAnalyst(StructuralHypergraphAnalysis, CommunityHypergraphAnalysis):
    """Classe che raccoglie metodi di analisi per una sequenza di RNA utilizzando un ipergrafo"""

    def __init__(self, HG: hnx.Hypergraph) -> None:
        if HG is None:
            raise Exception("None not valid")
        self.HG = HG
        self.__partitions: list = []
        self.__plotter = RnaStatsPlotter()

    def plot_hypergraph(self, size: tuple = (40, 40)) -> None:
        """Disegna un grafico che rappresenta l'ipergrafo costruito"""
        self.__plotter.plot_hypergraph(self.HG, size)

    def plot_structures(self):
        self.__plotter.plot_structures(self.secondary_structures())

    def secondary_structures(self) -> dict:
        """Restituisce le strutture secondarie rilevate"""
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
        :param plot: indica se fare il grafico della conduttanza
        :param plot_size: se viene richiesto il grafico, definisce la sua grandezza
        :return: la lista contenente la conduttanza di tutte le partizioni
        """
        conductances = [self.subset_conductance(subset) for subset in self.partitions()]
        if plot:
            self.__plotter.plot_partitions_conductance(conductances, size=plot_size)
        return conductances

    def s_between_centrality(self, s=1, edges=False, plot=False, plot_size=(20, 10)) -> dict:
        """
        Restituisce la n-between-centrality dei nucleotidi
        :param s: connectedness requirement
        :param plot: indica se fare stampare il grafico
        :param plot_size: se viene richiesto il grafico, definisce la sua grandezza
        :return: la n-between-centrality dei nucleotidi
        """
        centrality = hnx.algorithms.s_betweenness_centrality(self.HG, s=s, edges=edges)
        if plot:
            self.__plotter.plot_s_between_centrality(centrality,edges,   s=s, size=plot_size)
        return centrality

    def connection_differences(self, hypergraph: hnx.Hypergraph):
        """
        Restituisce le differenze di connessione nucleotide-nucleotide
        :param hypergraph : ipergrafo da mettere a confronto
        """
        if self.HG is hypergraph:
            return [],[]
        if len(self.HG.nodes) != len(hypergraph.nodes):
            raise Exception("Ipergrafi hanno un numero diverso di nodi")
        this_connections = [tuple(v) for k, v in self.HG.incidence_dict.items() if k[0:2] == 'db']
        other_connections = [tuple(v) for k, v in hypergraph.incidence_dict.items() if k[0:2] == 'db']

        t = set(this_connections)
        o = set(other_connections)
        new_connections = o - t
        #new_connections =  [item for other in other_connections for item in this_connections if item not in other_connections]
        if len(new_connections) == 0:
            return [],[]
        old_connections = set()
        for this in this_connections:
            for new in new_connections:
                if new[0] == this[0]:
                    old_connections.add(this)
        return new_connections, old_connections

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
        other_structures = RnaAnalyst(hypergraph).secondary_structures()

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
        other_structures = RnaAnalyst(hypergraph).secondary_structures()
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

    def get_nucleotide_sensibility_to_changes(self, start_temp: int, end_temp: int, plot=False,
                                              plot_size: tuple = (20, 10)) -> dict:
        """
        Restituisce un dizionario che indica, per ogni nucleotide, la sua sensibilità a cambiare struttura in un range di temperature
        :param start_temp: la temperatura iniziale
        :param end_temp: la temperatura finale
        :param plot: indica se fare il grafico della conduttanza
        :param plot_size: se viene richiesto il grafico, definisce la sua grandezza
        """
        self.THG.insert_temperature_range(start_temp, end_temp)
        counts = defaultdict(int)
        for temp in range(start_temp + 1, end_temp ):
            h1 = self.THG.get_hypergraph(temp)
            h2 = self.THG.get_hypergraph(temp +1 )
            if h1 is h2:
                continue
            st = RnaAnalyst(h1)
            elements = st.get_nucleotides_change_structure(h2)
            for elem in elements:
                counts[elem] += 1
        #for k, v in counts.items():
         #   counts[k] = v / (end_temp - start_temp)
        if plot:
            self.__plotter.plot_nucleotide_sensibility_to_changes(counts, size=plot_size)
        # TODO fare in modo di vedere effettivamente struttura di partenza e di arriv
        return counts

    def get_structure_differences(self, start_temp, end_temp, plot=False):
        """
        Restituisce un dizionario contenente il numero di strutture create o rimosse in un range di temperature dalla
        temperatura di partenza
        :param start_temp: la temperatura iniziale
        :param end_temp: la temperatura finale
        """
        self.THG.insert_temperature_range(start_temp, end_temp)
        diffs = {}
        h1 = self.THG.get_hypergraph(start_temp)
        for temp in range(start_temp + 1, end_temp + 1):
            h2 = self.THG.get_hypergraph(temp)
            st = RnaAnalyst(h1)
            elements = st.structure_differences(h2)
            diffs[temp] = elements

        if plot:
            self.__plotter.plot_structure_differences(diffs)
        return diffs

    def get_connection_differences(self, start_temp, end_temp, plot = False):
        self.THG.insert_temperature_range(start_temp, end_temp)
        diffs = {}
        h1 = self.THG.get_hypergraph(start_temp)
        st = RnaAnalyst(h1)
        for temp in range(start_temp + 1, end_temp + 1):
            h2 = self.THG.get_hypergraph(temp)
            elements = st.connection_differences(h2)
            diffs[temp] = elements

        if plot:
            self.__plotter.plot_connection_differences(diffs)
        return diffs


class RnaStatsPlotter:

    def plot_hypergraph(self, HG: hnx.Hypergraph, size) -> None:
        """
        Disegna un grafico che rappresenta l'ipergrafo costruito
        :param HG: l'ipergrafo da disegnare
        """
        if len(HG.nodes) > 250:
            plt.subplots(figsize=size)
            G = hmod.two_section(HG).to_networkx()
            nx.draw(G, node_size=50)
            plt.show()
        else:
            plt.subplots(figsize=size)
            hnx.draw(HG, **{'layout_kwargs': {'seed': 39}})
            plt.show()

    def plot_structures(self, structures):
        HG = hnx.Hypergraph(structures)
        hnx.draw(HG)

    def plot_partitions_conductance(self, conductances, size=(20, 10)) -> None:
        """Disegna un grafico che rappresenta la conduttanza delle partizioni"""
        plt.subplots(figsize=size)
        seq = []
        values = []
        for i, c in enumerate(conductances):
            seq.append(i)
            values.append(c)
        plt.bar(seq, values)
        plt.title("Partitions Conductance")
        plt.xlabel("Partition")
        plt.ylabel("Conductance score")
        plt.show()

    def plot_s_between_centrality(self, centrality, edges, s: int = 1, size=(20, 10)) -> None:
        """
        Disegna un grafico che rappresenta la n-between-centrality dei nucleotidi
        :param s: connectedness requirement
        """
        plt.subplots(figsize=size)
        if edges:
            centrality = {k: v for k, v in sorted(centrality.items())}
            plt.xticks(rotation=90)

        seq = list(centrality.keys())
        centr = list(centrality.values())

        plt.bar(seq, centr)
        plt.title(f"{s}-centrality")
        plt.xlabel("Nucleotides")
        plt.ylabel("Centrality")
        plt.show()


class TemperatureFoldingStatsPlotter:

    def plot_nucleotide_sensibility_to_changes(self, sensibilities, size=(20, 10)):
        ordered_keys = sorted(sensibilities.keys())
        ordered_values = [sensibilities[key] for key in ordered_keys]
        seq = list(ordered_keys)
        centr = list(ordered_values)
        plt.subplots(figsize=size)
        plt.bar(seq, centr)
        plt.title("Nucleotide Sensibility")
        plt.xlabel("Nucleotides")
        plt.ylabel("Sensibility score")
        plt.show()

    def plot_structure_differences(self, diffs):
        positives = defaultdict(list)
        negatives = defaultdict(list)
        for dic in diffs.values():
            for k, v in dic.items():
                if v > 0:
                    positives[k].append(v)
                else:
                    negatives[k].append(v)

        structures = sorted(positives.keys())
        positives = {k: (mean(positives[k]) if len(positives[k])>0 else 0) for k in structures }.values()
        negatives = {k: (mean(negatives[k]) if len(negatives[k])>0 else 0) for k in structures }.values()
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.bar(structures, negatives, color='r')
        ax.bar(structures, positives, color='b', bottom = 0)
        plt.title("Structure Differences")
        plt.xlabel("Structure")
        plt.ylabel("Structures added/removed")
        fig.show()

    def plot_connection_differences(self, diffs):
        def draw_arrow(plt, arr_start, arr_end):
            dx = arr_end[0] - arr_start[0]
            dy = arr_end[1] - arr_start[1]
            plt.arrow(arr_start[0], arr_start[1], dx, dy, head_width=0.03, head_length=0.03, length_includes_head=True, color='black')
        # TODO optimize
        new = sorted(list({tup for v in diffs.values() for tup in v[0]}), key=lambda x:x[0])
        old = sorted(list({tup for v in diffs.values() for tup in v[1]}), key=lambda x:x[0])
        print(f"{len(old)=}")
        print(f"{len(new)=}")

        for i in range(len(new)):
            draw_arrow(plt, old[i], new[i])
        plt.scatter(*zip(*new), edgecolors='r')
        plt.scatter(*zip(*old), edgecolors='b')

        plt.title("Connection Differences")
        plt.xlabel("Nucleotide")
        plt.ylabel("Number of changes")
        plt.show()
