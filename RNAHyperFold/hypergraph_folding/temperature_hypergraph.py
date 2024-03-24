import re
from abc import ABC, abstractmethod
from concurrent.futures.process import ProcessPoolExecutor

import hypernetx as hnx

from RNAHyperFold.incidence_producers.temperature_incidence_producer import (
    TemperatureIncidenceProducer,
)


class TemporalHypergraph(ABC):
    """Classe astratta che rappresenta un ipergrafo dinamico"""

    @abstractmethod
    def add_incidence_dict(self, incidence_dict: dict, time: int) -> None:
        pass

    @abstractmethod
    def get_time_hypergraph(self, time: int) -> hnx.Hypergraph:
        pass

    @abstractmethod
    def time_hypergraph_exists(self, time: int) -> bool:
        pass


class BasicTemporalHypergraph(TemporalHypergraph):
    """Ipergrafo dinamico standard"""

    def __init__(self):
        self.__temporal_hypergraph = {}

    def add_incidence_dict(self, incidence_dict, time):
        if incidence_dict is None or time is None:
            return
        self.__temporal_hypergraph[time] = hnx.Hypergraph(incidence_dict)
        pass

    def get_time_hypergraph(self, time):
        return self.__temporal_hypergraph[time]

    def time_hypergraph_exists(self, time):
        return time not in self.__temporal_hypergraph.keys()


class MemoryOptimizedFoldingHypergraph(TemporalHypergraph):
    """Ipergrafo dinamico ottimizzato, sotto il punto di vista della memoria, per i diversi folding dell'rna"""

    def __init__(self):
        self.__temporal_hypergraph = {}

    def add_incidence_dict(self, incidence_dict, time):
        found = False
        for temps, HG in self.__temporal_hypergraph.items():
            if incidence_dict == HG.incidence_dict:
                new_temps = (min(temps[0], time), max(temps[1], time))
                self.__temporal_hypergraph[new_temps] = self.__temporal_hypergraph[temps]
                del self.__temporal_hypergraph[temps]
                found = True
                break
        if not found:
            self.__temporal_hypergraph[(time, time)] = hnx.Hypergraph(incidence_dict)
        pass

    def get_time_hypergraph(self, time):
        for temps, HG in self.__temporal_hypergraph.items():
            if temps[0] <= time <= temps[1]:
                return HG
        return None

    def time_hypergraph_exists(self, time):
        return self.get_time_hypergraph(time) is not None


class SearchOptimizedFoldingHypergraph(TemporalHypergraph):
    """Ipergrafo dinamico ottimizzato, sotto il punto di vista della velocità di ricerca, per i diversi folding dell'rna"""

    def __init__(self):
        self.__temporal_hypergraph = {}
        self.__time_to_set = {}

    def add_incidence_dict(self, incidence_dict, time):
        found = False
        for temps, HG in self.__temporal_hypergraph.items():
            if incidence_dict == HG.incidence_dict:
                new_temps = (min(temps[0], time), max(temps[1], time))
                self.__temporal_hypergraph[new_temps] = self.__temporal_hypergraph[temps]
                del self.__temporal_hypergraph[temps]
                for t in range(new_temps[0], new_temps[1] + 1):
                    self.__time_to_set[t] = new_temps
                found = True
                break
        if not found:
            new_temp = (time, time)
            self.__temporal_hypergraph[new_temp] = hnx.Hypergraph(incidence_dict)
            self.__time_to_set[time] = new_temp
        pass

    def get_time_hypergraph(self, time):
        return self.__temporal_hypergraph[self.__time_to_set[time]]

    def time_hypergraph_exists(self, time):
        return time in self.__time_to_set.keys()


class SingleFoldingHypergraph(TemporalHypergraph):
    """Ipergrafo dinamico che utilizza un singolo ipergrafo per rappresentare ogni folding, TODO attualmente in sviluppo"""

    def __init__(self):
        self.__temporal_hypergraph = None
        self.__analyzed_temperatures = set()

    def add_incidence_dict(self, incidence_dict, time):
        self.__analyzed_temperatures.add(time)
        if not self.__temporal_hypergraph:
            self.__temporal_hypergraph = hnx.Hypergraph(incidence_dict)
            for h_arc in incidence_dict.keys():
                self.__temporal_hypergraph.properties['properties'][0][h_arc] = {'temperatures': set()}

        self.__update_hypergraph(incidence_dict, time)

    def __update_hypergraph(self, incidence_dict, time):
        hyperarc_set = {tuple(n) for n in incidence_dict.values()}
        for h_arc, nodes in incidence_dict.items():
            if h_arc.startswith('l'):
                self.__temporal_hypergraph.properties['properties'][0][h_arc]['temperatures'].add(time)
                continue
            if tuple(nodes) in hyperarc_set:
                for a, n in self.__temporal_hypergraph.incidence_dict.items():
                    if nodes == n:
                        self.__temporal_hypergraph.properties['properties'][0][a]['temperatures'].add(time)
                        break
            else:
                connection_name = h_arc.split('_', 1)[0]
                self.__temporal_hypergraph = self.__add_edge(self.__temporal_hypergraph, connection_name, nodes, time)

    def __add_edge(self, HG: hnx.Hypergraph, name_begins: str, edge: list, temperature: int):
        last_edge = max((int(match) for string in HG.incidence_dict.keys() for match in re.findall(r'\d+', string)))
        new_incidence_dict = HG.incidence_dict.copy()
        new_incidence_dict[f"{name_begins}_{last_edge}"] = edge
        HG = hnx.Hypergraph(new_incidence_dict)
        HG.properties['properties'][0][f"{name_begins}_{last_edge}"] = {'temperatures': {temperature}}
        return HG

    def get_time_hypergraph(self, time):
        time_incident_dict = {}
        for h_arc, nodes in self.__temporal_hypergraph.incidence_dict.items():
            for a, n in self.__temporal_hypergraph.incidence_dict.items():
                if nodes == n:
                    if time in self.__temporal_hypergraph.properties['properties'][0][a]['temperatures']:
                        time_incident_dict[a] = nodes
                    break
        return hnx.Hypergraph(time_incident_dict)

    def time_hypergraph_exists(self, time):
        return time in self.__analyzed_temperatures


class TemperatureFoldingHypergraph:
    """Classe che permette di computare e memorizzare i folding di diverse temperature"""

    def __init__(
            self,
            producer: TemperatureIncidenceProducer,
            temporal_hypergraph: TemporalHypergraph,
    ) -> None:
        self.__producer = producer
        self.temperature_HG = temporal_hypergraph
        self.__analyzed_temperatures = set()

    def insert_temperature(self, temperature: int) -> bool:
        """
        Computa il folding in una certa temperatura
        :return : True se il folding è stato computato, False se il folding è stato computato precedentemente
        """
        if temperature in self.__analyzed_temperatures:
            return False
        self.__analyzed_temperatures.add(temperature)
        incidence_dict = self.__producer.get_temperature_incidence_dict(temperature)
        self.temperature_HG.add_incidence_dict(incidence_dict, temperature)
        return True

    def insert_temperatures(self, temperatures: list[int]) -> None:
        for temp in temperatures.copy():
            if temp in self.__analyzed_temperatures:
                temperatures.remove(temp)
            else:
                self.__analyzed_temperatures.add(temp)
        with ProcessPoolExecutor() as executor:
            for i, incidence in enumerate(
                    executor.map(
                        self.__producer.get_temperature_incidence_dict, temperatures
                    )
            ):
                self.temperature_HG.add_incidence_dict(incidence, temperatures[i])

    def insert_temperature_range(
            self, start_temperature: int, end_temperature: int, step: int = 1
    ):
        temperatures = list(range(start_temperature, end_temperature + 1, step))
        self.insert_temperatures(temperatures)

    def get_hypergraph(self, temperature: int) -> hnx.Hypergraph:
        return self.temperature_HG.get_time_hypergraph(temperature)
