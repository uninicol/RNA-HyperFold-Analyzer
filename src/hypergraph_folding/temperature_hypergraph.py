from abc import ABC, abstractmethod
from concurrent.futures.process import ProcessPoolExecutor

import hypernetx as hnx

from incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


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
        self.__time_hypergraphs = {}

    def add_incidence_dict(self, incidence_dict, time):
        if incidence_dict is None or time is None:
            return
        self.__time_hypergraphs[time] = hnx.Hypergraph(incidence_dict)
        pass

    def get_time_hypergraph(self, time):
        return self.__time_hypergraphs[time]

    def time_hypergraph_exists(self, time):
        return time not in self.__time_hypergraphs.keys()


class MemoryOptimizedFoldingHypergraph(TemporalHypergraph):
    """Ipergrafo dinamico ottimizzato, sotto il punto di vista della memoria, per i diversi folding dell'rna"""

    def __init__(self):
        self.__time_hypergraphs = {}

    def add_incidence_dict(self, incidence_dict, time):
        found = False
        for temps, HG in self.__time_hypergraphs.items():
            if incidence_dict == HG.incidence_dict:
                new_temps = temps.union([time])
                self.__time_hypergraphs[new_temps] = self.__time_hypergraphs[temps]
                del self.__time_hypergraphs[temps]
                found = True
                break
        if not found:
            self.__time_hypergraphs[frozenset([time])] = hnx.Hypergraph(incidence_dict)
        pass

    def get_time_hypergraph(self, time):
        for temps, HG in self.__time_hypergraphs.items():
            if time in temps:
                return HG
        return None

    def time_hypergraph_exists(self, time):
        for temps, HG in self.__time_hypergraphs.items():
            if time in temps:
                return True
        return False


class SearchOptimizedFoldingHypergraph(TemporalHypergraph):
    """Ipergrafo dinamico ottimizzato, sotto il punto di vista della velocità di ricerca, per i diversi folding dell'rna"""

    def __init__(self):
        self.__time_hypergraphs = {}
        self.__time_to_set = {}

    def add_incidence_dict(self, incidence_dict, time):
        found = False
        for temps, HG in self.__time_hypergraphs.items():
            if incidence_dict == HG.incidence_dict:
                new_temps = temps.union([time])
                self.__time_hypergraphs[new_temps] = self.__time_hypergraphs[temps]
                del self.__time_hypergraphs[temps]
                for t in new_temps:
                    self.__time_to_set[t] = new_temps
                found = True
                break
        if not found:
            new_temp = frozenset([time])
            self.__time_hypergraphs[new_temp] = hnx.Hypergraph(incidence_dict)
            self.__time_to_set[time] = new_temp
        pass

    def get_time_hypergraph(self, time):
        return self.__time_hypergraphs[self.__time_to_set[time]]

    def time_hypergraph_exists(self, time):
        return time in self.__time_to_set.keys()


class TemperatureFoldingHypergraph:
    """Classe che permette di computare e memorizzare i folding di diverse temperature"""

    def __init__(self, producer: TemperatureIncidenceProducer, temporal_hypergraph: TemporalHypergraph) -> None:
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
            for i, incidence in enumerate(executor.map(self.__producer.get_temperature_incidence_dict, temperatures)):
                self.temperature_HG.add_incidence_dict(incidence, temperatures[i])

    def insert_temperature_range(self, start_temperature: int, end_temperature: int, step: int = 1):
        temperatures = list(range(start_temperature, end_temperature + 1, step))
        self.insert_temperatures(temperatures)

    def get_hypergraph(self, temperature: int) -> hnx.Hypergraph:
        return self.temperature_HG.get_time_hypergraph(temperature)
