from concurrent.futures.process import ProcessPoolExecutor

import hypernetx as hnx

from src.incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


class TemporalHypergraph:

    def __init__(self):
        self.time_hypergraphs = {}

    def add_incidence_dict(self, incidence_dict, time):
        if incidence_dict is None or time is None:
            return
        self.time_hypergraphs[time] = hnx.Hypergraph(incidence_dict)
        pass


class MemoryOptimizedFoldingHypergraph(TemporalHypergraph):

    def __init__(self):
        super().__init__()

    def add_incidence_dict(self, incidence_dict, time):
        found = False
        for k, v in self.time_hypergraphs.items():
            if incidence_dict == v.incidence_dict:
                s = k.union([time])
                self.time_hypergraphs[s] = self.time_hypergraphs[k]
                del self.time_hypergraphs[k]
                found = True
                break
        if not found:
            self.time_hypergraphs[frozenset([time])] = hnx.Hypergraph(incidence_dict)
        pass


class TemperatureFoldingHypergraph:

    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.__producer = producer
        self.__temperature_HG = MemoryOptimizedFoldingHypergraph()
        self.__analyzed_temperatures = set()

    def insert_temperature(self, temperature):
        if temperature in self.__analyzed_temperatures:
            return
        self.__analyzed_temperatures.add(temperature)
        incidence_dict = self.__producer.get_temperature_incidence_dict(temperature)
        self.__temperature_HG.add_incidence_dict(incidence_dict, temperature)

    def insert_temperatures(self, temperatures: list[int]):
        for temp in temperatures.copy():
            if temp in self.__analyzed_temperatures:
                temperatures.remove(temp)
            else:
                self.__analyzed_temperatures.add(temp)
        with ProcessPoolExecutor() as executor:
            for i, incidence in enumerate(executor.map(self.__producer.get_temperature_incidence_dict, temperatures)):
                self.__temperature_HG.add_incidence_dict(incidence, temperatures[i])

    def insert_temperature_range(self, start_temperature: int, end_temperature: int, step: int = 1):
        temperatures = list(range(start_temperature, end_temperature + 1, step))
        self.insert_temperatures(temperatures)
