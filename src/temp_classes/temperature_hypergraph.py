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


class TemperatureFoldingHypergraph:

    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.__producer = producer
        self.__temperature_HG = TemporalHypergraph()
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
