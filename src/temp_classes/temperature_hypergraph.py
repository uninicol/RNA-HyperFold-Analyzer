import re
from concurrent.futures.process import ProcessPoolExecutor

import hypernetx as hnx

from src.incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


class TemperatureHypergraph:
    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.incidence_dict = {}
        self.edge_properties = {}
        self.HG: hnx.Hypergraph = None  # TODO si puo sostituire con l'incidence dict e creare HG quando ne ho bisogno
        self.producer = producer
        self.analyzed_temperatures = set()

    def __insert_temperature_to(self, old_incidence_dict: dict, old_edge_properties: dict, temperature: int):
        if temperature in self.analyzed_temperatures:
            return
        self.analyzed_temperatures.add(temperature)
        temp_incidence_dict = self.producer.get_temperature_incidence_dict(temperature)
        print(self.producer.dotbracket)
        if old_incidence_dict is None or old_incidence_dict == {}:
            edge_properties = {edge: {"temp": {temperature}} for edge in temp_incidence_dict}
        else:
            new_edges = self.__get_edges_to_update(old_incidence_dict, temp_incidence_dict)
            temp_incidence_dict = self.__merge_incidence_dicts(old_incidence_dict, temp_incidence_dict)
            edge_properties = self.__update_edge_properties(old_edge_properties, new_edges, temperature)
        return temp_incidence_dict, edge_properties

    def insert_temperature(self, temperature: int):
        self.incidence_dict, self.edge_properties = self.__insert_temperature_to(self.incidence_dict,
                                                                                 self.edge_properties,
                                                                                 temperature)
        # self.HG = hnx.Hypergraph(incidence_dict, edge_properties=edge_properties)
        pass

    def insert_temperature_range(
            self, start_temperature: int, end_temperature: int, step: int = 1
    ):
        # incidence_dict = self.HG.incidence_dict
        # edge_properties = []  # TODO modificare i metodi di costruzione proprietà per permettere questo metodo (eliminare dipendenza con self.HG)
        # with ProcessPoolExecutor() as executor:
        #     for incidence, edge_p in executor.map(
        #             self.__insert_temperature, range(start_temperature, end_temperature + 1, step)
        #     ):
        #         incidence_dict.update(incidence)
        pass

    def insert_temperatures(self, temperatures: list[int]):
        pass

    def __merge_incidence_dicts(self, old_incidence, new_edges):
        # i link che connettono un nucleotide al successivo non cambiano
        incidence_dict = {}
        # incidence_dict.update(links)
        # incidence_dict.update(new_db_dict)
        # incidence_dict.update(new_structure_dict)
        return incidence_dict

    def __get_edges_to_update(self, old_incidence, new_incidence):
        links = {
            edge: temps for edge, temps in new_incidence.items() if edge.startswith("l")
        }
        new_db_dict = self.__get_dotbracket_to_update(old_incidence, new_incidence)
        new_structure_dict = self.__get_new_structures(old_incidence, new_incidence)
        new_incidence = {}
        new_incidence.update(links)
        new_incidence.update(new_db_dict)
        new_incidence.update(new_structure_dict)
        return new_incidence

    def __get_dotbracket_to_update(self, old_incidence, new_incidence):
        db_old = {
            edge: points
            for edge, points in old_incidence.items()
            if edge.startswith("db")
        }
        db_new = {
            edge: points  # TODO trovare il modo di non trasformarlo in list
            for edge, points in new_incidence.items()
            if edge.startswith("db")
        }

        new_db_dict = {}
        # gli archi in comune mantengono lo stesso nome
        for edge, points in db_old.items():
            if points in db_new.values():
                new_db_dict[edge] = points

        # Ricavo i set nuovi da aggiungere
        new_db_sets = [s for s in db_new.values() if s not in db_old.values()]
        # Trovo il numero dell'ultimo arco db
        last_db_number = find_max_numeric_value(
            db_old.keys()
        )  # TODO trovare il modo di farlo O(1)
        if last_db_number is None:
            last_db_number = 0
        else:
            last_db_number += 1

        for points in new_db_sets:
            new_db_dict[f"db_{last_db_number}"] = points
            last_db_number += 1
        return new_db_dict

    def __get_new_structures(self, old_incidence, new_incidence):
        str_old = {
            edge: temps
            for edge, temps in old_incidence.items()
            if not edge.startswith("db") and not edge.startswith("l")
        }
        str_new = {
            edge: list(temps)  # TODO trovare il modo di non trasformarlo in list
            for edge, temps in new_incidence.items()
            if not edge.startswith("db") and not edge.startswith("l")
        }
        # creo il dizionario nome_struttura->nodi
        new_str_to_add = {e[0]: s for e, s in str_new.items() if s not in str_old.values()}
        str_last = {edge[0]: find_max_numeric_value(
            [k for k in str_old.keys() if k.startswith(edge)]) for edge in new_str_to_add.keys()}
        new_structures = {}
        for edge, s in new_str_to_add.items():
            new_str_key = f"{edge}_{str_last[edge] + 1}"
            new_structures[new_str_key] = s
            str_last[edge] += 1
        return new_structures

    def __update_edge_properties(self, edge_properties, new_edges, temperature):
        for edge, prop in edge_properties.items():
            if edge in new_edges.keys():
                prop['temp'].add(temperature)
        return edge_properties


def add_temperature(temp_set, temperature):
    temp_set.add(temperature)
    return temp_set


def extract_numeric(s):
    numeric_parts = re.findall(r'\d+', s)
    return [int(num) for num in numeric_parts]


def find_max_numeric_value(string_list):
    all_numeric_values = [val for s in string_list for val in extract_numeric(s)]
    if all_numeric_values:
        return max(all_numeric_values)
    else:
        return None


class TemporalHypergraph:

    def __init__(self):
        self.time_graphs = {}
        pass

    def add_incidence_dict(self, incidence_dict, time):
        if incidence_dict is None:
            return
        if time in self.time_graphs.keys():
            return
        self.time_graphs[time] = hnx.Hypergraph(incidence_dict)
        pass


class TemperatureFoldingHypergraph:

    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.producer = producer
        self.temperature_HG = TemporalHypergraph()
        self.analyzed_temperatures = set()
        pass

    def generate_temperature_incidence_dict(self, temperature):
        return self.producer.get_temperature_incidence_dict(temperature)

    def insert_temperature(self, temperature):
        if temperature in self.analyzed_temperatures:
            pass
        self.analyzed_temperatures.add(temperature)
        incidence_dict = self.generate_temperature_incidence_dict(temperature)
        self.temperature_HG.add_incidence_dict(incidence_dict, temperature)
        pass

    def insert_temperatures(self, temperatures: list[int]):
        for temp in temperatures.copy():
            if temp in self.analyzed_temperatures:
                temperatures.remove(temp)
            else:
                self.analyzed_temperatures.add(temp)
        with ProcessPoolExecutor() as executor:
            for i, incidence in enumerate(executor.map(self.generate_temperature_incidence_dict, temperatures)):
                self.temperature_HG.add_incidence_dict(incidence, temperatures[i])
        pass

    def insert_temperature_range(self, start_temperature: int, end_temperature: int, step: int = 1):
        temperatures = list(range(start_temperature, end_temperature + 1, step))
        self.insert_temperatures(temperatures)
        pass
