import re

import hypernetx as hnx

from src.incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


class TemperatureHypergraph:
    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.HG: hnx.Hypergraph = None
        self.producer = producer
        self.analyzed_temperatures = set()

    def insert_temperature(self, temperature: int):
        if temperature in self.analyzed_temperatures:
            return
        self.analyzed_temperatures.add(temperature)
        incidence_dict = self.producer.get_temperature_incidence_dict(temperature)
        if self.HG is None:
            edge_properties = {edge: {"temp": {temperature}} for edge in incidence_dict}
        else:
            # TODO aggiornare anche l'incidence dict
            incidence_dict = self.__update_incidence_dict(incidence_dict)
            edge_properties = self.__update_edge_properties(incidence_dict)
            # incidence_dict.update(new_structures_connections)

        # edge_properties = {edge: {'temp': {temperature}} for edge in incidence_dict if
        #                        edge not in self.H.incidence_dict}
        #     edge_properties.update({edge: self.H.get_properties(edge)['properties']['temp'].add(temperature) for edge in
        #                             incidence_dict})
        self.HG = hnx.Hypergraph(incidence_dict, edge_properties=edge_properties)
        # print(self.H.edge_properties)
        pass

    def insert_temperature_range(
            self, start_temperature: int, end_temperature: int, step: int = 1
    ):
        # with Pool() as pool:
        #     for result in pool.imap(
        #         self.producer().temperature_folding, range(start_temperature, end_temperature + 1, step)
        #     ):
        pass

    def insert_temperatures(self, temperatures: list[int]):
        # with Pool() as pool:
        #     for result in pool.imap(
        #         self.producer().temperature_folding, temperatures
        #     ):
        pass

    def __update_incidence_dict(self, new_incidence):
        # i link che connettono un nucleotide al successivo non cambiano
        links = {
            edge: temps for edge, temps in new_incidence.items() if edge.startswith("l")
        }
        new_db_dict = self.__get_new_dotbracket(new_incidence)
        new_structure_dict = self.__get_new_structures(new_incidence)
        incidence_dict = {}
        incidence_dict.update(links)
        incidence_dict.update(new_db_dict)
        incidence_dict.update(new_structure_dict)
        return incidence_dict

    def __get_new_dotbracket(self, new_incidence):
        db_old = [
            temps
            for edge, temps in self.HG.incidence_dict.items()
            if edge.startswith("db")
        ]
        db_new = [
            list(temps)  # TODO trovare il modo di trasformarlo in list
            for edge, temps in new_incidence.items()
            if edge.startswith("db")
        ]
        new_db_sets = [s for s in db_new if s not in db_old]
        last_db_number = find_max_numeric_value(
            self.HG.incidence_dict.keys()
        )  # TODO trovare il modo di farlo O(1)
        last_db_number += 1
        new_db_dict = self.HG.incidence_dict
        for s in new_db_sets:
            new_db_dict[f"db_{last_db_number}"] = s
            last_db_number += 1
        return new_db_dict

    def __get_new_structures(self, new_incidence):
        str_old = {
            edge: temps
            for edge, temps in self.HG.incidence_dict.items()
            if not edge.startswith("db") and not edge.startswith("l")
        }
        str_new = {
            edge: list(temps)  # TODO trovare il modo di trasformarlo in list
            for edge, temps in new_incidence.items()
            if not edge.startswith("db") and not edge.startswith("l")
        }
        new_str_to_add = {e[0]: s for e, s in str_new.items() if s not in str_old.values()}
        str_last = {edge[0]: find_max_numeric_value(
            [k for k in str_old.keys() if k.startswith(edge[0])]) for edge in new_str_to_add.keys()}
        new_structures = {}
        for edge, s in new_str_to_add.items():
            new_str_key = f"{edge}_{str_last[edge] + 1}"
            new_structures[new_str_key] = s
            str_last[edge] += 1
        return new_structures

    def __update_edge_properties(self, incidence_dict):
        pass


def add_temperature(temp_set, temperature):
    temp_set.add(temperature)
    return temp_set


def extract_numeric(s):
    numeric_parts = re.findall(r'\d+', s)
    return [int(num) for num in numeric_parts]


def find_max_numeric_value(string_list):
    all_numeric_values = [val for s in string_list for val in extract_numeric(s)]
    if all_numeric_values:
        max_value = max(all_numeric_values)
        return max_value
    else:
        return None
