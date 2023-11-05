import re

import hypernetx as hnx

from src.incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


class TemperatureHypergraph:
    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.HG: hnx.Hypergraph = None  # TODO si puo sostituire con l'incidence dict e creare HG quando ne ho bisogno
        self.producer = producer
        self.analyzed_temperatures = set()

    def __insert_temperature(self, temperature: int):
        if temperature in self.analyzed_temperatures:
            return
        self.analyzed_temperatures.add(temperature)
        incidence_dict = self.producer.get_temperature_incidence_dict(temperature)
        if self.HG is None:
            edge_properties = {edge: {"temp": {temperature}} for edge in incidence_dict}
        else:
            incidence_dict = self.__merge_incidence_dicts(self.HG.incidence_dict, incidence_dict)
            edge_properties = self.__update_edge_properties(incidence_dict, temperature)  # TODO da modificare
        return incidence_dict, edge_properties

    def insert_temperature(self, temperature: int):
        incidence_dict, edge_properties = self.__insert_temperature(temperature)
        self.HG = hnx.Hypergraph(incidence_dict, edge_properties=edge_properties)
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

    def __merge_incidence_dicts(self, old_incidence, new_incidence):
        # i link che connettono un nucleotide al successivo non cambiano
        links = {
            edge: temps for edge, temps in new_incidence.items() if edge.startswith("l")
        }
        new_db_dict = self.__get_new_dotbracket(old_incidence, new_incidence)
        new_structure_dict = self.__get_new_structures(old_incidence, new_incidence)
        incidence_dict = {}
        incidence_dict.update(links)
        incidence_dict.update(new_db_dict)
        incidence_dict.update(new_structure_dict)
        return incidence_dict

    def __get_new_dotbracket(self, old_incidence, new_incidence):
        db_old = [
            temps
            for edge, temps in old_incidence.items()
            if edge.startswith("db")
        ]
        db_new = [
            list(temps)  # TODO trovare il modo di non trasformarlo in list
            for edge, temps in new_incidence.items()
            if edge.startswith("db")
        ]
        # Ricavo i set nuovi da aggiungere
        new_db_sets = [s for s in db_new if s not in db_old]
        # Trovo il numero dell'ultimo arco db
        last_db_number = find_max_numeric_value(
            old_incidence.keys()
        )  # TODO trovare il modo di farlo O(1)
        last_db_number += 1
        new_db_dict = old_incidence
        for s in new_db_sets:
            new_db_dict[f"db_{last_db_number}"] = s
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

    def __update_edge_properties(self, incidence_dict, temperature):
        edge_properties = {edge: {"temp": {temperature}} for edge in incidence_dict}
        for edge, s in self.HG.incidence_dict.items():
            edge_properties[edge]['temp'].update(self.HG.get_properties(edge)['properties']['temp'])
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
