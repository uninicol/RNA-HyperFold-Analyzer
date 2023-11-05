import hypernetx as hnx

from src.incidence_producers.temperature_incidence_producer import TemperatureIncidenceProducer


def find_unique_lists(dict1, dict2):
    # Get unique lists in each dictionary
    unique_lists_dict1 = set(tuple(lst) for lst in dict1.values())
    unique_lists_dict2 = set(tuple(lst) for lst in dict2.values())
    return unique_lists_dict1.symmetric_difference(unique_lists_dict2)


class TemperatureHypergraph:
    def __init__(self, producer: TemperatureIncidenceProducer) -> None:
        self.H: hnx.Hypergraph = None
        self.producer = producer
        self.analyzed_temperatures = set()

    def insert_temperature(self, temperature: int):
        if temperature in self.analyzed_temperatures:
            return
        self.analyzed_temperatures.add(temperature)
        incidence_dict = self.producer.get_temperature_incidence_dict(temperature)
        if self.H is None:
            edge_properties = {edge: {'temp': {temperature}} for edge in incidence_dict}
        else:
            # TODO aggiornare anche l'incidence dict
            # aggiungo i link che collegano un nucleotide al successivo
            links = {edge: temps for edge, temps in incidence_dict.items() if edge.startswith('l')}
            new_db_connections = find_unique_lists(
                {edge: temps for edge, temps in incidence_dict.items() if edge.startswith('db')},
                {edge: temps for edge, temps in self.H.incidence_dict.items() if edge.startswith('db')}
            )
            last_db = max((int(k[3:]) for k in self.H.incidence_dict.keys()))  # TODO trovare il modo di farlo O(1)
            new_structures_connections = [
                find_unique_lists(
                    {edge: temps for edge, temps in incidence_dict.items() if not edge.startswith(ignore)},
                    {edge: temps for edge, temps in self.H.incidence_dict.items() if not edge.startswith(ignore)}
                ) for ignore in ['l', 'db']
            ]
            incidence_dict = {}
            incidence_dict.update(links)
            incidence_dict.update(new_db_connections)
            # incidence_dict.update(new_structures_connections)

            # # i link in successione non cambiano mai
            # edge_properties = {
            #     edge: add_temperature(self.H.get_properties(edge)['properties']['temp'], temperature) for edge in
            #     incidence_dict if edge.startswith('l')}
            #
            # # ricavo le nuove connessioni date dal cambiamento del dotbracket
            # new_bd_connections = sets_not_in_both_dicts(
            #     {edge: temps for edge, temps in incidence_dict.items() if edge.startswith('db')},
            #     {edge: temps for edge, temps in self.H.incidence_dict.items() if edge.startswith('db')}
            # )
            # new_dotbracket_connections = find_different_lists(
            #     {edge: temps for edge, temps in incidence_dict.items() if edge.startswith('db')},
            #     {edge: temps for edge, temps in self.H.incidence_dict.items() if edge.startswith('db')})
            # # le nuove connessioni avranno lo stesso nome di quelle esistenti, bisogna quindi cambiare nome
            # db_edges = [edge for edge in self.H.incidence_dict.keys() if edge.startswith('db')]
            # last_edge = db_edges[len(db_edges) - 1][3:]
            # last_edge = int(last_edge)
            # for edge in new_dotbracket_connections.keys():
            #     last_edge += 1
            #     new_dotbracket_connections[f'db_{last_edge}'] = new_dotbracket_connections[edge]
            #     del new_dotbracket_connections[edge]
            # edge_properties.update(new_dotbracket_connections)

            # edge_properties.update(
            #     {edge: temps for edge, temps in self.H.incidence_dict.items() if edge.startswith('db')})
            # edge_properties.update({edge: temps for edge, temps in incidence_dict.items() if edge.startswith('db')})
            #
            # new_structures_connections = self.find_different_lists(
            #     {edge: temps for edge, temps in incidence_dict.items() if
            #      not edge.startswith('db') and not edge.startswith('l')},
            #     {edge: temps for edge, temps in self.H.incidence_dict.items() if
            #      not edge.startswith('db') and not edge.startswith('l')})
            # edge_properties.update(new_structures_connections)
        #
        # edge_properties = {edge: {'temp': {temperature}} for edge in incidence_dict if
        #                        edge not in self.H.incidence_dict}
        #     edge_properties.update({edge: self.H.get_properties(edge)['properties']['temp'].add(temperature) for edge in
        #                             incidence_dict})
        self.H = hnx.Hypergraph(incidence_dict, edge_properties=edge_properties)
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


def add_temperature(temp_set, temperature):
    temp_set.add(temperature)
    return temp_set


def find_different_lists(dict1, dict2):
    result = {}
    for key, list1 in dict1.items():
        is_different = True
        for _, list2 in dict2.items():
            if set(list1) == set(list2):
                is_different = False
                break
        if is_different:
            result[key] = list1
    return result
