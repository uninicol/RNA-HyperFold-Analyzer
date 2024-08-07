from collections import defaultdict
from collections import deque

import forgi

from RNAHyperFold.hypergraph_folding.rna_folder import RNAFolder
from RNAHyperFold.incidence_producers.connector import Connector
from RNAHyperFold.incidence_producers.temperature_incidence_producer import (
    TemperatureIncidenceProducer,
)


class ViennaIncidenceProducer(TemperatureIncidenceProducer, Connector):
    """
    Produce un dizionario di incidenza che rappresenta una sequenza di RNA come ipergrafo
    dato un file JSON di Forna.
    """

    def __init__(self, folder: RNAFolder) -> None:
        """
        Inizializza un'istanza della classe ViennaIncidenceProducer.

        Args:
            folder (RNAFolder): L'oggetto RNAFolder contenente la sequenza di RNA.
        """
        self.folder: RNAFolder = folder
        self.sequence: str = folder.sequence
        self.dotbracket: str = None
        self.incidence_dict: defaultdict = defaultdict(list)

    def get_temperature_incidence_dict(self, temperature: int) -> dict:
        """
        Restituisce il dizionario di incidenza per una data temperatura.

        Args:
            temperature (int): La temperatura per cui ottenere il dizionario di incidenza.

        Returns:
            dict: Il dizionario di incidenza.
        """
        self.folder.set_temperature(temperature)
        self.dotbracket = self.folder.get_dot_bracket()
        self.incidence_dict.clear()
        self.connect_to_next()
        self.dotbracket_connections()
        self.structure_connections()
        return self.incidence_dict

    def connect_to_next(self) -> None:
        """Collega ogni nucleotide con il suo successivo"""
        edge: int = 0
        for i in range(len(self.dotbracket) - 1):
            self.incidence_dict[f"l_{edge}"] = [i, i + 1]
            edge += 1

    def dotbracket_connections(self) -> None:
        """Collega i nucleotidi in base alla rappresentazione punto-parentesi"""
        edge: int = 0
        stack = deque()
        for i, value in enumerate(self.dotbracket):
            if value == "(":
                stack.append(i)
            elif value == ")":
                if len(stack) == 0:
                    raise (ValueError("Closing bracket not matching"))

                start = stack.pop()
                self.incidence_dict[f"db_{edge}"] = [start, i]
                # self.incidence_dict[f"db_{edge}"].append(i)
                edge += 1

    def structure_connections(self):
        """Collega le strutture rilevate da forna"""
        structures = self.get_structures()
        for struct, indexes in structures.items():
            self.incidence_dict[struct] = indexes

    def get_structures(self) -> dict:
        """
        Ottiene le strutture dell'RNA dal file Forna.

        Returns:
            dict: Un dizionario che rappresenta le strutture dell'RNA.
        """
        structures_dict = defaultdict(list)
        cg = forgi.load_rna(self.dotbracket, allow_many=False)
        structures = cg.to_element_string(with_numbers=True)
        structures = structures.split("\n")
        for i in range(len(structures[0])):
            structures_dict[f"{structures[0][i]}_{structures[1][i]}"].append(i)
        return structures_dict
        # TODO fare questa operazione con viennarna per non importare forgi
        # structures_dict = defaultdict(list)
        # structures = RNA.db_to_element_string(self.rna_dotbracket)
        # print(structures)
        # structures = structures.split("\n")
        # for i in range(len(structures[0])):
        #     structures_dict[f"{structures[0][i]}_{structures[1][i]}"].append(i)
        # return structures_dict
