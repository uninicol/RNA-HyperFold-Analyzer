from collections import defaultdict
from collections import deque

import forgi
from ViennaRNA import fold

from incidence_producers.Connector import Connector
from incidence_producers.IncidenceProducer import IncidenceProducer


class ViennaIncidenceProducer(IncidenceProducer, Connector):
    """Produce un dizionario di incidenza che rappresenta una sequenza di RNA come ipergrafo dato un file json di
    forna"""

    def __init__(self, rna_seq: str) -> None:
        self.rna_seq = rna_seq
        self.rna_dotbracket, _ = fold(rna_seq)
        self.incidence_dict: defaultdict = defaultdict(list)
        self.edge: int = 0

    def get_incidence_dict(self, node_with_nucleotide: bool = False) -> dict:
        """
        Restituisce il dizionario di incidenza
        :return: il dizionario di incidenza
        """
        self.connect_to_next()
        self.dotbracket_connections()
        self.structure_connections()
        if node_with_nucleotide:
            self.nodes_to_nucleotide_string()
        return self.incidence_dict

    def connect_to_next(self) -> None:
        """Collega ogni nucleotide con il suo successivo"""
        for i in range(len(self.rna_dotbracket) - 1):
            self.incidence_dict[f"l_{self.edge}"].append(i)
            self.incidence_dict[f"l_{self.edge}"].append(i + 1)
            self.edge += 1

    def dotbracket_connections(self) -> None:
        """Collega i nucleotidi in base alla rappresentazione punto-parentesi"""
        stack = deque()
        for i, value in enumerate(self.rna_dotbracket):
            if value == "(":
                stack.append(i)
            elif value == ")":
                if len(stack) == 0:
                    raise (ValueError("Closing bracket not matching"))

                start = stack.pop()
                self.incidence_dict[f"l_{self.edge}"].append(start)
                self.incidence_dict[f"l_{self.edge}"].append(i)
                self.edge += 1

    def structure_connections(self):
        """Collega le strutture rilevate da forna"""
        structures = self.get_structures()
        for i in range(len(structures[0])):
            self.incidence_dict[f"{structures[0][i]}_{structures[1][i]}"].append(i)

    def get_structures(self):
        cg = forgi.load_rna(self.rna_dotbracket, allow_many=False)
        structures = cg.to_element_string(with_numbers=True)
        structures = structures.split("\n")
        return structures

    def nodes_to_nucleotide_string(self) -> None:
        """Converte i nomi dei nodi nel formato "{indice}_{nucleotide corrispettivo}\" """
        for key in self.incidence_dict.keys():
            self.incidence_dict[key] = [
                f"{i}_{self.rna_seq[i]}" for i in self.incidence_dict[key]
            ]
