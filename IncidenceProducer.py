import json
from collections import defaultdict
from collections import deque


class IncidenceProducer:
    def __init__(self, forna_file_path):
        with open(forna_file_path, "r") as json_file:
            _, self.molecule = json.load(json_file)["rnas"].popitem()
        self.incidence_dict = defaultdict(list)
        self.edge = 0

    def connect_to_next(self):
        for i in range(len(self.molecule["dotbracket"]) - 1):
            self.incidence_dict[f"l_{self.edge}"].append(i)
            self.incidence_dict[f"l_{self.edge}"].append(i + 1)
            self.edge += 1

    def dotbracket_connections(self):
        stack = deque()
        for i, value in enumerate(self.molecule["dotbracket"]):
            if value == "(":
                stack.append(i)
            elif value == ")":
                if len(stack) == 0:
                    raise (ValueError("Closing bracket not matching"))

                start = stack.pop()
                self.incidence_dict[f"l_{self.edge}"].append(start)
                self.incidence_dict[f"l_{self.edge}"].append(i)
                self.edge += 1

    def structure_connections(self, structures):
        struct_counter = defaultdict(int)
        for struct in structures:
            nucleotide = f"{struct[0]}_{struct_counter[struct[0]]}"  # {structure name letter}_{number of structure}
            # nelle strutture i nucleotidi sono numerati da 1 a n, in alcuni casi con 0 e n+1 che vengono scartati
            self.incidence_dict[nucleotide].extend(
                [i - 1 for i in struct[2] if 0 < i <= len(self.molecule["seq"])]
            )
            struct_counter[struct[0]] += 1
        pass

    def nodes_to_nucleotide_string(self):
        for key in self.incidence_dict.keys():
            self.incidence_dict[key] = [
                f"{i}_{self.molecule['seq'][i]}" for i in self.incidence_dict[key]
            ]

    def get_incidence_dict(self):
        self.connect_to_next()
        self.dotbracket_connections()
        structures = self.molecule["elements"]
        self.structure_connections(structures)
        self.nodes_to_nucleotide_string()
        return self.incidence_dict
