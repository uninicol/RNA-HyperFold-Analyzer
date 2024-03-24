import forgi

from RNAHyperFold.incidence_producers.forna_incidence_producer import FornaIncidenceProducer


class ForgiIncidenceProducer(FornaIncidenceProducer):
    def __init__(self, forna_file_path: str) -> None:
        super().__init__(forna_file_path)

    def structure_connections(self):
        structures = self.get_structures()
        for i in range(len(structures[0])):
            self.incidence_dict[f"{structures[0][i]}_{structures[1][i]}"].append(i)

    def get_structures(self):
        cg = forgi.load_rna(self.molecule["dotbracket"], allow_many=False)
        structures = cg.to_element_string(with_numbers=True)
        structures = structures.split("\n")
        return structures
