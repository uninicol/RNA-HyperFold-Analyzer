import subprocess

from FornaIncidenceProducer import FornaIncidenceProducer


class ForgiIncidenceProducer(FornaIncidenceProducer):
    def __init__(self, forna_file_path: str) -> None:
        super().__init__(forna_file_path)

    def structure_connections(self, structures):
        structures = self.get_structures()
        for i in range(len(structures[0])):
            self.incidence_dict[f"{structures[1][i]}_{structures[2][i]}"].append(i)

    def get_structures(self):
        command = f"python3 forgi/rnaConvert.py {self.molecule['dotbracket']} -T element_string".split()
        result = subprocess.run(command, stdout=subprocess.PIPE)
        result = result.stdout.decode("utf-8")
        result = result.split("\n")
        return result
