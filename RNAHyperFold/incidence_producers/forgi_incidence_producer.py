import forgi

from RNAHyperFold.incidence_producers.forna_incidence_producer import (
    FornaIncidenceProducer,
)


class ForgiIncidenceProducer(FornaIncidenceProducer):
    """
    Classe che estende FornaIncidenceProducer per produrre dizionari di incidenza
    utilizzando la libreria forgi.
    """

    def __init__(self, forna_file_path: str) -> None:
        """
        Inizializza un'istanza della classe ForgiIncidenceProducer.

        Args:
            forna_file_path (str): Il percorso del file Forna.
        """
        super().__init__(forna_file_path)

    def structure_connections(self) -> None:
        """
        Collega le strutture utilizzando un dizionario di incidenza.
        """
        structures = self.get_structures()
        for i in range(len(structures[0])):
            self.incidence_dict[f"{structures[0][i]}_{structures[1][i]}"].append(i)

    def get_structures(self) -> list[str]:
        """
        Ottiene le strutture dell'RNA dal file Forna.

        Returns:
            list[str]: Una lista di stringhe che rappresentano le strutture dell'RNA.
        """
        cg = forgi.load_rna(self.molecule["dotbracket"], allow_many=False)
        structures = cg.to_element_string(with_numbers=True)
        structures = structures.split("\n")
        return structures
