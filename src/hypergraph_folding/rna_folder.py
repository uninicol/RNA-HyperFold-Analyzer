from ViennaRNA import RNA, fold


class RNAFolder:
    """Classe che permette di computare dei folding di sequenze di rna"""
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def temperature_folding(self, temperature):
        RNA.cvar.temperature = temperature
        dot_bracket, _ = fold(self.sequence)
        return dot_bracket
