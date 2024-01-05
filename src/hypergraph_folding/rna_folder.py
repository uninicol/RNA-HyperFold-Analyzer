from ViennaRNA import RNA, fold


class RNAFolder:
    """Classe di configurazione che permette di computare dei folding di sequenze di rna"""

    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def set_temperature(self, temperature):
        RNA.cvar.temperature = temperature

    def get_dot_bracket(self):
        dot_bracket, _ = fold(self.sequence)
        return dot_bracket
