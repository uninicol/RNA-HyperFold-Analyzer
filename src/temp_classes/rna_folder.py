from ViennaRNA import RNA, fold


class RNAFolder:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence
        pass

    def temperature_folding(self, temperature):
        RNA.cvar.temperature = temperature
        dot_bracket, _ = fold(self.sequence)
        return dot_bracket
