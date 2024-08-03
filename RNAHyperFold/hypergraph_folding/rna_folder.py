from ViennaRNA import RNA, fold


class RNAFolder:
    """Classe di configurazione che permette di computare dei folding di sequenze di rna"""

    def __init__(self, sequence: str) -> None:
        """
        Inizializza un'istanza della classe RNAFolder.

        Args:
            sequence (str): La sequenza di RNA da foldare.
        """
        self.sequence: str = sequence

    def set_temperature(self, temperature: int) -> None:
        """
        Imposta la temperatura per il calcolo del folding.

        Args:
            temperature (int): La temperatura da impostare.
        """
        RNA.cvar.temperature = temperature

    def get_dot_bracket(self) -> str:
        """
        Restituisce la rappresentazione dot-bracket della sequenza di RNA.

        Returns:
            str: La rappresentazione dot-bracket della sequenza di RNA.
        """
        dot_bracket, _ = fold(self.sequence)
        return dot_bracket
