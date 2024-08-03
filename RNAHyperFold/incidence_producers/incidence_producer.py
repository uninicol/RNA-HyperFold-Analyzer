from abc import ABC, abstractmethod


class IncidenceProducer(ABC):
    """Classe astratta che definisce un produttore di dizionari di incidenza."""

    @abstractmethod
    def get_incidence_dict(self) -> dict:
        """
        Restituisce il dizionario di incidenza.

        Returns:
            dict: Il dizionario di incidenza.
        """
        pass
