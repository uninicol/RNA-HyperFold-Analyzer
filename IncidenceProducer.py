from abc import ABC, abstractmethod


class IncidenceProducer(ABC):
    @abstractmethod
    def get_incidence_dict(self) -> dict:
        """
        Restituisce il dizionario di incidenza
        :return: il dizionario di incidenza
        """
        pass
