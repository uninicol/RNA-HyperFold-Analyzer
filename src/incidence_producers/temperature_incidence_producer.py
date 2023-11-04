from abc import abstractmethod

from src.incidence_producers.incidence_producer import IncidenceProducer


class TemperatureIncidenceProducer(IncidenceProducer):
    @abstractmethod
    def get_temperature_incidence_dict(self, temperature: int) -> dict:
        """
        Restituisce il dizionario di incidenza
        :return: il dizionario di incidenza
        """
        pass

    def get_incidence_dict(self) -> dict:
        return self.get_temperature_incidence_dict(temperature=37)
