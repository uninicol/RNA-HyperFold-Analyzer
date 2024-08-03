from abc import abstractmethod

from RNAHyperFold.incidence_producers.incidence_producer import IncidenceProducer


class TemperatureIncidenceProducer(IncidenceProducer):
    """
    Classe astratta che definisce un produttore di dizionari di incidenza
    basati sulla temperatura.
    """

    @abstractmethod
    def get_temperature_incidence_dict(self, temperature: int) -> dict:
        """
        Restituisce il dizionario di incidenza per una data temperatura.

        Args:
            temperature (int): La temperatura per cui ottenere il dizionario di incidenza.

        Returns:
            dict: Il dizionario di incidenza.
        """
        pass

    def get_incidence_dict(self) -> dict:
        """
        Restituisce il dizionario di incidenza per la temperatura predefinita di 37°C.

        Returns:
            dict: Il dizionario di incidenza.
        """
        return self.get_temperature_incidence_dict(temperature=37)
