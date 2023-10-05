from abc import ABC, abstractmethod


class Connector(ABC):
    """Permette di costruire i collegamenti dell'rna"""

    @abstractmethod
    def connect_to_next(self) -> None:
        """Collega ogni nucleotide con il suo successivo"""
        pass

    @abstractmethod
    def dotbracket_connections(self) -> None:
        """Collega i nucleotidi in base alla rappresentazione punto-parentesi"""
        pass

    @abstractmethod
    def structure_connections(self) -> None:
        """Collega le strutture con un iperarco"""
        pass
