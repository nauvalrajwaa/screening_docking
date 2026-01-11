from abc import ABC, abstractmethod

class Docker(ABC):
    def __init__(self, binary_path):
        self.binary_path = binary_path

    @abstractmethod
    def dock(self, ligand_path, receptor_path, center, size, output_path, **kwargs):
        pass
