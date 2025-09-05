import abc
from simservice.PySimService import PySimService
from typing import Any, Dict, List, Tuple


class PlanarSheetSimService(PySimService, abc.ABC):

    def __init__(self,
                 num_cells_x: int,
                 num_cells_y: int,
                 cell_radius: float):
        self.num_cells_x = num_cells_x
        self.num_cells_y = num_cells_y
        self.cell_radius = cell_radius

    @abc.abstractclassmethod
    def init_arginfo(cls) -> List[Tuple[str, str]]:
        """
        Returns information about implementation initialization positional arguments

        * description
        * type
        """
        raise NotImplementedError

    @abc.abstractclassmethod
    def init_kwarginfo(cls) -> List[Tuple[str, str, str, bool, Any]]:
        """
        Returns information about implementation initialization keyword arguments

        * keyword
        * description
        * type
        * optional flag
        * default value if optional, otherwise None
        """
        raise NotImplementedError

    @classmethod
    def default_init_kwarginfo(cls) -> List[Tuple[str, str, str, bool, Any]]:
        return [
            ('num_cells_x', 'Number of cells along the X-axis', int.__name__, False, None),
            ('num_cells_y', 'Number of cells along the Y-axis', int.__name__, False, None),
            ('cell_radius', 'Approximate radius of each cell', float.__name__, False, None)
        ]

    @abc.abstractmethod
    def neighbor_surface_areas(self) -> Dict[int, Dict[int, float]]:
        raise NotImplementedError

    @abc.abstractmethod
    def num_cells(self) -> int:
        raise NotImplementedError

    @abc.abstractmethod
    def cell_spatial_data(self):
        raise NotImplementedError

    def cell_volumes(self) -> Dict[int, float]:
        raise NotImplementedError

    def set_cell_volume_targets(self, _targets: Dict[int, float]) -> None:
        raise NotImplementedError

    def divide_cells(self, _ids: Dict[int,Tuple[int,int]]) -> Dict[int, int]:
        """
        input from mother: (daughter1_id, daughter2_id)
        """
        raise NotImplementedError
