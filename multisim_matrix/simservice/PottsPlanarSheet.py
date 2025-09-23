import numpy as np
from multisim_matrix.simservice.PlanarSheetSimService import PlanarSheetSimService
from typing import Dict, Optional, Tuple


from cc3d.core.simservice.CC3DSimService import CC3DSimService
from cc3d.core import PyCoreSpecs as pcs
from cc3d.core.iterators import CellList, CellNeighborListFlex
from cc3d.core.PySteppables import MitosisSteppableBase


def core_specs(num_cells_x: int,
               num_cells_y: int,
               cell_radius: int):
    area = np.pi * cell_radius * cell_radius
    cell_len = int(np.sqrt(area))
    area = int(area)

    dim_x = num_cells_x * cell_len
    dim_y = num_cells_y * cell_len

    return [
        pcs.PottsCore(dim_x=dim_x,
                      dim_y=dim_y,
                      fluctuation_amplitude=10,
                      neighbor_order=2),
        pcs.CellTypePlugin('TypeA'),
        pcs.ContactPlugin(4,
                          pcs.ContactEnergyParameter('Medium', 'TypeA', 5),
                          pcs.ContactEnergyParameter('TypeA', 'TypeA', 5)),
        pcs.VolumePlugin(pcs.VolumeEnergyParameter('TypeA', area, 5)),
        pcs.UniformInitializer(pcs.UniformInitializerRegion((0, 0, 0),
                                                            (dim_x, dim_y, 1),
                                                            width=cell_len,
                                                            cell_types=['TypeA'])),
        pcs.NeighborTrackerPlugin()
    ]


DEF_OUTPUT_FREQUENCY = 0
DEF_SCREENSHOT_OUTPUT_FREQUENCY = 0
DEF_OUTPUT_DIR = None
DEF_OUTPUT_FILE_CORE_NAME = None


class PottsPlanarSheet(CC3DSimService, PlanarSheetSimService):
    """Maps internal and external cell IDs such that an external ID of 0 is the first cell"""

    def __init__(self,
                 num_cells_x: int,
                 num_cells_y: int,
                 cell_radius,
                 output_frequency=DEF_OUTPUT_FREQUENCY,
                 screenshot_output_frequency=DEF_SCREENSHOT_OUTPUT_FREQUENCY,
                 output_dir=DEF_OUTPUT_DIR,
                 output_file_core_name=DEF_OUTPUT_FILE_CORE_NAME):

        PlanarSheetSimService.__init__(self, num_cells_x, num_cells_y, cell_radius)
        CC3DSimService.__init__(self,
                                output_frequency=output_frequency,
                                screenshot_output_frequency=screenshot_output_frequency,
                                output_dir=output_dir,
                                output_file_core_name=output_file_core_name,
                                sim_name='PlanarSheet')

        self.register_specs(core_specs(num_cells_x, num_cells_y, cell_radius))

        self.register_steppable(MitosisSteppableBase)
        self.mitosis_steppable = None

        self._cell_id_map: Optional[Dict[str, int]] = None
        self._cell_id_map_inv: Optional[Dict[int, str]] = None

    def start(self) -> bool:
        from cc3d import CompuCellSetup
        steppable_registry = CompuCellSetup.persistent_globals.steppable_registry
        self.mitosis_steppable = steppable_registry.getSteppablesByClassName('MitosisSteppableBase')[0]

        result = super().start()

        cinv = PottsPlanarSheet._get_cell_inventory()
        if cinv is None:
            return result
        for cell in CellList(cinv):
            self._process_new_cell(cell)

        return result

    @staticmethod
    def _get_simulator():
        from cc3d.CompuCellSetup import persistent_globals as pg
        if pg.simulator is None:
            return None
        return pg.simulator

    @staticmethod
    def _get_potts():
        sim = PottsPlanarSheet._get_simulator()
        if sim is None:
            return None
        return sim.getPotts()

    @staticmethod
    def _get_cell_field():
        potts = PottsPlanarSheet._get_potts()
        if potts is None:
            return None
        return potts.getCellFieldG()

    @staticmethod
    def _get_cell_inventory():
        potts = PottsPlanarSheet._get_potts()
        if potts is None:
            return None
        return potts.getCellInventory()

    @staticmethod
    def _get_neighbor_tracker_plugin():
        from cc3d.cpp import CompuCell
        return CompuCell.getNeighborTrackerPlugin()

    @staticmethod
    def _get_cell_by_id(_id: int):
        cinv = PottsPlanarSheet._get_cell_inventory()
        if cinv is None:
            return None
        return cinv.attemptFetchingCellById(_id)

    def cell_spatial_data(self):
        cell_field = self._get_cell_field()
        dim = cell_field.getDim()
        x = np.zeros((dim.x, dim.y), dtype=int)
        for i in range(dim.x):
            for j in range(dim.y):
                cell = cell_field[i, j, 0]
                if cell is not None:
                    x[i, j] = cell.id
        return x, dim.x, dim.y

    # PlanarSheetSimService interface

    @classmethod
    def init_arginfo(cls):
        return []

    @classmethod
    def init_kwarginfo(cls):
        """
        Returns information about implementation initialization keyword arguments

        * keyword
        * description
        * type
        * optional flag
        * default value if optional, otherwise None
        """
        result = cls.default_init_kwarginfo()
        result.extend([
            ('output_frequency', 'Steps per data output', int.__name__, True, DEF_OUTPUT_FREQUENCY),
            ('screenshot_output_frequency', 'Steps per screenshot output', int.__name__, True,
             DEF_SCREENSHOT_OUTPUT_FREQUENCY),
            ('output_dir', 'Data output directory', str.__name__, True, DEF_OUTPUT_DIR),
            ('output_file_core_name', 'Data output file core name', str.__name__, True, DEF_OUTPUT_FILE_CORE_NAME)
        ])
        return result

    def _set_id_map(self, _cell, _external_id: str):
        if self._cell_id_map is None:
            self._cell_id_map = dict()
        if self._cell_id_map_inv is None:
            self._cell_id_map_inv = dict()

        if _cell.id in self._cell_id_map_inv:
            self._cell_id_map.pop(self._cell_id_map_inv[_cell.id])

        self._cell_id_map[_external_id] = _cell.id
        self._cell_id_map_inv[_cell.id] = _external_id

    def _process_new_cell(self, _cell, new_id: str = None):
        if new_id is None:
            new_id = str(_cell.id - 1)
        self._set_id_map(_cell, new_id)

    def _neighbor_surface_areas(self, _cell) -> Dict[str, float]:
        cinv = PottsPlanarSheet._get_cell_inventory()
        result = {}
        if cinv is None:
            return result

        neighbor_tracker_plugin = PottsPlanarSheet._get_neighbor_tracker_plugin()
        if neighbor_tracker_plugin is None:
            return result

        for nbs, csa in CellNeighborListFlex(neighbor_tracker_plugin, _cell):
            if nbs:
                result[self._cell_id_map_inv[nbs.id]] = float(csa)
        return result

    def neighbor_surface_areas(self) -> Dict[str, Dict[str, float]]:
        result = {}
        cinv = PottsPlanarSheet._get_cell_inventory()
        if cinv is None:
            return result
        for cell in CellList(cinv):
            result[self._cell_id_map_inv[cell.id]] = self._neighbor_surface_areas(cell)
        return result

    def num_cells(self) -> int:
        potts = PottsPlanarSheet._get_potts()
        if potts is None:
            return 0
        return potts.getNumCells()

    def cell_volumes(self) -> Dict[str, float]:
        result = {}
        cinv = PottsPlanarSheet._get_cell_inventory()
        if cinv is None:
            return result
        for cell in CellList(cinv):
            result[self._cell_id_map_inv[cell.id]] = float(cell.volume)
        return result

    def set_cell_volume_targets(self, _targets: Dict[str, float]) -> None:
        for cell_id, cell_volume in _targets.items():
            cell = self._get_cell_by_id(self._cell_id_map[cell_id])
            if cell is not None:
                cell.targetVolume = cell_volume

    def divide_cells(self, _ids: Dict[str, Tuple[str, str]]) -> Dict[str, str]:
        result = {}
        for cell_id, (new_parent_id, child_id) in _ids.items():
            cell = self._get_cell_by_id(self._cell_id_map[cell_id])
            if cell is None:
                continue
            self.mitosis_steppable.divide_cell_random_orientation(cell)
            new_cell = self._get_cell_by_id(self.mitosis_steppable.child_cell)
            self._process_new_cell(new_cell, new_id=child_id)
            self._set_id_map(cell, new_parent_id)
            result[self._cell_id_map_inv[cell.id]] = self._cell_id_map_inv[new_cell.id]
        return result


def test():
    sim = PottsPlanarSheet(10, 10, 3)
    sim.run()
    sim.init()
    sim.start()
    sim.visualize()
    input('Press any key to continue...')
    print('Done!')


if __name__ == '__main__':
    test()
