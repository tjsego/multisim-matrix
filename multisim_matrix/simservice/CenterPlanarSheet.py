import numpy as np
from simservice.PySimService import PySimService
from multisim_matrix.simservice.PlanarSheetSimService import PlanarSheetSimService
import tissue_forge as tf
from typing import Dict, Optional

# Total area: pi * d
# Nominal neighborhood: 6 neighbors perfectly in contact
# Note that TF searches for neighbors from the surface of a particle of interest

neighbor_cutoff_cd = 1.5


def neighbor_area(cell_diameter: float, dist: float):
    mag = np.pi * cell_diameter / 6

    cutoff = neighbor_cutoff_cd * cell_diameter
    if dist < cell_diameter:
        return mag
    elif dist > cutoff:
        return 0.0

    cf = 1.0 - (dist - cell_diameter) / (cutoff - cell_diameter)
    return mag * cf * cf


DEF_STEP_SIZE = 1.0
DEF_DT = 0.01
DEF_SHOW = False


class CenterPlanarSheet(PlanarSheetSimService):

    def __init__(self,
                 num_cells_x: int,
                 num_cells_y: int,
                 cell_radius,
                 step_size=DEF_STEP_SIZE,
                 dt=DEF_DT,
                 show=DEF_SHOW):

        PySimService.__init__(self,
                              sim_name='PlanarSheet')
        PlanarSheetSimService.__init__(self,
                                       num_cells_x,
                                       num_cells_y,
                                       cell_radius)

        self._step_size = step_size
        self._dt = dt
        self._show = show
        self._cell_type: Optional[tf.ParticleType] = None

        self._cell_id_map: Optional[Dict[str, int]] = None
        self._cell_id_map_inv: Optional[Dict[int, str]] = None

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
            ('step_size', 'Simulation time per service step', float.__name__, True, DEF_STEP_SIZE),
            ('dt', 'Simulation time per simulation step', float.__name__, True, DEF_DT),
            ('show', 'Flag to show simulation in real-time', bool.__name__, True, DEF_SHOW)
        ])
        return result

    def _neighbor_surface_areas(self, _ph: tf.ParticleHandle) -> Dict[str, float]:
        if not self._cell_type:
            return {}
        cell_diameter = _ph.radius * 2
        return {self._cell_id_map_inv[nh.id]: neighbor_area(cell_diameter, _ph.relativePosition(nh.position).length())
                for nh in _ph.neighbors(distance=neighbor_cutoff_cd * cell_diameter - _ph.radius)}

    def _set_id_map(self, _cell, _external_id: str):
        if self._cell_id_map is None:
            self._cell_id_map = dict()
        if self._cell_id_map_inv is None:
            self._cell_id_map_inv = dict()

        if _cell.id in self._cell_id_map_inv:
            self._cell_id_map.pop(self._cell_id_map_inv[_cell.id])

        self._cell_id_map[_external_id] = _cell.id
        self._cell_id_map_inv[_cell.id] = _external_id

    def _process_new_cell(self, _ph: tf.ParticleHandle, new_id: str = None):
        if new_id is None:
            new_id = str(_ph.id)
        self._set_id_map(_ph, new_id)

        return new_id

    def neighbor_surface_areas(self) -> Dict[str, Dict[int, float]]:
        return {self._cell_id_map_inv[ph.id]: self._neighbor_surface_areas(ph) for ph in tf.Universe.particles}

    def num_cells(self) -> int:
        return len(tf.Universe.particles)

    def cell_spatial_data(self):
        pos_x = []
        pos_y = []
        cell_ids = []
        for ph in tf.Universe.particles:
            p_pos_x, p_pos_y = ph.position.xy().as_list()
            pos_x.append(p_pos_x)
            pos_y.append(p_pos_y)
            cell_ids.append(self._cell_id_map_inv[ph.id])
        return pos_x, pos_y, cell_ids, *tf.Universe.dim.xy().as_list(), self._cell_type.radius

    # PySimService interface

    def _run(self) -> None:
        pass

    def _init(self) -> bool:

        # Initialize the simulation

        pad = 2.0 * self.cell_radius

        dim_x = 2.0 * pad + (self.num_cells_x - 1) * self.cell_radius * 2
        dim_y = 2.0 * pad + ((self.num_cells_y - 1) * 2 + 1) * 2 * self.cell_radius / np.sqrt(3)
        min_dim = min(dim_x, dim_y)
        min_cells = 3
        len2cells = min_dim / min_cells

        dim = [dim_x, dim_y, len2cells]

        tf.init(windowless=not self._show,
                dim=dim,
                bc={'x': tf.BOUNDARY_FREESLIP, 'y': tf.BOUNDARY_FREESLIP, 'z': tf.BOUNDARY_FREESLIP},
                cells=[int(d / len2cells) + 1 for d in dim],
                dt=self._dt,
                cutoff=2 * self.cell_radius * 1.5)

        # Create a cell type

        self._cell_type = tf.ParticleType()
        self._cell_type.thisown = 0
        self._cell_type.radius = self.cell_radius
        self._cell_type.dynamics = tf.Overdamped
        self._cell_type.frozen_z = True

        # Add potentials

        pot_contact = tf.Potential.morse(d=1E-4,
                                         a=3,
                                         r0=2 * self._cell_type.radius,
                                         min=1E-6,
                                         max=4 * self._cell_type.radius,
                                         shifted=False)
        pot_contact.thisown = 0
        tf.bind.types(pot_contact, self._cell_type, self._cell_type)

        # Add some noise

        rforce = tf.Force.random(1E-2, 0.0)
        rforce.thisown = 0
        tf.bind.force(rforce, self._cell_type)

        # Initialize the population

        for i in range(self.num_cells_x):
            for j in range(self.num_cells_y):
                pos_x = pad + i * self.cell_radius * 2
                pos_y = pad + (j * 2 + (i % 2)) * 2 * self.cell_radius / np.sqrt(3)
                self._process_new_cell(self._cell_type([pos_x, pos_y, tf.Universe.center[2]]))

        if tf.err_occurred():
            print(tf.err_get_all())

        return True

    def _start(self) -> bool:
        return True

    def _step(self) -> bool:
        return tf.step(self._step_size) == 0

    def _finish(self) -> None:
        pass


def test():
    sim = CenterPlanarSheet(10, 10, 1.0, show=True)
    # sim = CenterPlanarSheet(4, 9, 1.0, show=True)
    # sim = CenterPlanarSheet(9, 4, 1.0, show=True)
    sim._init()
    tf.system.camera_view_top()
    tf.show()


if __name__ == '__main__':
    test()
