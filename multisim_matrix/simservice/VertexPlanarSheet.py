import numpy as np
from simservice.PySimService import PySimService
from multisim_matrix.simservice.PlanarSheetSimService import PlanarSheetSimService
import tissue_forge as tf
from tissue_forge.models.vertex import solver as tfvs
from typing import Dict, Optional, Tuple

DEF_STEP_SIZE = 1.0
DEF_DT = 0.01
DEF_SHOW = False


class VertexPlanarSheet(PlanarSheetSimService):

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
        self._cell_type: Optional[tfvs.SurfaceType] = None
        self._cell_id_map: Optional[Dict[str, int]] = None
        self._cell_id_map_inv: Optional[Dict[int, str]] = None

        self._init_area_constraint = None
        self._init_perimeter_constraint = None
        # todo: refine how constraints are defined

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
            ('step_size', 'Simulation time per service step', float.__name__, True, DEF_STEP_SIZE),
            ('dt', 'Simulation time per simulation step', float.__name__, True, DEF_DT),
            ('show', 'Flag to show simulation in real-time', bool.__name__, True, DEF_SHOW)
        ])
        return result

    def _neighbor_surface_areas(self, _sh: tfvs.SurfaceHandle) -> Dict[str, float]:
        if not self._cell_type:
            return {}
        result = {}
        vertices = [v for v in _sh.vertices]
        vertices.append(vertices[0])
        for i in range(len(vertices) - 1):
            va: tfvs.VertexHandle = vertices[i]
            vb = vertices[i + 1]
            dist = va.position.relative_to(vb.position, tf.Universe.dim, False, False, False).length()
            for s in va.shared_surfaces(vb):
                if s == _sh:
                    continue
                try:
                    result[self._cell_id_map_inv[s.id]] += dist
                except KeyError:
                    result[self._cell_id_map_inv[s.id]] = dist
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

    def _process_new_cell(self, _sh: tfvs.SurfaceHandle, new_id: str = None):
        if new_id is None:
            new_id = str(len(self._cell_type))
        self._set_id_map(_sh, new_id)

        if self._init_area_constraint is None:
            self._init_area_constraint = _sh.area
        if self._init_perimeter_constraint is None:
            self._init_perimeter_constraint = _sh.perimeter

        area_constraint = tfvs.SurfaceAreaConstraint(0.1, self._init_area_constraint)
        area_constraint.thisown = 0
        tfvs.bind.surface(area_constraint, _sh)

        perim_constraint = tfvs.PerimeterConstraint(1.0, self._init_perimeter_constraint)
        perim_constraint.thisown = 0
        tfvs.bind.surface(perim_constraint, _sh)

        return new_id

    def neighbor_surface_areas(self) -> Dict[str, Dict[str, float]]:
        result = {}
        for sh in self._cell_type:
            sh_id = self._cell_id_map_inv[sh.id]
            result[sh_id] = self._neighbor_surface_areas(sh)
        return result

    def num_cells(self) -> int:
        return len(self._cell_type)

    def cell_spatial_data(self):
        points = []
        cell_ids = []
        for sh in self._cell_type:
            points.append([v.position.xy().as_list() for v in sh.vertices])
            cell_ids.append(self._cell_id_map_inv[sh.id])
        return points, *tf.Universe.dim.xy().as_list(), cell_ids

    def cell_volumes(self) -> Dict[str, float]:
        result = {}
        for sh in self._cell_type:
            cell_id = self._cell_id_map_inv[sh.id]
            result[cell_id] = sh.area
        return result

    def set_cell_volume_targets(self, _targets: Dict[str, float]) -> None:
        for cell_id, cell_area in _targets.items():
            sh_id = self._cell_id_map[cell_id]
            sh = tfvs.SurfaceHandle(sh_id)
            sac: tfvs.SurfaceAreaConstraint = sh.surface_area_constraints[0]
            sac.constr = cell_area

    def divide_cells(self, _ids: Dict[str, Tuple[str, str]]) -> Dict[str, str]:
        result = {}
        for cell_id, (new_parent_id, child_id) in _ids.items():
            sh_id = self._cell_id_map[cell_id]
            sh = tfvs.SurfaceHandle(sh_id)
            ang = np.random.random() * 2 * np.pi
            new_sh = sh.split(sh.centroid, tf.FVector3(np.cos(ang), np.sin(ang), 0.0))
            self._process_new_cell(new_sh, new_id=child_id)
            self._set_id_map(sh, new_parent_id)
            result[self._cell_id_map_inv[sh.id]] = self._cell_id_map_inv[new_sh.id]
        return result

    # PySimService interface

    def _run(self) -> None:
        pass

    def _init(self) -> bool:

        # Initialize the simulation

        pad_x = 3 / 2 * self.cell_radius * 2
        pad_y = np.sqrt(3) * self.cell_radius * 2
        dim_x = 2 * pad_x + (self.num_cells_x - 1 / 2) * self.cell_radius * 3 / 2
        dim_y = 2 * pad_y + (2 * self.num_cells_y - 1 / 2) * self.cell_radius * np.sqrt(3) / 2

        min_dim = min(dim_x, dim_y)
        min_cells = 3
        len2cells = min_dim / min_cells

        dim = [dim_x, dim_y, min_dim]

        tf.init(
            windowless=not self._show,
            dim=dim,
            bc={'x': tf.BOUNDARY_FREESLIP, 'y': tf.BOUNDARY_FREESLIP, 'z': tf.BOUNDARY_FREESLIP},
            cells=[int(d / len2cells) + 1 for d in dim],
            dt=self._dt,
            cutoff=self.cell_radius)

        # Create a cell type

        self._cell_type = tfvs.SurfaceType()

        # Initialize the population

        start_pos = tf.FVector3([pad_x, pad_y, tf.Universe.center[2]])
        tfvs.create_hex2d_mesh(self._cell_type, start_pos, self.num_cells_x, self.num_cells_y, self.cell_radius)

        for sh in self._cell_type:
            self._process_new_cell(sh)

        # 2D simulation

        vtype = tfvs.MeshParticleType_get()
        vtype.frozen_z = True

        # Add some noise

        rforce = tf.Force.random(1E0, 0.0)
        rforce.thisown = 0
        tf.bind.force(rforce, vtype)

        if tf.err_occurred():
            print(tf.err_get_all())

        tf.show()

        return True

    def _start(self) -> bool:
        return True

    def _step(self) -> bool:
        return tf.step(self._step_size) == 0

    def _finish(self) -> None:
        pass


def test():
    sim = VertexPlanarSheet(10, 10, 1.0, show=True)
    # sim = VertexPlanarSheet(4, 9, 1.0, show=True)
    # sim = VertexPlanarSheet(9, 4, 1.0, show=True)
    sim._init()
    tf.system.camera_view_top()
    tf.show()


if __name__ == '__main__':
    test()
