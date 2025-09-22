import json
import os

# DO NOT REMOVE THESE SEEMINGLY UNUSED IMPORTS!
import multisim_matrix.simservice.CenterPlanarSheetFactory
import multisim_matrix.simservice.VertexPlanarSheetFactory
import multisim_matrix.simservice.PottsPlanarSheetFactory
import multisim_matrix.simservice.MaBoSSDeltaNotchFactory
import multisim_matrix.simservice.RoadRunnerDeltaNotchFactory


from multisim_matrix.vivarium.CenterPlanarProcess import CenterPlanarProcess
from multisim_matrix.vivarium.MaBoSSDeltaNotchProcess import MaBoSSDeltaNotchProcess
from multisim_matrix.vivarium.PottsPlanarProcess import PottsPlanarProcess
from multisim_matrix.vivarium.RoadRunnerDeltaNotchProcess import RoadRunnerDeltaNotchProcess
from multisim_matrix.vivarium.VertexPlanarProcess import VertexPlanarProcess
from multisim_matrix.vivarium.cell_connector import CellConnector
from multisim_matrix.vivarium.MultiCellRenderer import MCCenterRenderer2D, MCPottsRenderer2D, MCVertexRenderer2D
from multisim_matrix.vivarium.grow_divide import GrowDivide

__processes__ = [
    CenterPlanarProcess,
    MaBoSSDeltaNotchProcess,
    PottsPlanarProcess,
    RoadRunnerDeltaNotchProcess,
    VertexPlanarProcess,
    CellConnector,
    MCCenterRenderer2D,
    MCPottsRenderer2D,
    MCVertexRenderer2D,
    GrowDivide
]


def register_processes(core):
    [core.register_process(p.__name__, p) for p in __processes__]
    return core


def register_types(core):
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schemas', 'type_info.json'), 'r') as f:
        for d in json.load(f):
            core.register(*d)

    register_processes(core)
    return core
