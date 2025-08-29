'''
Composites of Delta-Notch

TODO -- need to pass in ids for the cells so that we can synchronize.
'''
from process_bigraph import ProcessTypes, Composite, default, gather_emitter_results
from process_bigraph.emitter import emitter_from_wires
from bigraph_schema.registry import deep_merge_copy

from multisim_matrix import register_processes, register_types
import os


renderer_registry = {
    'CenterPlanarProcess': 'MCCenterRenderer2D',
    'PottsPlanarProcess': 'MCPottsRenderer2D',
    'VertexPlanarProcess': 'MCVertexRenderer2D'
}


def get_renderer_address(_address: str):
    domain, mcsim = _address.split(':')
    return f'{domain}:{renderer_registry[mcsim]}'


def render_output_dir(_root_dir: str, _mc_address: str, _sc_address: str):
    return os.path.join(_root_dir, f'{_mc_address.split(":")[1]}/{_sc_address.split(":")[1]}')


def run_composites(core):

    # TODO -- maka this work:
    # subcellular_processes = core.query('subcellular')  # TODO -- how do we get the list of possible subcellular processes?
    # multicellular_processes = core.query('multicellular')
    # TODO -- consider making subcellular simulators typical processes; startup for potentially 100s of services will be expensive without adding much value
    num_cells_x = 20 #4
    num_cells_y = 20 #4
    cell_radius = 5
    n_initial_cells = num_cells_x * num_cells_y
    step_size = 1.0  # time of one step in process bigraph engine
    dt = 0.5         # the period of time step according to tissue forge
    interval = 100.   # the total time to run the simulation
    assert step_size >= dt, 'The time step of the process bigraph engine must be greater than or equal to the time step of tissue forge'

    multicellular_startup_settings = {
        # 'local:PottsPlanarProcess': {},
        'local:CenterPlanarProcess': {
            'step_size': step_size,
            'dt': dt,
            'num_cells_x': 2,
            'num_cells_y': 2,
            'cell_radius': 5,
        },
        'local:VertexPlanarProcess': {
            'step_size': step_size,
            'dt': dt,
            'num_cells_x': 2,
            'num_cells_y': 2,
            'cell_radius': 5,
        },
    }
    subcellular_startup_settings = {
        'local:MaBoSSDeltaNotchProcess': {
            'time_step': dt * 1E0,
            'time_tick': dt * 1E0,
            'discrete_time': False,
            'seed': -1
        },
        'local:RoadRunnerDeltaNotchProcess': {
            'step_size': dt * 5,
            'num_steps': int(step_size/dt),
            'stochastic': False,
            'seed': 0
        },
    }

    # general config settings
    multicell_config = {
        'num_cells_x': num_cells_x,
        'num_cells_y': num_cells_y,
        'cell_radius': cell_radius
    }
    subcellular_config = {}
    grow_divide_config = {
        'threshold': 10.0,   # this is the volume threshold for division
        'growth_rate': 0.1,
    }

    fig_root_dir = '../_figs'
    print(os.path.abspath(fig_root_dir))

    # go through all the combinations of multicellular and subcellular processes
    for multicell_address, multicell_settings in multicellular_startup_settings.items():
        for subcell_address, subcell_settings in subcellular_startup_settings.items():

            # merge specific simulator settings with the general settings
            multicell_config_merged = deep_merge_copy(
                {'simservice_config': multicell_config}, multicell_settings)
            subcellular_config_merged = deep_merge_copy(
                subcellular_config, subcell_settings)

            # make the document
            document = {
                # 'neighborhood_surface_areas_store': {
                #     f'{n}': {
                #         f'{m}': 0.0 for m in range(n_initial_cells) if m != n
                #     } for n in range(n_initial_cells)
                # },
                'tissue': {
                    '_type': 'process',
                    'address': f'{multicell_address}',
                    'config': multicell_config_merged,
                    'inputs': {},
                    'outputs': {
                        'neighborhood_surface_areas': ['neighborhood_surface_areas'],
                        'cell_spatial_data': ['cell_spatial_data']
                        # this is a map from each cell to ids of its neighbors and their common surface area
                    }
                },
                'cells': {},
                'cell connector': {
                    '_type': 'step',
                    'address': 'local:CellConnector',
                    'config': {
                        'cells_count': num_cells_x * num_cells_y,
                        'read_molecules': ['delta'],  # TODO -- this will tell the connector what molecule id to read
                    },
                    'inputs': {
                        'connections': ['neighborhood_surface_areas'],  # this gives the connectivity and surface area
                        'cells': ['cells']  # it sees the cells so it can read their delta values
                    },
                    'outputs': {
                        'cells': ['cells']  # this updates the total delta values seen by each cell
                    }
                },
                'renderer': {
                    '_type': 'step',
                    'address': get_renderer_address(multicell_address),
                    'config': {
                        'render_specs': {
                            'dpi': 300,
                            'file_extensions': ['.png'],
                            'figure_height': 3.0,
                            'figure_width': 3.0
                        },
                        'output_dir': render_output_dir(fig_root_dir, multicell_address, subcell_address)
                    },
                    'inputs': {
                        'cells': ['cells'],
                        'cell_spatial_data': ['cell_spatial_data']
                    }
                },
                'emitter': emitter_from_wires({
                    'cells': ['cells'],
                    'neighborhood_surface_areas': ['neighborhood_surface_areas']}),
            }

            composition = {
                'cells': {
                    '_type': 'map',
                    '_value': {
                        'cell_process': {
                            '_type': 'process',
                            'address': default('string', f'{subcell_address}'),
                            'config': default('quote', subcellular_config_merged),
                            'inputs': default('tree[wires]', {
                                'delta_neighbors': ['delta_neighbors'],
                                'delta': ['delta'],
                                'notch': ['notch']
                            }),
                            'outputs': default('tree[wires]', {
                                'delta': ['delta'],  # this has to be called delta store for the connector to read it
                                'notch': ['notch']
                            })
                        },
                        'grow_divide_process': {
                            '_type': 'process',
                            'address': 'local:gd_process',
                            'config': {**grow_divide_config,
                                       **{'cell_id': 1,
                                          'cell_schema': {
                                              'volume': 'cell_volume',
                                              'notch': 'concentration',
                                          }}
                                       },
                            'inputs': {
                                'activator': ['notch'],
                                'trigger': ['volume']
                            },
                            'outputs': {
                                'target': ['volume'],
                                'environment': ['..', ]  # point IN the environment to the map of cells
                            },
                            'interval': 1.0
                        },
                    },

                }
            }

            # TODO -- set initial state

            # import ipdb; ipdb.set_trace()

            # make the composite
            print(f'Building composite with {multicell_address} and {subcell_address}')
            sim = Composite(
                config={
                    'state': document,
                    'composition': composition},
                core=core
            )

            # run the simulation
            print(f'Running composite with {multicell_address} and {subcell_address}')
            sim.run(interval=interval)

            # retrieve the results
            results = gather_emitter_results(sim)

            import ipdb; ipdb.set_trace()

            # print the results
            # TODO -- is the emitter not wired to the right location
            print(f'Results: {results[("emitter",)]}')


if __name__ == '__main__':
    core = ProcessTypes()
    register_types(core)

    run_composites(core)
