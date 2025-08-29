from process_bigraph import ProcessTypes, Process, Composite, default
from process_bigraph.emitter import emitter_from_wires, gather_emitter_results

from multisim_matrix.vivarium import register_processes, register_types

# concentration = {
#     '_inherit': 'float',
#     '_divide': 'divide_any',
#     '_description': 'concentrations',
# }


class GrowDivide(Process):
    """
    """

    config_schema = {
        'threshold': 'float',
        'growth_rate': 'float',
        'cell_id': 'string',
        'cell_schema': 'schema',
        'divisions': {
            '_type': 'integer',
            '_default': 2
        }
    }

    # def initialize(self, config):
    #     """
    #     When this process is created,
    #     """
    #     pass

    def inputs(self):
        return {
            'activator': 'concentration',  # this determines growth rate (notch)
            'trigger': 'float'     # This triggers division when passed a threshold
        }

    def outputs(self):
        return {
            'target': 'float',             # this is the mass of the cell
            'environment': {
                '_type': 'map',
                '_value': self.config['cell_schema']},  # we need somewhere to divide into
        }

    def update(self, state, interval):

        activator_level = state['activator']
        trigger_level = state['trigger']

        divide = {}
        if trigger_level >= self.config['threshold']:
            # trigger division of self
            mother = self.config['cell_id']
            daughters = [(
                f'{mother}_{i}', {
                    'cell_id': f'{mother}_{i}'})
                for i in range(self.config['divisions'])]
            divide = {
                '_react': {
                    'divide': {
                        'mother': mother,
                        'daughters': daughters}}}

        # grow
        updated_target = self.config['growth_rate'] * activator_level * interval

        return {
            'target': updated_target,
            'environment': divide
        }

import pprint


def run_process(core):

    cell_schema = {
        'mass': 'cell_mass',
        'notch': 'concentration',
    }

    config = {
        'threshold': 2.0,   # this is the mass threshold for division
        'growth_rate': 0.1,
        'cell_id': '0',
        'cell_schema': cell_schema,
    }

    cell_id = '0'
    cell_state = {
        'mass': 1.0,
        'notch': 1.0,  # 1.0 concentration to trigger division
        'grow': {
            '_type': 'process',
            'address': 'local:gd_process',
            'config': config,
            'inputs': {
                'activator': ['notch'],
                'trigger': ['mass']
            },
            'outputs': {
                'target': ['mass'],
                'environment': ['..',]  # point IN the environment to the map of cells
            },
            'interval': 1.0
        },
    }

    environment_state = {
        'environment': {
            cell_id: cell_state
        },
        # add an emitter
        'emitter': emitter_from_wires({
            'time': ['global_time'],
            'environment': ['environment']
        })
    }

    # make the composite
    composite = Composite({
        'state': environment_state,
    }, core=core)

    duration = 20.0
    composite.run(duration)

    results = gather_emitter_results(composite)[('emitter',)]

    print(f'RESULTS: {pprint.pformat(results)}')


if __name__ == '__main__':
    core = ProcessTypes()
    register_types(core)
    register_processes(core)
    core.register_process(
        'gd_process',
        GrowDivide)

    run_process(core)