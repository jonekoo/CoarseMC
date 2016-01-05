#-*- coding: utf-8 -*-
import sys
import json

def main():
    filename = sys.argv[1]
    parameter_ds = convert_parameters(sys.argv[2])
    
    f = open(filename, 'r')
    i_configuration = 0

    rodgroup = {}
    rodgroup['type'] = 'rod'
    rodgroup['name'] = 'gb'
    rodgroup['description'] = ['x', 'y', 'z', 'ux', 'uy', 'uz']

    pointgroup = {}
    pointgroup['type'] = 'point'
    pointgroup['name'] = 'lj'
    pointgroup['description'] = ['x', 'y', 'z']
    #print '[',
    json_ds = []
    for line in f:
        if line.startswith('#configuration'):
            json_d = {}
            rodgroup['coordinates'] = []
            pointgroup['coordinates'] = []
            groups = ['gb']
        elif line.startswith('#end_configuration'):
            particle_groups = [rodgroup]
            rodgroup['n'] = len(rodgroup['coordinates'])
            if (len(pointgroup['coordinates']) > 0):
                pointgroup['n'] = len(pointgroup['coordinates'])
                particle_groups.append(pointgroup)
                groups.append('lj')
            json_d['groups'] = groups
            json_d['particle_groups'] = particle_groups
            json_d.update(parameter_ds[i_configuration])
            #if i_configuration > 0:
            #    print ','
            #print json.dumps(json_d, sort_keys=False),
            json_ds.append(json_d.copy())
            i_configuration += 1
        elif line.startswith('gb'):
            rodgroup['coordinates'].append(
                [float(item) for item in line.split()[1:]])
        elif line.startswith('lj'):
            pointgroup['coordinates'].append(
                [float(item) for item in line.split()[1:4]])
        elif line.startswith('cylindrical'):
            records = line.split()
            json_d['box'] = {'type': 'cylindrical',
                             'diameter': float(records[1]),
                             'length': float(records[3]),
                             'periodic': True}
        elif line.startswith('rectangular'):
            records = line.split()
            json_d['box'] = {'type': 'rectangular',
                             'lx': float(records[1]),
                             'ly': float(records[2]),
                             'lz': float(records[3]),
                             'xperiodic': True,
                             'yperiodic': True,
                             'zperiodic': True}
    #print ']'
    if len(json_ds) > 1:
        print json.dumps(json_ds, sort_keys=False)
    elif len(json_ds) == 1:
        print json.dumps(json_ds[0], sort_keys=False)
        

def convert_parameters(parameterfile):
    # Convert parameters:
    f = open(parameterfile)
    parameter_ds = []
    
    parameter_d = None
    for line in f:
        if 'mc engine parameters' in line:
            if parameter_d is not None:
                edit_parameters(parameter_d)
                parameter_ds.append(parameter_d)
            parameter_d = {}
        elif line.startswith('$'):
            records = line.split()
            if records[1] == 'T':
                val = True
            elif records[1] == 'F':
                val = False
            elif records[1] == 'NaN':
                val = None
            else:
                # convert to float, int etc.
                val = json.loads(records[1])
            parameter_d[records[0][1:]] = val
    edit_parameters(parameter_d)
    parameter_ds.append(parameter_d)
    return parameter_ds

def edit_parameters(parameter_d):
    if 'scaling_type' in parameter_d:
        parameter_d['scaling_types'] = parameter_d['scaling_type'].split()
        del parameter_d['scaling_type']
        parameter_d['ensemble'] = 'npt'
    else:
        parameter_d['ensemble'] = 'nvt'
    if 'pt_ratio' in parameter_d:
        parameter_d['pt_period'] = parameter_d['pt_ratio']
        del parameter_d['pt_ratio']
    edit_pair_interactions(parameter_d)
    edit_particlewall(parameter_d)


def edit_pair_interactions(parameter_d):
    # Construct pair interactions:
    parameter_d["pair_interactions"] = []
    gb_interaction = {}
    gblj_interaction = {}
    lj_interaction = {}
    for key in parameter_d.keys():
        if key.startswith('gb_'):
            gb_interaction[key] = parameter_d[key]
            del parameter_d[key]
        elif key.startswith('gblj_'):
            gblj_interaction[key] = parameter_d[key]
            del parameter_d[key]
        elif key.startswith('lj_'):  
            lj_interaction[key] = parameter_d[key]
            del parameter_d[key]
    if len(gb_interaction.keys()) > 0:
        gb_interaction['type'] = 'gayberne'
        gb_interaction['participants'] = ['gb', 'gb']
        gb_interaction['r_cutoff'] = parameter_d['r_cutoff']
        gb_interaction['gb_hardcore'] = 0.6
        parameter_d['pair_interactions'].append(gb_interaction)
    if len(gblj_interaction.keys()) > 0:
        gblj_interaction['type'] = 'gblj'
        gblj_interaction['participants'] = ['gb', 'lj']
        gblj_interaction['r_cutoff'] = parameter_d['r_cutoff']
        gblj_interaction['gblj_hardcore'] = 0.3
        parameter_d['pair_interactions'].append(gblj_interaction)
    if len(lj_interaction.keys()) > 0:
        lj_interaction['type'] = 'lj'
        lj_interaction['participants'] = ['lj', 'lj']
        lj_interaction['r_cutoff'] = parameter_d['r_cutoff']
        parameter_d['pair_interactions'].append(lj_interaction)


def edit_particlewall(parameter_d):
    #Construct wall interaction if any:
    particlewall_keys = ['sigwall_LJ', 'Kw', 'Kw_LJ', 'sigwall',
                         'alpha_A', 'alpha_B', 'alpha_LJ', 'LJ_dist',
                         'wall_density', 'is_uniform_alignment', 'epswall',
                         'epswall_LJ']
    if parameter_d['is_wall_on']:
        single_interactions = []
        particlewall_d = {}
        for key in particlewall_keys:
            if key in parameter_d:
                particlewall_d[key] = parameter_d[key]
        if 'Kw' in particlewall_d:
            particlewall_d['epswall'] = particlewall_d['Kw'] / 8.0
        if 'Kw_LJ' in particlewall_d:
            particlewall_d['epswall_LJ'] = particlewall_d['Kw_LJ'] / \
            (particlewall_d['Kw'] * (particlewall_d['sigwall_LJ'] / \
                                     particlewall_d['sigwall'])**3)
        temp = particlewall_d.pop('Kw', None)
        temp = particlewall_d.pop('Kw_LJ', None)
        particlewall_d['type'] = 'ljwall_interaction'
        particlewall_d['rodwall'] = 'ljdimer-wall'
        particlewall_d['pointwall'] = 'ljcylinder'
        particlewall_d['wall_density'] = 1.0
        parameter_d['single_interactions'] = [particlewall_d, 
                                              particlewall_d.copy()]
        parameter_d['single_interactions'][0]['participant'] = 'gb'
        parameter_d['single_interactions'][1]['participant'] = 'lj'
            
    for key in particlewall_keys:
        temp = parameter_d.pop(key, None)
    temp = parameter_d.pop('is_wall_on')
        
    #return parameter_ds   
    #print json.dumps(parameter_d, sort_keys=False, indent=4,
    #                 separators=(',', ': '))

    
                     
if __name__=='__main__':
    main()
