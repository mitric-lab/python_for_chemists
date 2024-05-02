#!/usr/bin/env python

def parse_mol(mol_lines):
    symbols = []
    coords = []
    bonds = {}
    properties = {}
    for i, line in enumerate(mol_lines):
        if i == 0:
            # Title line
            title = line.strip()
        elif i == 1:
            # Program / file timestamp line
            stamp = line.strip()
        elif i == 2:
            # Comment line
            comment = line.strip()
        elif i == 3:
            # Counts line
            tmp = line.split()
            n_atoms = int(tmp[0])
            n_bonds = int(tmp[1])
        elif line.startswith('M'):
            tmp = line.split()
            if tmp[1] == 'END':
                # END line
                break
            else:
                # Property line
                properties[tmp[1]] = ' '.join(tmp[2:])
        else:
            # Atom or bond block
            tmp = line.split()
            try:
                int(tmp[0])
            except ValueError:
                # Atom block
                symbols.append(tmp[3].capitalize())
                coords.append([float(tmp[0]), float(tmp[1]), float(tmp[2])])
            else:
                # Bond block
                idx1 = int(tmp[0]) - 1
                idx2 = int(tmp[1]) - 1
                btype = int(tmp[2])
                btype = btype if btype <= 4 else 1

                bname = f'{symbols[idx1]}{symbols[idx2]}{btype}'
                if bname in bonds:
                    bonds[bname] += 1
                else:
                    bonds[bname] = 1
    
    mol = {
        'title': title,
        'stamp': stamp,
        'comment': comment,
        'symbols': symbols,
        'coords': coords,
        'bonds': bonds,
        'properties': properties,
    }
    return mol
            
