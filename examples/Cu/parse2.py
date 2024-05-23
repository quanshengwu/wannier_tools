import re

def parse_fortran_namelist(file_content):
    sections = re.split(r'&(\w+)|(\w+)', file_content)
    parameters = {}
    for i in range(1, len(sections), 3):
        section_name = sections[i] or sections[i+1]
        section_content = sections[i+2]
        if section_name == 'LATTICE':
            parameters['LATTICE'] = parse_lattice_section(section_content)
        elif section_name == 'ATOM_POSITIONS':
            parameters['ATOM_POSITIONS'] = parse_atom_positions_section(section_content)
        elif section_name == 'PROJECTORS':
            parameters['PROJECTORS'] = parse_projectors_section(section_content)
        elif section_name == 'SURFACE':
            parameters['SURFACE'] = parse_surface_section(section_content)
        elif section_name == 'SELECTEDBANDS':
            parameters['SELECTEDBANDS'] = parse_selectedbands_section(section_content)
        elif section_name == 'KCUBE_BULK':
            parameters['KCUBE_BULK'] = parse_kcube_bulk_section(section_content)
        else:
            parameters[section_name] = parse_key_value_section(section_content)
    return parameters

def parse_key_value_section(section):
    # Remove comments
    section = re.sub(r'!.*', '', section)
    # Parse parameters
    pattern = r"(\w+)\s*=\s*([^\s/]+)"
    matches = re.findall(pattern, section)
    section_parameters = {key: parse_value(value) for key, value in matches}
    return section_parameters

def parse_lattice_section(section):
    lines = section.split('\n')[1:]  # Skip the 'Angstrom' line
    matrix = [list(map(float, line.split())) for line in lines if line.strip()]
    return matrix

def parse_atom_positions_section(section):
    lines = [line for line in section.split('\n') if line.strip()]  # Filter out empty lines
    number_of_atoms = int(lines[0].split()[1])  # Number of atoms is on the first line
    coordinate_type = lines[1].strip()  # Coordinate type is on the second line
    atom_positions = [line.split() for line in lines[2:2+number_of_atoms]]  # Atom positions start from the third line
    return {'number_of_atoms': number_of_atoms, 'coordinate_type': coordinate_type, 'atom_positions': atom_positions}


def parse_projectors_section(section):
    lines = section.split('\n')[1:]  # Skip the first line
    number_of_projectors = int(lines[0].strip())
    projectors = [line.split() for line in lines[1:number_of_projectors+1] if line.strip()]
    return {'number_of_projectors': number_of_projectors, 'projectors': projectors}

def parse_surface_section(section):
    lines = section.split('\n')[1:]  # Skip the first line
    matrix = [list(map(int, line.split())) for line in lines if line.strip()]
    return matrix

def parse_selectedbands_section(section):
    lines = section.split('\n')[1:]  # Skip the first line
    number_of_bands = int(lines[0].strip())
    bands = [int(line.strip()) for line in lines[1:number_of_bands+1] if line.strip()]
    return {'number_of_bands': number_of_bands, 'bands': bands}

def parse_kcube_bulk_section(section):
    lines = section.split('\n')[1:]  # Skip the first line
    vectors = [list(map(float, line.split())) for line in lines if line.strip()]
    return vectors

def parse_value(value):
    if ',' in value:  # Parse tuples like '0, 90'
        return tuple(map(float, value.split(',')))
    elif '.' in value or 'e' in value.lower():
        return float(value)
    elif value.lower() in ['t', 'f']:
        return value.lower() == 't'
    else:
        return int(value)

# Open and read the file content
with open('wt.in', 'r') as file:
    file_content = file.read()

parameters = parse_fortran_namelist(file_content)
print(parameters)


