import re

def parse_parameters_section(file_content):
    pattern = re.compile(r'&PARAMETERS\n(.*?)\n/', re.DOTALL)
    parameters_section = pattern.findall(file_content)[0]
    parameter_lines = parameters_section.split('\n')
    parameters = {}
    for line in parameter_lines:
        line = line.split('!')[0].strip()  # Ignore comments
        if '=' in line:
            key, value = line.split('=')
            key = key.upper().strip()  # Convert keys to uppercase
            value = value.strip()
            if '.' in value or 'e' in value.lower():  # The value is a float number
                parameters[key] = float(value)
            elif ',' in value:  # The value is a tuple of integers
                parameters[key] = tuple(int(i) for i in value.split(','))
            else:  # The value is an integer
                parameters[key] = int(value)
    return parameters

# Usage:
with open('wt.in', 'r') as f:
    file_content = f.read()

parameters = parse_parameters_section(file_content)
print(parameters)

