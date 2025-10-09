#!/usr/bin/env python3
import re

# Read the file
with open('RodneyHQMOM.jl/src/simulation_runner.jl', 'r') as f:
    lines = f.readlines()

# Remove debug blocks
output_lines = []
skip_until_end = 0
i = 0
while i < len(lines):
    line = lines[i]
    
    # Check if this is a DEBUG comment line
    if '# DEBUG:' in line or '⚠️' in line or 'ALREADY CORRUPTED' in line:
        # Skip this line and look for the matching if/end block
        if i+1 < len(lines) and 'if rank == 0' in lines[i+1]:
            # Skip until we find the matching end
            skip_until_end = 1
            i += 1
            while i < len(lines) and skip_until_end > 0:
                if lines[i].strip().startswith('if ') or lines[i].strip().startswith('for '):
                    skip_until_end += 1
                elif lines[i].strip() == 'end':
                    skip_until_end -= 1
                i += 1
            continue
    else:
        output_lines.append(line)
    
    i += 1

# Write the cleaned file
with open('RodneyHQMOM.jl/src/simulation_runner.jl', 'w') as f:
    f.writelines(output_lines)

print("Debug statements removed successfully!")

