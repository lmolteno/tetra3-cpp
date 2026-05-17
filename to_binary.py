import numpy as np
import struct
import sys

def convert_npz_to_binary(npz_file, output_prefix):
    data = np.load(npz_file)
    # Save stars - ensure float32 and correct order
    stars = data['star_table'].astype(np.float32)
    print(stars[0])
    with open('tetra3_db_stars.bin', 'wb') as f:
        f.write(struct.pack('I', stars.shape[0]))
        f.write(stars.tobytes())
    
    # Save patterns - ensure uint16
    patterns = data['pattern_catalog']#.astype(np.uint16)  
    print(patterns[0])
    with open('tetra3_db_patterns.bin', 'wb') as f:
        f.write(struct.pack('I', patterns.shape[0]))
        f.write(patterns.tobytes())

    print(data['props_packed'])
# Usage
convert_npz_to_binary(sys.argv[1], 'tetra3_db')
