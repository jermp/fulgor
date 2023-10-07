from sklearn.cluster import KMeans
import sys
import struct
import numpy as np

sketches_filename = sys.argv[1]
num_clusters = int(sys.argv[2])

sketches = []
with open(sketches_filename, 'rb') as f:
    num_bytes_per_sketch = struct.unpack('Q', f.read(8))[0]
    num_sketches = struct.unpack('Q', f.read(8))[0]
    for i in range(num_sketches):
        sketch = np.empty(num_bytes_per_sketch, dtype=np.int8)
        for j in range(num_bytes_per_sketch):
            sketch[j] = struct.unpack('B', f.read(1))[0]
        sketches.append(sketch)

# print(sketches)

kmeans = KMeans(
    n_clusters = num_clusters,
    random_state = 13,
    n_init = 'auto' # or 10, depending on the sklearn version
).fit(sketches)

for x in kmeans.labels_:
    print(x, end=' ')