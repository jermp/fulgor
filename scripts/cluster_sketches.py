from sklearn.cluster import KMeans
import sys
import struct
import numpy as np

sketches_filename = sys.argv[1]
num_clusters = int(sys.argv[2])

sketches = []
set_sizes = []
with open(sketches_filename, 'rb') as f:
    num_bytes_per_sketch = struct.unpack('I', f.read(4))[0]
    num_sketches = struct.unpack('I', f.read(4))[0]
    for i in range(num_sketches):
        set_size = struct.unpack('I', f.read(4))[0]
        sketch = np.empty(num_bytes_per_sketch, dtype=np.int8)
        for j in range(num_bytes_per_sketch):
            sketch[j] = struct.unpack('B', f.read(1))[0]
        sketches.append(sketch)
        set_sizes.append(set_size)

# print(sketches)

kmeans = KMeans(
    n_clusters = num_clusters,
    random_state = 13,
    n_init='auto' # or 10, depending on the sklearn version
).fit(sketches)

for i in range(len(kmeans.labels_)):
    print(set_sizes[i], kmeans.labels_[i], end=' ')

# for x in kmeans.labels_:
#     print(x, end=' ')