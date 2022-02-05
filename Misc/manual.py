#!/usr/bin/env python

import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

# Can build a record manually or use translate_record() as below
features=[
    GraphicFeature(start=1000, end=1020, strand=+1, color="#ffd700",
                   label="Small feature"),
    GraphicFeature(start=1020, end=1500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),
    GraphicFeature(start=1400, end=1700, strand=-1, color="#cffccc",
                   label="Gene 2"),
    GraphicFeature(start=1600, end=1900, strand=+1, color="#ccccff",
                   label="Gene 3")
]
print(features)
record = GraphicRecord(sequence_length=2000, features=features)
ax1, _ = record.plot(draw_line = True, with_ruler = True)
plt.show()


