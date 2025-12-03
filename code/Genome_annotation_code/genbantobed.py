from Bio import SeqIO
import pandas as pd

data = []
for record in SeqIO.parse("f2.gb", "genbank"):
    for feature in record.features:
        if feature.type == "gene":
            row = {
                "molecule": record.id,
                "gene": feature.qualifiers.get("gene", [""])[0],
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "strand": "forward" if feature.location.strand > 0 else "reverse",
                "orientation": 1 if feature.location.strand > 0 else 0
            }
            data.append(row)

df = pd.DataFrame(data)
df.to_csv("output_genes.csv", index=False)
