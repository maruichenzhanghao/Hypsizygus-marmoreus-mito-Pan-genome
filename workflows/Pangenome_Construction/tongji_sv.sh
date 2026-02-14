(base) maxinxin@maxinxin-X870E-AORUS-ELITE-WIFI7:~/workspace/01.mito/03.pangenome.cactus/mito_results815$ awk '$1=="S" {cnt++} END {print "Nodes (S lines):", cnt+0}' mito_pg.sv.gfa
Nodes (S lines): 217
(base) maxinxin@maxinxin-X870E-AORUS-ELITE-WIFI7:~/workspace/01.mito/03.pangenome.cactus/mito_results815$ awk '$1=="L" {cnt++} END {print "Edges (L lines):", cnt+0}' mito_pg.sv.gfa
Edges (L lines): 293
(base) maxinxin@maxinxin-X870E-AORUS-ELITE-WIFI7:~/workspace/01.mito/03.pangenome.cactus/mito_results815$ awk '$1=="P" {cnt++} END {print "Paths (P lines):", cnt+0}' mito_pg.sv.gfa
Paths (P lines): 0
(base) maxinxin@maxinxin-X870E-AORUS-ELITE-WIFI7:~/workspace/01.mito/03.pangenome.cactus/mito_results815$ grep -E "^S" mito_pg.sv.gfa | grep -E "LN:i:[0-9]+" | sed -E 's/.*LN:i:([0-9]+).*/\1/' | awk '{sum+=$1} END {print "Total length from LN:i:", sum+0}'
Total length from LN:i: 220364
