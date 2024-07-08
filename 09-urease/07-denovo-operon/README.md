# 1. Get seq ids for ureABC
```sh
cd ~/rna_dohmain/09-urease/07-denovo-operon
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureA | sed 's/\t.*//' > ureA.seqs
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureB | sed 's/\t.*//' > ureB.seqs
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureC | sed 's/\t.*//' > ureC.seqs
wc -l *seqs
  #  397 ureA.seqs
  #  314 ureB.seqs
  #  418 ureC.seqs
  # 1129 total
```
# 2. Find common seqids for ureABC
```sh
python3 common_elements.py
```
# 3. Find how many copies different sequences have
```sh
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureA.ids"> ureA.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureB.ids"> ureB.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureC.ids"> ureC.counts
paste -d "," ureABC.seqs ureA.counts ureB.counts ureC.counts >ureABC.counts