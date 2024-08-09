cfg.txt  merged.gtf
GL000008.2      .       1000044-997804-997645-997433

GL000008.2      .       1000071-997804
GL000008.2      .       1000071-997804-997645-997433
GL000008.2      .       100009837-100003023-100002941-99972134-99971979-99969237-99969114-99965628
GL000008.2      .       100009837-99972134
GL000008.2      .       100009837-99972134-99971979-99969237
GL000008.2      .       100009837-99972134-99971979-99969237-99969114-99957205

Why do these have '.' strand and no tss start / tes start?
```bash
(pt_snakemake) bsc083001@glogin1:~/mele_lab/projects/240706_pantranscriptome_cerberus_gtf_merge/merge_only_ics/data/espresso/raw$
echo "" > espresso_coords.gtf

for f in *gtf
do
echo $f
grep 1000045 $f | head >> espresso_coords.gtf
done
```
