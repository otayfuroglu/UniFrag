# Fragment Size Elimination Report (Mg MOFs)

This report documents the post-processing step to filter out large molecular fragments with more than **200** atoms.

## Filter Parameters

- **Target Directory:** `runUniFrag/mg_cr_cifs_noduplicated`
- **Max Atoms Threshold:** `200`
- **Original ExtXYZ Backup:** `fragments_collection_original.extxyz`
- **Original CSV Backup:** `fragmentation_summary_original.csv`

## Statistics Summary

- **Total Frames Loaded:** 101
- **Frames Retained (<= 200 atoms):** 89
- **Frames Eliminated (> 200 atoms):** 12
- **Filtered Collection Path:** `fragments_collection.extxyz`
- **Eliminated Collection Path:** `eliminated_fragments.extxyz`

## Eliminated Fragments Detail

| Rank | Fragment Label | Parent REFCODE | Type | Atom Count | Formula | Elimination Reason |
| :---: | :--- | :---: | :---: | :---: | :--- | :--- |
| 1 | `DAJWETFragMof` | `DAJWET` | Normal | **556** | `C288 H192 Mg3 N24 O49` | CSV normal_atoms count (556) > 200 |
| 2 | `ORUKETFragMof` | `ORUKET` | Normal | **381** | `C168 H116 Mg N12 O72 P12` | CSV normal_atoms count (381) > 200 |
| 3 | `EQERAUFragMof` | `EQERAU` | Normal | **378** | `C168 H140 Mg4 O58 P8` | CSV normal_atoms count (378) > 200 |
| 4 | `OBIBAFFragMof` | `OBIBAF` | Normal | **326** | `C140 H138 Mg4 O44` | CSV normal_atoms count (326) > 200 |
| 5 | `EBIMOTFragMof` | `EBIMOT` | Normal | **325** | `C136 H122 Mg3 O64` | CSV normal_atoms count (325) > 200 |
| 6 | `KAPRIGFragMof` | `KAPRIG` | Normal | **315** | `C172 H108 Mg O34` | CSV normal_atoms count (315) > 200 |
| 7 | `MUDLONFragMof` | `MUDLON` | Normal | **309** | `C160 H108 Mg N8 O32` | CSV normal_atoms count (309) > 200 |
| 8 | `DIHHIPFragMof` | `DIHHIP` | Normal | **284** | `C132 H100 Mg4 O48` | CSV normal_atoms count (284) > 200 |
| 9 | `RAVWAOFragMof` | `RAVWAO` | Normal | **279** | `C132 H108 Mg3 O36` | CSV normal_atoms count (279) > 200 |
| 10 | `VALXOYFragMof` | `VALXOY` | Normal | **232** | `C96 H72 Mg3 N12 O49` | CSV normal_atoms count (232) > 200 |
| 11 | `AVIPAXFragMof` | `AVIPAX` | Normal | **222** | `C84 H74 Mg4 N12 O48` | CSV normal_atoms count (222) > 200 |
| 12 | `UDUSARFragMof` | `UDUSAR` | Normal | **211** | `C96 H76 Mg3 O36` | CSV normal_atoms count (211) > 200 |
