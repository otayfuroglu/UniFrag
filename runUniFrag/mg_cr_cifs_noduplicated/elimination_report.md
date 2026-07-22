# Fragment Size Elimination Report (Mg MOFs)

This report documents the post-processing step to filter out large molecular fragments with more than **200** atoms.

## Filter Parameters

- **Target Directory:** `runUniFrag/mg_cr_cifs_noduplicated`
- **Max Atoms Threshold:** `200`
- **Original ExtXYZ Backup:** `fragments_collection_original.extxyz`
- **Original CSV Backup:** `fragmentation_summary_original.csv`

## Statistics Summary

- **Total Frames Loaded:** 139
- **Frames Retained (<= 200 atoms):** 113
- **Frames Eliminated (> 200 atoms):** 26
- **Filtered Collection Path:** `fragments_collection.extxyz`
- **Eliminated Collection Path:** `eliminated_fragments.extxyz`

## Eliminated Fragments Detail

| Rank | Fragment Label | Parent REFCODE | Type | Atom Count | Formula | Elimination Reason |
| :---: | :--- | :---: | :---: | :---: | :--- | :--- |
| 1 | `DAJWETFragMof` | `DAJWET` | Normal | **556** | `C288 H192 Mg3 N24 O49` | CSV normal_atoms count (556) > 200 |
| 2 | `ORUKETFragMof` | `ORUKET` | Normal | **381** | `C168 H116 Mg N12 O72 P12` | CSV normal_atoms count (381) > 200 |
| 3 | `EQERAUFragMof` | `EQERAU` | Normal | **378** | `C168 H140 Mg4 O58 P8` | CSV normal_atoms count (378) > 200 |
| 4 | `WEHNAAFragMof` | `WEHNAA` | Normal | **355** | `C72 H174 Mg N12 O72 P24` | CSV normal_atoms count (355) > 200 |
| 5 | `OBIBAFFragMof` | `OBIBAF` | Normal | **324** | `C140 H136 Mg4 O44` | CSV normal_atoms count (324) > 200 |
| 6 | `KAPRIGFragMof` | `KAPRIG` | Normal | **323** | `C176 H112 Mg O34` | CSV normal_atoms count (323) > 200 |
| 7 | `EBIMOTFragMof` | `EBIMOT` | Normal | **315** | `C136 H112 Mg3 O64` | CSV normal_atoms count (315) > 200 |
| 8 | `MUDLONFragMof` | `MUDLON` | Normal | **309** | `C160 H108 Mg N8 O32` | CSV normal_atoms count (309) > 200 |
| 9 | `RAVWAOFragMof` | `RAVWAO` | Normal | **279** | `C132 H108 Mg3 O36` | CSV normal_atoms count (279) > 200 |
| 10 | `DIHHIPFragMof` | `DIHHIP` | Normal | **276** | `C132 H92 Mg4 O48` | CSV normal_atoms count (276) > 200 |
| 11 | `VALXOYFragMof` | `VALXOY` | Normal | **232** | `C96 H72 Mg3 N12 O49` | CSV normal_atoms count (232) > 200 |
| 12 | `AVIPAXFragMof` | `AVIPAX` | Normal | **220** | `C84 H72 Mg4 N12 O48` | CSV normal_atoms count (220) > 200 |
| 13 | `UDURUKFragMof` | `UDURUK` | Normal | **211** | `C96 H76 Mg3 O36` | CSV normal_atoms count (211) > 200 |
| 14 | `XUFYAAFragMof` | `XUFYAA` | Normal | **201** | `C108 H68 Mg O24` | CSV normal_atoms count (201) > 200 |
| 15 | `ORUKETFragMofOnlyLinker` | `ORUKET` | Normal | **102** | `C42 H36 N3 O18 P3` | CSV normal_atoms count (381) > 200 |
| 16 | `DAJWETFragMofOnlyLinker` | `DAJWET` | Normal | **94** | `C48 H34 N4 O8` | CSV normal_atoms count (556) > 200 |
| 17 | `KAPRIGFragMofOnlyLinker` | `KAPRIG` | Normal | **82** | `C44 H30 O8` | CSV normal_atoms count (323) > 200 |
| 18 | `MUDLONFragMofOnlyLinker` | `MUDLON` | Normal | **82** | `C40 H32 N2 O8` | CSV normal_atoms count (309) > 200 |
| 19 | `OBIBAFFragMofOnlyLinker` | `OBIBAF` | Normal | **74** | `C33 H32 O9` | CSV normal_atoms count (324) > 200 |
| 20 | `WEHNAAFragMofOnlyLinker` | `WEHNAA` | Normal | **66** | `C12 H36 N2 O12 P4` | CSV normal_atoms count (355) > 200 |
| 21 | `DIHHIPFragMofOnlyLinker` | `DIHHIP` | Normal | **48** | `C22 H18 O8` | CSV normal_atoms count (276) > 200 |
| 22 | `EQERAUFragMofOnlyLinker` | `EQERAU` | Normal | **48** | `C21 H19 O7 P` | CSV normal_atoms count (378) > 200 |
| 23 | `RAVWAOFragMofOnlyLinker` | `RAVWAO` | Normal | **48** | `C22 H20 O6` | CSV normal_atoms count (279) > 200 |
| 24 | `EBIMOTFragMofOnlyLinker` | `EBIMOT` | Normal | **41** | `C17 H16 O8` | CSV normal_atoms count (315) > 200 |
| 25 | `VALXOYFragMofOnlyLinker` | `VALXOY` | Normal | **40** | `C16 H14 N2 O8` | CSV normal_atoms count (232) > 200 |
| 26 | `AVIPAXFragMofOnlyLinker` | `AVIPAX` | Normal | **19** | `C7 H7 N O4` | CSV normal_atoms count (220) > 200 |
