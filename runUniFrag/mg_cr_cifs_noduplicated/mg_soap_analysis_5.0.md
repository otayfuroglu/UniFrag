# Mg-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Mg (`Mg`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **5.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9921`
* **Median SOAP Cosine Similarity**: `0.9955`
* **Average SOAP Fingerprint RMSD**: `0.0665`
* **Median SOAP Fingerprint RMSD**: `0.0584`
* **Total Mg Centers Analyzed**: `535`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 460 | **85.98%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 75 | **14.02%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 0 | **0.00%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Mg SOAP PCA and UMAP Map](mg_soap_distribution_5.0.png)

## Poorly Represented Mg Environments (Bottom 25 Worst Matches)

These Mg centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Mg Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `WEHNAA` | 1 | `0.9584` | `0.1520` | `WEHNAAFragMofMin` |
| 2 | `WEHNAA` | 0 | `0.9584` | `0.1520` | `WEHNAAFragMofMin` |
| 3 | `LIQHAX` | 0 | `0.9677` | `0.1261` | `NUDMUWFragMofMin` |
| 4 | `LIQHAX` | 2 | `0.9677` | `0.1261` | `NUDMUWFragMofMin` |
| 5 | `NUDLIJ` | 8 | `0.9713` | `0.1184` | `NUDMUWFragMof` |
| 6 | `NUDLIJ` | 10 | `0.9713` | `0.1184` | `NUDMUWFragMof` |
| 7 | `NUDLIJ` | 11 | `0.9713` | `0.1184` | `NUDMUWFragMof` |
| 8 | `NUDLIJ` | 9 | `0.9713` | `0.1184` | `NUDMUWFragMof` |
| 9 | `NUDLIJ` | 1 | `0.9719` | `0.1176` | `NUDMUWFragMof` |
| 10 | `NUDLIJ` | 2 | `0.9719` | `0.1176` | `NUDMUWFragMof` |
| 11 | `NUDLIJ` | 3 | `0.9719` | `0.1176` | `NUDMUWFragMof` |
| 12 | `NUDLIJ` | 0 | `0.9719` | `0.1176` | `NUDMUWFragMof` |
| 13 | `BAKYOE` | 4 | `0.9722` | `0.1520` | `BAKYOEFragMof` |
| 14 | `BAKYOE` | 6 | `0.9722` | `0.1520` | `BAKYOEFragMof` |
| 15 | `BAKYOE` | 7 | `0.9722` | `0.1520` | `BAKYOEFragMof` |
| 16 | `BAKYOE` | 5 | `0.9722` | `0.1520` | `BAKYOEFragMof` |
| 17 | `LIQHAX` | 1 | `0.9723` | `0.1211` | `NUDMUWFragMofMin` |
| 18 | `LIQHAX` | 3 | `0.9723` | `0.1211` | `NUDMUWFragMofMin` |
| 19 | `NUDMUW` | 3 | `0.9735` | `0.1136` | `NUDMUWFragMof` |
| 20 | `NUDMUW` | 2 | `0.9735` | `0.1136` | `NUDMUWFragMof` |
| 21 | `NUDMUW` | 4 | `0.9735` | `0.1136` | `NUDMUWFragMof` |
| 22 | `NUDMUW` | 5 | `0.9735` | `0.1136` | `NUDMUWFragMof` |
| 23 | `NUDNEH` | 0 | `0.9740` | `0.1125` | `NUDMUWFragMof` |
| 24 | `NUDNEH` | 2 | `0.9740` | `0.1125` | `NUDMUWFragMof` |
| 25 | `NUDNEH` | 3 | `0.9740` | `0.1125` | `NUDMUWFragMof` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high. This indicates that the local coordination environment of Mg (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 5.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent structures are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `5.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some metal centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of metal coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Mg environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
