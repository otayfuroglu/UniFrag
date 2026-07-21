# Mg-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Mg (`Mg`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **6.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9890`
* **Median SOAP Cosine Similarity**: `0.9940`
* **Average SOAP Fingerprint RMSD**: `0.0950`
* **Median SOAP Fingerprint RMSD**: `0.0918`
* **Total Mg Centers Analyzed**: `551`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 464 | **84.21%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 83 | **15.06%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 4 | **0.73%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Mg SOAP PCA and UMAP Map](mg_soap_distribution_6.0.png)

## Poorly Represented Mg Environments (Bottom 25 Worst Matches)

These Mg centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Mg Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `BAXSIE` | 15 | `0.8120` | `0.2927` | `RATVEPFragMofMin` |
| 2 | `BAXSIE` | 14 | `0.8120` | `0.2927` | `RATVEPFragMofMin` |
| 3 | `BAXSIE` | 12 | `0.8120` | `0.2927` | `RATVEPFragMofMin` |
| 4 | `BAXSIE` | 13 | `0.8120` | `0.2927` | `RATVEPFragMofMin` |
| 5 | `BAKYOE` | 4 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 6 | `BAKYOE` | 6 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 7 | `BAKYOE` | 5 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 8 | `BAKYOE` | 7 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 9 | `BAXSIE` | 10 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 10 | `BAXSIE` | 4 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 11 | `BAXSIE` | 5 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 12 | `BAXSIE` | 11 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 13 | `BAXSIE` | 8 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 14 | `BAXSIE` | 7 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 15 | `BAXSIE` | 6 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 16 | `BAXSIE` | 9 | `0.9623` | `0.1822` | `BAXSIEFragMof` |
| 17 | `QIWPET` | 4 | `0.9663` | `0.3083` | `EQERAUFragMofMin` |
| 18 | `BAXSIE` | 0 | `0.9673` | `0.2298` | `EQERAUFragMofMin` |
| 19 | `BAXSIE` | 1 | `0.9673` | `0.2298` | `EQERAUFragMofMin` |
| 20 | `BAXSIE` | 2 | `0.9673` | `0.2298` | `EQERAUFragMofMin` |
| 21 | `BAXSIE` | 3 | `0.9673` | `0.2298` | `EQERAUFragMofMin` |
| 22 | `WAMRIN` | 1 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 23 | `WAMRIN` | 2 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 24 | `WAMRIN` | 0 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 25 | `WAMRIN` | 3 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high. This indicates that the local coordination environment of Mg (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 6.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent structures are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `6.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some metal centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of metal coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Mg environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
