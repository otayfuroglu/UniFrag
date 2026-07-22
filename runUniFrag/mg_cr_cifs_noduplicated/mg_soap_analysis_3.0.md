# Mg-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Mg (`Mg`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **3.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9992`
* **Median SOAP Cosine Similarity**: `0.9997`
* **Average SOAP Fingerprint RMSD**: `0.0079`
* **Median SOAP Fingerprint RMSD**: `0.0069`
* **Total Mg Centers Analyzed**: `535`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 531 | **99.25%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 4 | **0.75%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 0 | **0.00%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Mg SOAP PCA and UMAP Map](mg_soap_distribution_3.0.png)

## Poorly Represented Mg Environments (Bottom 25 Worst Matches)

These Mg centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Mg Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `AVIPAX` | 2 | `0.9769` | `0.0380` | `LIQHAXFragMof` |
| 2 | `AVIPAX` | 1 | `0.9769` | `0.0380` | `LIQHAXFragMof` |
| 3 | `AVIPAX` | 0 | `0.9769` | `0.0380` | `LIQHAXFragMof` |
| 4 | `AVIPAX` | 3 | `0.9769` | `0.0380` | `LIQHAXFragMof` |
| 5 | `LIQHAX` | 0 | `0.9820` | `0.0572` | `NUDMUWFragMof` |
| 6 | `LIQHAX` | 2 | `0.9820` | `0.0572` | `NUDMUWFragMof` |
| 7 | `LIQHAX` | 1 | `0.9856` | `0.0381` | `XESKAJFragMofMin` |
| 8 | `LIQHAX` | 3 | `0.9856` | `0.0381` | `XESKAJFragMofMin` |
| 9 | `BAKYOE` | 0 | `0.9956` | `0.0288` | `DUWRAQFragMof` |
| 10 | `BAKYOE` | 1 | `0.9956` | `0.0288` | `DUWRAQFragMof` |
| 11 | `BAKYOE` | 2 | `0.9956` | `0.0288` | `DUWRAQFragMof` |
| 12 | `BAKYOE` | 3 | `0.9956` | `0.0288` | `DUWRAQFragMof` |
| 13 | `KAPRIG` | 1 | `0.9963` | `0.0242` | `EQERAUFragMofMin` |
| 14 | `KAPRIG` | 0 | `0.9963` | `0.0242` | `EQERAUFragMofMin` |
| 15 | `TAGVAB` | 13 | `0.9977` | `0.0161` | `TAGVABFragMof` |
| 16 | `TAGVAB` | 15 | `0.9977` | `0.0161` | `TAGVABFragMof` |
| 17 | `TAGVAB` | 14 | `0.9977` | `0.0161` | `TAGVABFragMof` |
| 18 | `TAGVAB` | 16 | `0.9977` | `0.0161` | `TAGVABFragMof` |
| 19 | `TAGVAB` | 12 | `0.9977` | `0.0161` | `TAGVABFragMof` |
| 20 | `TAGVAB` | 17 | `0.9977` | `0.0161` | `TAGVABFragMof` |
| 21 | `HAFVUH` | 1 | `0.9978` | `0.0220` | `HAFVUHFragMofMin` |
| 22 | `HAFVUH` | 2 | `0.9978` | `0.0220` | `HAFVUHFragMofMin` |
| 23 | `HAFVUH` | 4 | `0.9978` | `0.0220` | `HAFVUHFragMofMin` |
| 24 | `HAFVUH` | 0 | `0.9978` | `0.0220` | `HAFVUHFragMofMin` |
| 25 | `HAFVUH` | 5 | `0.9978` | `0.0220` | `HAFVUHFragMofMin` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high. This indicates that the local coordination environment of Mg (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 3.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent structures are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `3.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some metal centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of metal coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Mg environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
