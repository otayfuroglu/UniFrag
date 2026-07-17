# Mg-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Mg (`Mg`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **6.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9903`
* **Median SOAP Cosine Similarity**: `0.9927`
* **Average SOAP Fingerprint RMSD**: `0.0988`
* **Median SOAP Fingerprint RMSD**: `0.0920`
* **Total Mg Centers Analyzed**: `535`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 468 | **87.48%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 67 | **12.52%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 0 | **0.00%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Mg SOAP PCA and UMAP Map](mg_soap_distribution_6.0.png)

## Poorly Represented Mg Environments (Bottom 25 Worst Matches)

These Mg centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Mg Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `BAKYOE` | 4 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 2 | `BAKYOE` | 6 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 3 | `BAKYOE` | 5 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 4 | `BAKYOE` | 7 | `0.9529` | `0.2277` | `HIBGEFFragMof` |
| 5 | `QIWPET` | 4 | `0.9663` | `0.3083` | `EQERAUFragMofMin` |
| 6 | `WAMRIN` | 1 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 7 | `WAMRIN` | 2 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 8 | `WAMRIN` | 0 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 9 | `WAMRIN` | 3 | `0.9689` | `0.1787` | `NUDMOQFragMofMin` |
| 10 | `LIQHAX` | 0 | `0.9700` | `0.1376` | `LIQHAXFragMof` |
| 11 | `LIQHAX` | 2 | `0.9700` | `0.1376` | `LIQHAXFragMof` |
| 12 | `HIBGEF` | 0 | `0.9714` | `0.3227` | `DAFYANFragMof` |
| 13 | `HIBGEF` | 1 | `0.9714` | `0.3227` | `DAFYANFragMof` |
| 14 | `HIBGEF` | 2 | `0.9714` | `0.3227` | `DAFYANFragMof` |
| 15 | `HIBGEF` | 3 | `0.9714` | `0.3227` | `DAFYANFragMof` |
| 16 | `XEHSOT` | 0 | `0.9715` | `0.3221` | `DAFYANFragMof` |
| 17 | `XEHSOT` | 2 | `0.9715` | `0.3221` | `DAFYANFragMof` |
| 18 | `XEHSOT` | 1 | `0.9715` | `0.3221` | `DAFYANFragMof` |
| 19 | `XEHSOT` | 3 | `0.9715` | `0.3221` | `DAFYANFragMof` |
| 20 | `XEHSIN` | 3 | `0.9716` | `0.3233` | `DAFYANFragMof` |
| 21 | `XEHSIN` | 0 | `0.9716` | `0.3233` | `DAFYANFragMof` |
| 22 | `XEHSIN` | 1 | `0.9716` | `0.3233` | `DAFYANFragMof` |
| 23 | `XEHSIN` | 2 | `0.9716` | `0.3233` | `DAFYANFragMof` |
| 24 | `XEHRUY` | 1 | `0.9717` | `0.3204` | `DAFYANFragMof` |
| 25 | `XEHRUY` | 2 | `0.9717` | `0.3204` | `DAFYANFragMof` |

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
