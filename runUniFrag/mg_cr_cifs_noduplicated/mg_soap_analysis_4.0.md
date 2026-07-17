# Mg-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Mg (`Mg`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **4.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9972`
* **Median SOAP Cosine Similarity**: `0.9984`
* **Average SOAP Fingerprint RMSD**: `0.0268`
* **Median SOAP Fingerprint RMSD**: `0.0268`
* **Total Mg Centers Analyzed**: `535`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 531 | **99.25%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 4 | **0.75%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 0 | **0.00%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Mg SOAP PCA and UMAP Map](mg_soap_distribution_4.0.png)

## Poorly Represented Mg Environments (Bottom 25 Worst Matches)

These Mg centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Mg Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `AVIPAX` | 0 | `0.9677` | `0.0974` | `AVIPAXFragMofMin` |
| 2 | `AVIPAX` | 1 | `0.9677` | `0.0974` | `AVIPAXFragMofMin` |
| 3 | `AVIPAX` | 2 | `0.9677` | `0.0974` | `AVIPAXFragMofMin` |
| 4 | `AVIPAX` | 3 | `0.9677` | `0.0974` | `AVIPAXFragMofMin` |
| 5 | `LIQHAX` | 0 | `0.9807` | `0.0667` | `LIQHAXFragMof` |
| 6 | `LIQHAX` | 2 | `0.9807` | `0.0667` | `LIQHAXFragMof` |
| 7 | `KAPRIG` | 0 | `0.9893` | `0.0583` | `EQERAUFragMofMin` |
| 8 | `KAPRIG` | 1 | `0.9893` | `0.0583` | `EQERAUFragMofMin` |
| 9 | `MOPQIT` | 6 | `0.9899` | `0.0767` | `AVIPAXFragMofMin` |
| 10 | `MOPQIT` | 5 | `0.9899` | `0.0767` | `AVIPAXFragMofMin` |
| 11 | `MOPQIT` | 7 | `0.9899` | `0.0767` | `AVIPAXFragMofMin` |
| 12 | `MOPQIT` | 4 | `0.9899` | `0.0767` | `AVIPAXFragMofMin` |
| 13 | `HIBGEF` | 3 | `0.9899` | `0.0849` | `AVIPAXFragMofMin` |
| 14 | `HIBGEF` | 2 | `0.9899` | `0.0849` | `AVIPAXFragMofMin` |
| 15 | `HIBGEF` | 0 | `0.9899` | `0.0849` | `AVIPAXFragMofMin` |
| 16 | `HIBGEF` | 1 | `0.9899` | `0.0849` | `AVIPAXFragMofMin` |
| 17 | `XEHSEJ` | 0 | `0.9900` | `0.0793` | `AVIPAXFragMofMin` |
| 18 | `XEHSEJ` | 1 | `0.9900` | `0.0793` | `AVIPAXFragMofMin` |
| 19 | `XEHSEJ` | 2 | `0.9900` | `0.0793` | `AVIPAXFragMofMin` |
| 20 | `XEHSEJ` | 3 | `0.9900` | `0.0793` | `AVIPAXFragMofMin` |
| 21 | `XEHSAF` | 0 | `0.9901` | `0.0830` | `AVIPAXFragMofMin` |
| 22 | `XEHSAF` | 1 | `0.9901` | `0.0830` | `AVIPAXFragMofMin` |
| 23 | `XEHSAF` | 2 | `0.9901` | `0.0830` | `AVIPAXFragMofMin` |
| 24 | `XEHSAF` | 3 | `0.9901` | `0.0830` | `AVIPAXFragMofMin` |
| 25 | `XEHSOT` | 0 | `0.9901` | `0.0843` | `AVIPAXFragMofMin` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high. This indicates that the local coordination environment of Mg (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 4.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent structures are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `4.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some metal centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of metal coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Mg environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
