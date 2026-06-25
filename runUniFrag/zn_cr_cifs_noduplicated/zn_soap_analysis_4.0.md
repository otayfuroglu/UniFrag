# Zn-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Zinc (`Zn`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **4.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9976`
* **Median SOAP Cosine Similarity**: `0.9995`
* **Average SOAP Fingerprint RMSD**: `0.0120`
* **Median SOAP Fingerprint RMSD**: `0.0084`
* **Total Zn Centers Analyzed**: `9176`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 9030 | **98.41%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 134 | **1.46%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 12 | **0.13%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Zn SOAP PCA and UMAP Map](file:///Users/omert/.gemini/antigravity/brain/153d2da2-7e4a-4474-b36a-6be8db573d0d/zn_soap_distribution_4.0.png)

## Poorly Represented Zn Environments (Bottom 25 Worst Matches)

These Zn centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Zn Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `REDYEG` | 0 | `0.8253` | `0.1392` | `REDYEGFragMof` |
| 2 | `REDYEG` | 1 | `0.8253` | `0.1392` | `REDYEGFragMof` |
| 3 | `OQOXAV` | 3 | `0.8328` | `0.1680` | `OQOXAVFragMof` |
| 4 | `OQOXAV` | 4 | `0.8328` | `0.1680` | `OQOXAVFragMof` |
| 5 | `OQOXAV` | 2 | `0.8328` | `0.1680` | `OQOXAVFragMof` |
| 6 | `OQOXAV` | 5 | `0.8328` | `0.1680` | `OQOXAVFragMof` |
| 7 | `RUMRUO` | 1 | `0.8331` | `0.1675` | `OQOXAVFragMof` |
| 8 | `RUMRUO` | 0 | `0.8331` | `0.1675` | `OQOXAVFragMof` |
| 9 | `REDYEG` | 2 | `0.8813` | `0.1106` | `REDYEGFragMof` |
| 10 | `REDYEG` | 4 | `0.8813` | `0.1106` | `REDYEGFragMof` |
| 11 | `REDYEG` | 3 | `0.8813` | `0.1106` | `REDYEGFragMof` |
| 12 | `REDYEG` | 5 | `0.8813` | `0.1106` | `REDYEGFragMof` |
| 13 | `MOMRIR` | 8 | `0.9110` | `0.1388` | `CACBUGFragMofMin` |
| 14 | `MOMRIR` | 6 | `0.9110` | `0.1388` | `CACBUGFragMofMin` |
| 15 | `MOMRIR` | 7 | `0.9110` | `0.1388` | `CACBUGFragMofMin` |
| 16 | `MOMRIR` | 9 | `0.9110` | `0.1388` | `CACBUGFragMofMin` |
| 17 | `MOMRIR` | 14 | `0.9131` | `0.1387` | `CACBUGFragMofMin` |
| 18 | `MOMRIR` | 15 | `0.9131` | `0.1387` | `CACBUGFragMofMin` |
| 19 | `WAWGOQ` | 1 | `0.9299` | `0.0955` | `MUDTAJFragMofMin` |
| 20 | `WAWGOQ` | 5 | `0.9299` | `0.0955` | `MUDTAJFragMofMin` |
| 21 | `WAWGOQ` | 0 | `0.9299` | `0.0955` | `MUDTAJFragMofMin` |
| 22 | `WAWGOQ` | 4 | `0.9299` | `0.0955` | `MUDTAJFragMofMin` |
| 23 | `XATZOK` | 0 | `0.9301` | `0.0958` | `MUDTAJFragMofMin` |
| 24 | `XATZOK` | 3 | `0.9301` | `0.0958` | `MUDTAJFragMofMin` |
| 25 | `XAGDAL` | 2 | `0.9490` | `0.1421` | `NINTEMFragMof` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high (above 0.98). This indicates that the local coordination environment of Zinc (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 4.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent MOFs are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `4.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some Zn centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of Zn-coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Zinc environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
