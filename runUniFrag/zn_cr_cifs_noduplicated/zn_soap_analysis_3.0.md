# Zn-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Zinc (`Zn`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **3.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9992`
* **Median SOAP Cosine Similarity**: `0.9999`
* **Average SOAP Fingerprint RMSD**: `0.0037`
* **Median SOAP Fingerprint RMSD**: `0.0018`
* **Total Zn Centers Analyzed**: `9176`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 9112 | **99.30%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 60 | **0.65%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 4 | **0.04%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Zn SOAP PCA and UMAP Map](file:///Users/omert/.gemini/antigravity/brain/153d2da2-7e4a-4474-b36a-6be8db573d0d/zn_soap_distribution_3.0.png)

## Poorly Represented Zn Environments (Bottom 25 Worst Matches)

These Zn centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Zn Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `OQOXAV` | 4 | `0.8999` | `0.0673` | `KUSTEAFragMof` |
| 2 | `OQOXAV` | 2 | `0.8999` | `0.0673` | `KUSTEAFragMof` |
| 3 | `OQOXAV` | 3 | `0.8999` | `0.0673` | `KUSTEAFragMof` |
| 4 | `OQOXAV` | 5 | `0.8999` | `0.0673` | `KUSTEAFragMof` |
| 5 | `RUMRUO` | 0 | `0.9001` | `0.0672` | `KUSTEAFragMof` |
| 6 | `RUMRUO` | 1 | `0.9001` | `0.0672` | `KUSTEAFragMof` |
| 7 | `REDYEG` | 0 | `0.9262` | `0.0615` | `REDYEGFragMof` |
| 8 | `REDYEG` | 1 | `0.9262` | `0.0615` | `REDYEGFragMof` |
| 9 | `MOMRIR` | 15 | `0.9492` | `0.0659` | `LEDLAJFragMofMin` |
| 10 | `MOMRIR` | 14 | `0.9492` | `0.0659` | `LEDLAJFragMofMin` |
| 11 | `MOMRIR` | 7 | `0.9507` | `0.0638` | `LEDLAJFragMofMin` |
| 12 | `MOMRIR` | 9 | `0.9507` | `0.0638` | `LEDLAJFragMofMin` |
| 13 | `MOMRIR` | 6 | `0.9507` | `0.0638` | `LEDLAJFragMofMin` |
| 14 | `MOMRIR` | 8 | `0.9507` | `0.0638` | `LEDLAJFragMofMin` |
| 15 | `JISMEG` | 5 | `0.9576` | `0.0664` | `BIBQOVFragMof` |
| 16 | `JISMEG` | 3 | `0.9576` | `0.0664` | `BIBQOVFragMof` |
| 17 | `JISMEG` | 7 | `0.9576` | `0.0664` | `BIBQOVFragMof` |
| 18 | `JISMEG` | 1 | `0.9576` | `0.0664` | `BIBQOVFragMof` |
| 19 | `REDYEG` | 2 | `0.9628` | `0.0417` | `REDYEGFragMof` |
| 20 | `REDYEG` | 4 | `0.9628` | `0.0417` | `REDYEGFragMof` |
| 21 | `REDYEG` | 3 | `0.9628` | `0.0417` | `REDYEGFragMof` |
| 22 | `REDYEG` | 5 | `0.9628` | `0.0417` | `REDYEGFragMof` |
| 23 | `VEMDAT` | 0 | `0.9743` | `0.0272` | `RUGZAWFragMofMin` |
| 24 | `VEMDAT` | 1 | `0.9743` | `0.0272` | `RUGZAWFragMofMin` |
| 25 | `VEMDAT` | 2 | `0.9743` | `0.0272` | `RUGZAWFragMofMin` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high (above 0.98). This indicates that the local coordination environment of Zinc (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 3.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent MOFs are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `3.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some Zn centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of Zn-coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Zinc environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
