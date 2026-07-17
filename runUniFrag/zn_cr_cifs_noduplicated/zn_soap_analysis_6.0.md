# Zn-Centered SOAP Fingerprint Environment Analysis

This report evaluates the preservation of local chemical environments around Zn (`Zn`) centers in the fragment library using high-dimensional **SOAP (Smooth Overlap of Atomic Positions)** fingerprints. 

SOAP compares environments within a continuous **6.0 Å cutoff**, capturing details like geometry, coordination shells, and local density.

## Executive Summary

* **Average SOAP Cosine Similarity**: `0.9908`
* **Median SOAP Cosine Similarity**: `0.9957`
* **Average SOAP Fingerprint RMSD**: `0.0564`
* **Median SOAP Fingerprint RMSD**: `0.0452`
* **Total Zn Centers Analyzed**: `9176`

| Similarity Range | Category | Count | Percentage | Description |
| :--- | :--- | :---: | :---: | :--- |
| **$\ge 0.98$** | Highly Represented | 7734 | **84.29%** | Local environment is almost perfectly preserved in the fragment library. |
| **$[0.90, 0.98)$** | Moderately Represented | 1424 | **15.52%** | Local environment is structurally similar, with minor variations (e.g. capped bonds, minor coordinates shift). |
| **$< 0.90$** | Poorly Represented / Missing | 18 | **0.20%** | Environment has significant structural/coordination divergence in the fragment library. |

## PCA and UMAP Environment Distribution Map

Below are 2D PCA and UMAP projections of the SOAP fingerprints. The close overlap between parent crystal centers (blue) and fragment library centers (green) visually demonstrates excellent chemical coverage:

![Zn SOAP PCA and UMAP Map](zn_soap_distribution_6.0.png)

## Poorly Represented Zn Environments (Bottom 25 Worst Matches)

These Zn centers in parent crystals have the lowest similarity scores to any fragment in the library, highlighting potential distortions caused by capping or trimming:

| Rank | Parent REFCODE | Zn Index | Max Cosine Similarity | Fingerprint RMSD | Best Matching Fragment |
| :---: | :---: | :---: | :---: | :---: | :--- |
| 1 | `REDYEG` | 1 | `0.8010` | `0.1998` | `REDYEGFragMof` |
| 2 | `REDYEG` | 0 | `0.8010` | `0.1998` | `REDYEGFragMof` |
| 3 | `REDYEG` | 2 | `0.8105` | `0.2146` | `REDYEGFragMof` |
| 4 | `REDYEG` | 4 | `0.8105` | `0.2146` | `REDYEGFragMof` |
| 5 | `REDYEG` | 5 | `0.8105` | `0.2146` | `REDYEGFragMof` |
| 6 | `REDYEG` | 3 | `0.8105` | `0.2146` | `REDYEGFragMof` |
| 7 | `FIDMEL` | 2 | `0.8752` | `0.3001` | `FIDMELFragMof` |
| 8 | `FIDMEL` | 0 | `0.8768` | `0.2951` | `FIDMELFragMof` |
| 9 | `FIDMEL` | 3 | `0.8769` | `0.2973` | `FIDMELFragMof` |
| 10 | `FIDMEL` | 5 | `0.8771` | `0.2959` | `FIDMELFragMof` |
| 11 | `FIDMEL` | 4 | `0.8778` | `0.2915` | `FIDMELFragMof` |
| 12 | `FIDMEL` | 1 | `0.8782` | `0.2942` | `FIDMELFragMof` |
| 13 | `RUMRUO` | 0 | `0.8857` | `0.3204` | `OQOXAVFragMof` |
| 14 | `RUMRUO` | 1 | `0.8857` | `0.3204` | `OQOXAVFragMof` |
| 15 | `OQOXAV` | 3 | `0.8863` | `0.3245` | `OQOXAVFragMof` |
| 16 | `OQOXAV` | 2 | `0.8863` | `0.3245` | `OQOXAVFragMof` |
| 17 | `OQOXAV` | 5 | `0.8863` | `0.3245` | `OQOXAVFragMof` |
| 18 | `OQOXAV` | 4 | `0.8863` | `0.3245` | `OQOXAVFragMof` |
| 19 | `XATZOK` | 0 | `0.9147` | `0.6029` | `CEFKACFragMof` |
| 20 | `XATZOK` | 3 | `0.9147` | `0.6029` | `CEFKACFragMof` |
| 21 | `IYOQAQ` | 0 | `0.9253` | `0.1905` | `DAFSOVFragMof` |
| 22 | `IYOQAQ` | 1 | `0.9253` | `0.1905` | `DAFSOVFragMof` |
| 23 | `FIDMEL` | 7 | `0.9297` | `0.2571` | `FIDMELFragMof` |
| 24 | `LACJAC` | 0 | `0.9327` | `0.1754` | `WEVDUAFragMofMin` |
| 25 | `WAWGOQ` | 5 | `0.9345` | `0.5459` | `CEFKACFragMof` |

## Discussion & Chemical Analysis

1. **High Overall Similarity**:
   The median similarity of SOAP descriptors is extremely high. This indicates that the local coordination environment of Zn (including coordination shell composition, distance distribution, and local symmetry) is well preserved by the UniFrag extraction algorithm within the 6.0 Å sphere.
   
2. **Periodic vs Non-Periodic Context**:
   Because SOAP descriptors for parent structures are calculated with `periodic=True` (capturing atoms extending outside the unit cell boundaries) while fragments are computed with `periodic=False` (treating them as isolated molecules), some divergence is expected. The fact that the overlap is so tight demonstrates that the `6.0 Å` extraction shell captures almost all relevant local chemical details.
   
3. **Capping Effects**:
   Capped terminals (like O-H, C-H) introduce small hydrogen atoms at boundaries that were originally occupied by other framework atoms. This contributes to moderate similarity values ($0.90 - 0.98$) for some metal centers located close to linker cut sites.

4. **PCA vs UMAP Projection**:
   * **PCA** shows the global directions of largest linear variance, capturing the primary geometric axes of metal coordination variations across the dataset.
   * **UMAP** preserves non-linear local neighborhood structures. The tight grouping and consistent overlap in UMAP space further verify that the fragment library does not form isolated topological clusters detached from the parent distributions, but rather covers the continuous space of parent environments.

## Conclusion

The SOAP fingerprint analysis confirms that **the fragment library provides exceptional, continuous structural coverage of the local Zn environments** in the parent crystal structures, making the generated fragments highly representative models for downstream Quantum Chemical (QM) calculations.
