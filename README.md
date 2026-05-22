# Codes for the project *Crystallography, Group Cohomology, and Lieb–Schultz–Mattis Constraints* 
# 

[![arXiv](https://img.shields.io/badge/arXiv-2410.03607-b31b1b.svg)](https://arxiv.org/abs/2410.03607)

## Content

This repository contains the following three parts

### GAP codes

The GAP sources are organised as a proper GAP package, living in
[`gap_codes/SpaceGroupCohomology/`](gap_codes/SpaceGroupCohomology/).  The
package depends on [HAP](https://gap-packages.github.io/hap/) (≥ 1.30) and
loads it automatically.

**Install** by symlinking the package directory into GAP's package search
path (typically `~/.gap/pkg/`):

```bash
mkdir -p ~/.gap/pkg
ln -s "$(pwd)/gap_codes/SpaceGroupCohomology" ~/.gap/pkg/SpaceGroupCohomology
```

**Use** the single command:

```gap
gap> LoadPackage("SpaceGroupCohomology");
true
gap> SpaceGroupCohomologyRingGapInterface(IT);   # 1 <= IT <= 230
```

`SpaceGroupCohomologyRingGapInterface(IT)` prints the mod-2 cohomology ring
of space group No. `IT` together with its LSM anomaly classes.

The package layout, dependencies and self-test are documented in
[`gap_codes/SpaceGroupCohomology/README.md`](gap_codes/SpaceGroupCohomology/README.md).

> **Note:** The older flat files `gap_codes/SpaceGroupCohomologyData.gi` and
> `gap_codes/SpaceGroupCohomologyFunctions.gi` are now thin deprecation
> shims that point to the new package.  They will be removed in a future
> release; please update any scripts that `Read` them directly.

### Mathematica codes

Space_Group_Cohomology_Data.nb contains the standard inhomogeneous functions for the 1-,2-, and 3-cocycles (except for No. 225, 227, and 229)
Execute the "Preparations" section first.

### Sage codes

Codes to obtain the explicit cochain representatives for degree-3 generators of several complicated space groups


## Reference

If this repository is useful for your research, please consider citing the [Scipost](https://www.scipost.org/SciPostPhys.18.5.161) article:

```bibtex
@Article{10.21468/SciPostPhys.18.5.161,
	title={{Crystallography, group cohomology, and Lieb–Schultz–Mattis constraints}},
	author={Chunxiao Liu and Weicheng Ye},
	journal={SciPost Phys.},
	volume={18},
	pages={161},
	year={2025},
	publisher={SciPost},
	doi={10.21468/SciPostPhys.18.5.161},
	url={https://scipost.org/10.21468/SciPostPhys.18.5.161},
}
```
