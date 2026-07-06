# SpaceGroupCohomology

A GAP package for computing mod-2 cohomology rings and
Lieb-Schultz-Mattis anomaly classes for the 230 crystallographic space
groups.

Reference:

Chunxiao Liu and Weicheng Ye, *Crystallography, group cohomology, and
Lieb-Schultz-Mattis constraints*, SciPost Phys. **18**, 161 (2025)
([arXiv:2410.03607](https://arxiv.org/abs/2410.03607),
[SciPost](https://scipost.org/10.21468/SciPostPhys.18.5.161)).

## What This Package Does

For a space group specified by its International Tables number `IT`, the
main interface

```gap
SpaceGroupCohomologyRingGapInterface(IT);
```

prints:

- a presentation of the mod-2 cohomology ring,
- the named generators used in the presentation,
- the relations found up to the required degree, and
- the LSM anomaly classes associated with the irreducible Wyckoff
  positions recorded in the package data.

The computation is based on HAP's `ResolutionSpaceGroup`, together with
package code for cup products, mod-2 linear algebra, generator/relations
selection, point-group bookkeeping, and LSM evaluation.

## Installation

Copy or symlink this directory into a GAP package search path, for example:

```bash
mkdir -p ~/.gap/pkg
ln -s /path/to/SpaceGroupCohomology ~/.gap/pkg/SpaceGroupCohomology
```

The package depends on:

- GAP,
- HAP, version at least 1.30,
- polymake on the `PATH`.

HAP is loaded automatically by GAP. Polymake is used through HAP's
crystallographic-complex machinery to compute fundamental domains and face
lattices.

On macOS, for example, polymake can be installed with:

```bash
brew install polymake
```

## Basic Usage

```gap
gap> LoadPackage("SpaceGroupCohomology");
true

gap> SpaceGroupCohomologyRingGapInterface(98);
```

The input must be an integer `IT` with `1 <= IT <= 230`.

The interface chooses the resolution length automatically:

- length 9 for groups currently known to require degree-7/8 relation checks,
- length 7 for groups `IT <= 220`,
- length 6 for groups `IT > 220`.

## Current Computational Pipeline

The public ring interface now uses

```gap
SGC_ResolutionSpaceGroupSearchHomotopy(G, n)
```

for all space groups. This builds HAP's small-memory
`ResolutionSpaceGroup(G, n)` resolution and attaches a mod-2 contracting
homotopy by finite GF(2) filling. This is the path used by cup products in
the mod-2 cohomology-ring calculation.

This avoids the memory bottleneck of `ResolutionAlmostCrystalGroup`, while
still providing the homotopy data needed by HAP's cup-product machinery.

The package also contains an experimental integral homotopy constructor,

```gap
SGC_ResolutionSpaceGroup(G, n)
```

which builds a `ResolutionSpaceGroup`-compatible free resolution with an
integral contracting homotopy coming from cellular contraction data. This is
useful for testing and future development, but the mod-2 ring interface does
not currently rely on it: some Voronoi-shelling cases have genuine local
shelling obstructions, whereas the mod-2 filling homotopy is sufficient for
the present cohomology-ring calculations.

## Recent Fixes

This version includes several fixes and robustness improvements:

- The public ring interface no longer uses a special case for space group
  No. 98. It uses the same `ResolutionSpaceGroup` plus mod-2 filling
  homotopy path as the other groups.
- The point-group normal-form code now validates the enumerated point-group
  words by comparing against the actual finite point-group size and checking
  uniqueness. This fixes the earlier failures for groups `75..88`, where
  the second point-group generator has order 4 even though the chosen normal
  form only uses exponent `0` or `1`.
- The HAP model used in the integral homotopy code for the relevant
  `SL(2,O-2)` complex was updated to the variant with the stabilizer data
  expected by current HAP.
- The cellular old-cell mask in the experimental Voronoi contraction code is
  closed downward under faces, preventing non-subcomplex masks.
- Mod-2 linear algebra can use sparse wrappers and, optionally, an external
  binary accelerator.

## Optional Linear-Algebra Accelerator

Large mod-2 linear-algebra computations can optionally be offloaded to a
small bundled C++ program.

Build it with:

```bash
make -C cpp
make -C cpp test
```

If `bin/sgclinalg` is present, matrices whose size exceeds
`SGC_LINALG_THRESHOLD` are handled through a sparse file interface. Smaller
matrices remain on GAP-native routines.

Useful controls:

```gap
SGC_LINALG_THRESHOLD := 10^6;  # default-style threshold
SGC_LINALG_THRESHOLD := 0;     # force external offload when available
```

The environment variable `SGC_LINALG_PATH` can point to an external
`sgclinalg` binary outside the package directory.

## Package Layout

```text
SpaceGroupCohomology/
├── PackageInfo.g       package metadata
├── init.g              declarations of package globals
├── read.g              loads data, linalg, homotopy, and ring code
├── gap/
│   ├── data.gi         PGGens230, IWP, GENNAMES, funcs230
│   ├── linalg.gi       mod-2 linear-algebra wrappers and sparse helpers
│   ├── homotopy.gi     ResolutionSpaceGroup homotopy constructors
│   └── functions.gi    cohomology rings, point-group bookkeeping, LSM code
├── cpp/
│   └── ...             optional C++ linear-algebra accelerator source
├── bin/
│   └── sgclinalg       optional compiled accelerator, if built
└── tst/
    ├── smoke.tst       end-to-end ring computations for sample groups
    ├── linalg.tst      mod-2 linear-algebra tests
    └── homotopy.tst    integral homotopy consistency tests
```

## Testing

Run the package smoke test with:

```gap
gap> TestPackage("SpaceGroupCohomology");
```

or run individual tests:

```gap
gap> Test("tst/smoke.tst");
gap> Test("tst/linalg.tst");
gap> Test("tst/homotopy.tst");
```

The tests cover:

- data-table loading for all 230 groups,
- end-to-end ring computations for representative groups,
- mod-2 linear-algebra correctness, with and without the external binary,
- consistency of the experimental integral homotopy on a sample group.

Note: GAP transcript tests are sensitive to exact whitespace. Long LSM
output lines may produce harmless formatting diffs on some terminals even
when the mathematical output is unchanged.

## Developer Notes

The main public function is:

```gap
SpaceGroupCohomologyRingGapInterface(IT)
```

Important internal entry points include:

```gap
Mod2RingGenerators(R, deg, spacedim)
Mod2RingGensAndRels(IT, spacedim, R, gens)
SGC_ResolutionSpaceGroupSearchHomotopy(G, n)
SGC_ResolutionSpaceGroup(G, n)
SGC_HomotopyCheckMod2(R, k, letter)
SGC_TotalHomotopyCheck(R, k, letter)
```

For ring calculations, prefer `SGC_ResolutionSpaceGroupSearchHomotopy`.
Use `SGC_ResolutionSpaceGroup` when explicitly testing or developing the
integral cellular contraction.

## Citation

If this package is useful for your research, please cite:

```bibtex
@Article{10.21468/SciPostPhys.18.5.161,
  title={{Crystallography, group cohomology, and Lieb-Schultz-Mattis constraints}},
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
