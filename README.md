# PDL::IO::NYHead

Read the **New York Head** forward model (`sa_nyhead.mat`, MATLAB v7.3 / HDF5)
into [PDL](https://pdl.perl.org). A focused reader that exposes electrodes,
leadfield, cortical/scalp surfaces, the Harvard-Oxford atlas, coordinate
transforms, and nearest-vertex lookup with a small, PDL-native API.

This distribution was split out of the `PDL-EEG` project so it can be reused
independently. It has no dependency on `PDL-EEG`.

## Why not PDL::IO::Matlab?

`PDL::IO::Matlab` (matio-backed) can open MAT v7.3 files, but its Perl API only
surfaces real numeric arrays: it does not traverse the nested `sa` struct and
cannot decode MATLAB cell-of-strings (the electrode and atlas names, stored as
HDF5 object references). This distribution reads numeric data with
`PDL::IO::HDF5` and recovers the string fields via `h5dump`
(`PDL::IO::NYHead::H5Dump`), which is a genuinely general MATLAB-v7.3
cell-of-strings helper.

## Requirements

- PDL, PDL::IO::HDF5
- `h5dump` (from the HDF5 tools) on `PATH` — **only** needed for the label
  fields (`electrode_labels`, `ho_labels`). Everything numeric works without it.
  MacPorts: `port install hdf5`.
- `PDL::Graphics::Cairo` (giza) — **only** for the MRI viewers in `examples/`
  (`nyhead_ortho_gs.pl` needs the interactive giza backend; the static
  `nyhead_ortho.pl` / `nyhead_mri_panels.pl` render PNGs via Cairo).

## Install

    perl Makefile.PL
    make
    make test
    make install

## Synopsis

    use PDL::IO::NYHead;
    my $ny = PDL::IO::NYHead->new('sa_nyhead.mat');

    my $labels = $ny->electrode_labels;       # 231 names
    my $pos    = $ny->electrode_pos;          # (231,3) MNI mm
    my $L      = $ny->leadfield;              # (231,74382) normal-oriented CAR

    my $ctx    = $ny->surface('cortex75K');   # { vc, tri, tri0 }
    my $lo     = $ny->cortex('cortex2K');     # reconstructed from cortex75K

    my $areas  = $ny->ho_labels;              # 97 Harvard-Oxford area names
    print $ny->area_of_vertex(1234), "\n";    # e.g. 'Left Frontal Pole'

    my $iv = $ny->nearest_vertex($pos, 'head'); # each electrode -> scalp vertex

    # export a surface to OBJ / USDA (USDA is Z-up MNI mm)
    $ny->write_mesh('cortex75K', 'obj',  'cortex.obj');
    $ny->write_mesh('head',      'usda', 'head.usda');
    $ny->write_mesh('cortex10K', 'usda', 'cortex.usda', axes => 90); # + XYZ axis triad

## Axis convention

h5dump reports shapes `{A,B,...}`; PDL reads them reversed as `(...,B,A)`. The
New York Head stores geometry as `{3,N}`, so `vc` / `normals` / `locs` / `tri`
come back as `(N,3)` / `(Nf,3)` directly (no transpose). `tri` is 1-based
(kept for OBJ export; `tri0` is the 0-based copy). `V_fem_normal` is
`(electrode, source)`.

USDA export carries `upAxis = "Z"` and `metersPerUnit = 0.001` (the New York
Head is Z-up, millimetres), so a surface stands upright in usdview / Keynote.
`write_mesh($surf, 'usda', $out, axes => <mm>)` adds an origin axis triad at the
MNI origin (red = +X, green = +Y, blue = +Z); the default (`axes => 0`) writes
the mesh alone. `examples/nyhead_export.pl --format usda --axes <mm>` is the CLI.

The MRI viewers map MNI mm ↔ voxel with the `sa` affines
(`mri2mni` / `mni2mri`), auto-orienting each to math convention by the
`[0,0,0,1]` bottom-row invariant (robust to the HDF5 axis reversal). The
background polarity (New York Head pads air at the bright extreme) is detected
automatically and only the border-connected external background is blanked, so
internal air cavities survive for the cyan overlay. The default grayscale window
is chosen by maximising the displayed-histogram entropy (auto-contrast).

## See also

`examples/nyhead_export.pl` — export any surface to OBJ or USDA
(`--format obj|usda`, `--surf head|cortexNK`, `--axes <mm>` for a USDA axis
triad).

`examples/nyhead_mri_to_nii.pl` — export the New York Head MRI volume
(`/sa/mri`) to NIfTI-1. Writes an MNI `sform` (`sform_code = 4`) and streams the
voxels through PDL's data pointer (`get_dataref`, no intermediate Perl list);
`--dtype float32` halves the file. The written `.nii` reads back identically in
both nibabel and `PDL::IO::BIDS` (affine and voxel↔world round-trip verified).

`examples/nyhead_ear_usda.pl` — write the scalp, the three BEM shells and a set
of electrode / marker spheres to one `.usda` for usdview / Keynote / Blender.
The skin (`/sa/head`) is a see-through `BasisCurves` wireframe cage, the BEM
shells (`/#refs#/b,c,d`) are solid opaque meshes, the electrodes translucent
spheres, and named target / candidate points opaque spheres, with an optional
XYZ axis triad. `upAxis = "Z"`, `metersPerUnit = 0.001`.

`examples/check_nyhead.pl` — runs the full reader against a real
`sa_nyhead.mat` and prints a summary with PASS/FAIL structural checks.

`examples/nyhead_ortho_gs.pl` — an interactive 2×2 ortho viewer (axial /
coronal / sagittal + an intensity histogram) over the New York Head MRI,
rendered through the giza backend of `PDL::Graphics::Cairo`. Clicking a slice
moves the linked crosshair in all three panes and prints the clicked point as a
ready-to-paste `--eec X,Y,Z` (MNI mm) on the terminal. Two native giza sliders
set the grayscale window level and width; internal air (external ear canal,
sinuses, mastoid) is separated from the external background by connected
components and overlaid in cyan. Electrodes come from `--nyhead <labels>`
(via `electrode_labels` / `electrode_pos`) or explicit `--elec 'LABEL=x,y,z'`.

`examples/nyhead_ortho.pl` — the same 2×2 ortho view as a static PNG.
`--click PANEL,fx,fy` (fx,fy = image fraction) resolves a click to its MNI mm
and prints the `--eec X,Y,Z` line, without opening a window.

`examples/nyhead_mri_panels.pl` — static axial / coronal / sagittal panels
(1×3 PNG) with electrode markers. For each electrode it reads the raw
`/sa/mri` intensity at that point, its percentile rank within the head, a
classification (air / bright / mid / dark tissue) and the distance to the
nearest air voxel — the numeric readout for judging where a clicked point sits
relative to bone and the ear-canal air.
