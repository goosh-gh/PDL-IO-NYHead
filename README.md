# PDL::IO::NYHead

Read the **New York Head** forward model (`sa_nyhead.mat`, MATLAB v7.3 / HDF5)
into [PDL](https://pdl.perl.org). A focused reader that exposes electrodes,
leadfield, cortical/scalp surfaces, the Harvard-Oxford atlas, coordinate
transforms, and nearest-vertex lookup with a small, PDL-native API.

This distribution was split out of the `PDL-EEG` project so it can be reused
independently. It has no dependency on `PDL-EEG`.

Version: 0.04 (`$PDL::IO::NYHead::VERSION`).

## Why not PDL::IO::Matlab?

`PDL::IO::Matlab` (matio-backed) can open MAT v7.3 files, but its Perl API only
surfaces real numeric arrays: it does not traverse the nested `sa` struct and
cannot decode MATLAB cell-of-strings (the electrode and atlas names, stored as
HDF5 object references). This distribution reads numeric data with
`PDL::IO::HDF5` and recovers the string fields with `PDL::IO::NYHead::H5Dump`,
a general MATLAB-v7.3 cell-of-strings helper: it uses `h5dump` when it is found,
and otherwise reads the object references directly with the low-level functions
of `PDL::IO::HDF5`.

## Requirements

- PDL, PDL::IO::HDF5
- `h5dump` (from the HDF5 tools) is optional. The label fields
  (`electrode_labels`, `ho_labels`) are read with it when it is on `PATH` (or
  given with `h5dump => PATH`), and with the low-level functions of
  `PDL::IO::HDF5` otherwise. `PDL::IO::NYHead->new($file, strings => 'h5dump')`
  or `strings => 'lowlevel'` fixes the choice; `$ny->strings_method` reports it.
  Both give the same names on the bundled fixture (`t/04_lowlevel.t`).
- `PDL::Graphics::Cairo` (giza) — **only** for the MRI viewers:
  `tools/nyhead_ortho_gs.pl` needs the interactive giza backend; the static
  `examples/nyhead_ortho.pl` / `examples/nyhead_mri_panels.pl` render PNGs via
  Cairo.

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

    # leadfield rows at own points, interpolated from the rows of the electrodes
    my ($rows, $W) = $ny->interp_rows(pdl([[0, -113, -114.5]]), kernel => 'r3');
    # $rows (np,74382), $W (231,np); k => N uses the N nearest electrodes,
    # exclude => [names] leaves electrodes out

    # export a surface to OBJ / USDA (USDA is Z-up MNI mm)
    $ny->write_mesh('cortex75K', 'obj',  'cortex.obj');
    $ny->write_mesh('head',      'usda', 'head.usda');
    $ny->write_mesh('cortex10K', 'usda', 'cortex.usda', axes => 90); # + XYZ axis triad

## Interpolating leadfield rows

`PDL::IO::NYHead::Interp` interpolates the leadfield row at any point from the
rows of the electrodes: a 3-D polyharmonic radial basis function (`phi(r) = r`
or `r^3`) with a linear term, so the weights depend only on the positions, sum
to 1 and reproduce linear functions exactly, and an average-referenced row stays
average-referenced. The weights `W (ne,np)` are applied to all sources with one
matrix product (`interp_rows($L, $W)`). `loo` gives the leave-one-out error of
each electrode's row (correlation over the sources and RMS error relative to the
RMS of the row); `gap => G` leaves out, with each electrode, the electrodes
closer to it than G mm, so that every electrode is interpolated from the
spacing of a given point. With all electrodes it uses the closed form of
Rippa (1999), extended to the left-out set (one inverse of the system instead
of one solve per electrode; equal to the explicit solves to 1e-14 in
`t/05_interp.t`). `t/05_interp.t` checks it on
current-dipole fields over a sphere with a negative control (values unrelated
to position).

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

## Tools

`tools/nyhead_ortho_gs.pl` — an interactive 2×2 ortho viewer (axial / coronal /
sagittal) over the New York Head MRI, rendered through the giza backend of
`PDL::Graphics::Cairo`. Clicking a slice moves the linked crosshair in all three
panes and prints the clicked point as a ready-to-paste `--eec X,Y,Z` (MNI mm) on
the terminal. The viewer's two native giza sliders set the grayscale window
level (bottom) and width (right); `--zoom N` sets the initial magnification.
Zoom is data-space: the full-resolution volume is re-sliced around the crosshair
so more voxels appear instead of a magnified subsample, and while zoomed the
fourth pane switches from the intensity histogram to a 3-plane overview showing
the current zoom box. Internal air (external ear canal, sinuses, mastoid) is
separated from the external background by connected components and overlaid in
cyan. Electrodes — from `--nyhead <labels>` (via `electrode_labels` /
`electrode_pos`) or explicit `--elec 'LABEL=x,y,z'` — are drawn as colored open
circles, and each pane's bottom-right corner lists the 3-D straight-line
distance from both the original `--eec` point and the current click point to
each reference electrode, updated on every click while the original is retained.
A headless `--self-test --probe` path renders to PNG without opening a window.

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
`sa_nyhead.mat` and prints a summary with PASS/FAIL structural checks
(`--strings lowlevel` reads the names without `h5dump`).

`examples/nyhead_electrodes_3d.pl` — shows the electrodes with their labels
in the GS3D window of `PDL::Graphics::Cairo` (giza-server) over the scalp mesh
(`--surface cortex10K` etc. for a cortex instead) and lists them with their
`locs_3D_orig` coordinates: `--mark REGEX` colours electrodes by name, `--below
MM` lists the electrodes behind y = 0 and lower than Iz + MM, `--only REGEX`
limits the electrodes drawn. `--points TSV` / `--point NAME=X,Y,Z` add your own
labelled points (e.g. neck electrode sites placed with
`nyhead_neck_electrodes.pl --points-out`), each joined by a line to its nearest
NYHead electrode (`--no-link` to leave out), together with the back-midline
profile (green on the mesh, orange for its straight extension). `--list-only`
prints the tables without opening a window.

`examples/nyhead_neck_electrodes.pl` — for the neck electrodes (`Nk1`..`Nk4` by
default): position, nearest head-mesh vertex and distance, height above the
lowest head-mesh vertex, the leadfield row RMS relative to the median of all
electrodes, and the correlation of the rows (across cortical sources) with
chosen electrodes and with the nearest non-neck electrode. `--point NAME=X,Y,Z`
prints the nearest NYHead electrodes to a point; with `--rename TSV --out TSV`
they are added as columns to a channel table, and `--points-out TSV` writes
the points and the profile for `nyhead_electrodes_3d.pl --points`. `--arc 'NAME=MM,...'` places points
on the back midline of the head mesh (the cut by the plane x = 0) at lengths
measured along it downwards from Iz (`--from`), as with a tape measure along the
neck. The head mesh ends in the upper neck and is closed by a cap: the profile is
cut where it first runs forward at less than 45° from the horizontal and is
extended straight below that along the direction of its last 20 mm; points on
the extension are marked. The neck electrodes are projected onto the profile
(length along it and distance from it). For pairs of neck electrodes the script
prints the correlation of the leadfield rows and the RMS of their difference,
next to the same RMS for neighbouring scalp electrodes. `--selftest` checks the
profile, the cut, the extension and the projection on a synthetic cone with a
cap.

`--along 'NAME,...'` (with `--arc`) projects NYHead electrodes onto the back
midline profile and lists them, Iz and the `--arc` points by their length along
it, with the skin distance (difference of the lengths) between every two of
them.

`examples/nyhead_interp_rows.pl` — leadfield rows at own points (`--point
NAME=X,Y,Z`, or `--points TSV` such as the `--points-out` file above)
interpolated from the 231 electrode rows (`--method r3:all` by default; `r` or
`r3`, all electrodes or the N nearest). It prints, for each point, the nearest
electrodes, whether the point lies below the lowest electrode or more than
45 mm from any electrode (then marked as extrapolated), and the largest
weights; the leave-one-out errors of every electrode for the `--methods`, by
the distance to the nearest remaining electrode and one by one for `--show`;
for each point the leave-one-out errors with a gap equal to the point's
distance to its nearest electrode (RMS error relative to the median row RMS,
median and 90th percentile, over all electrodes and over the lower back of the
head), taken as the expected error of its row; and with `--pairs 'A-B,...'`
(electrodes and points mixed) the distance, the correlation of the rows and the
RMS of their difference relative to the median row RMS, with the expected error
of the difference and its range over the `--methods`, next to the same RMS for
neighbouring scalp electrodes. `--selftest`
runs it on a synthetic sphere.

`examples/nyhead_ortho.pl` — the same 2×2 ortho view as a static PNG.
`--click PANEL,fx,fy` (fx,fy = image fraction) resolves a click to its MNI mm
and prints the `--eec X,Y,Z` line, without opening a window.

`examples/nyhead_mri_panels.pl` — static axial / coronal / sagittal panels
(1×3 PNG) with electrode markers. For each electrode it reads the raw
`/sa/mri` intensity at that point, its percentile rank within the head, a
classification (air / bright / mid / dark tissue) and the distance to the
nearest air voxel — the numeric readout for judging where a clicked point sits
relative to bone and the ear-canal air.
