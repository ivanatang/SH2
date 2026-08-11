"""
sh2_peptide_movie.py

Render a movie of the SH2 superbinder + pTyr peptide production trajectory
(50 ns, PBC-corrected and fit by postprocess_SH2.sh), for visualizing whether
the peptide stays engaged in the engineered binding pocket.

Adapted from biosensors' pymol_renders/gate_latch_movie.py -- same design
(fixed camera set once, drift removed via intra_fit before any camera call,
non-ray-traced per-frame PNGs stitched with ffmpeg, simulation-time label
placed via a fixed-in-model-space pseudoatom) -- but two things differ
because this is a different kind of system:

  1. No cartoon for the peptide. The GROMACS topology (built in
     SH2_superbinder_pTyr_peptide_parameterization.ipynb) writes the entire
     11-residue peptide as a SINGLE PDB residue (chain B, resn LIG, resi 110
     in the exported frame -- confirmed via grep on PL_only_frame0.pdb) to
     satisfy an earlier GROMACS export requirement (see notebook). PyMOL's
     cartoon algorithm traces backbone per-residue, so a lumped single-residue
     "peptide" cannot be drawn as a ribbon -- it would just show as a blob or
     nothing. The peptide is rendered as sticks instead (confirmed via a local
     bond-perception test: all 208 peptide atoms resolve bonds correctly via
     PyMOL's default distance-based bonding, so sticks render properly).

  2. No gate/latch loops to exclude from the fit. SH2 domains don't have the
     PYR1-style gate-loop/latch/recoil architecture that gate_latch_movie.py's
     FLEXIBLE_EXCLUDE_RESI carves out. The fit here uses ALL protein chain-A
     CA atoms as the reference -- the whole point is to hold the SH2 domain
     itself rock-steady so that any visible motion is the peptide's, since
     the peptide (not a protein loop) is the thing whose dynamics we care
     about (does it stay in the pocket or not).

Highlighted elements (see postprocess_SH2.sh / phosphate_pocket_mindist.xvg
for the numeric version of what this movie shows visually):
  - Phosphate group (P, O1P, O2P, O3P) of the pTyr residue -- the primary
    specificity contact for SH2 domains.
  - The three engineered superbinder mutations (chain A resi 40/45/63 in
    this system's numbering = V183/A188/L206 in the paper's UniProt
    numbering, per af3/README.md's AF3_resi = UniProt_resi - 143 offset,
    confirmed by residue identity: VAL/ALA/LEU at 40/45/63 in
    sh2_superbinder_pTyr_dodecahedron_protein.itp).

Run non-interactively with the local PyMOL build:

    /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/SH2/md_prep/pymol_renders/sh2_peptide_movie.py

Smoke test (renders only a handful of frames, fast, writes to a
"_smoketest" suffixed mp4 so it never clobbers the real deliverable):

    SH2_MOVIE_MAX_FRAMES=25 /opt/homebrew/bin/pymol -cq /Users/ivanatang/Developer/SH2/md_prep/pymol_renders/sh2_peptide_movie.py

Other environment variable overrides (all optional):
    SH2_MOVIE_STRIDE        raw-trajectory-frame stride between rendered
                             output frames (default 3; PL_only.xtc has 1001
                             raw frames at 50 ps spacing spanning t=0-50 ns,
                             so stride 3 gives 334 output frames, i.e. ~33.4 s
                             of video at 10 fps -- same ~30-40 s pacing
                             convention established in gate_latch_movie.py)
    SH2_MOVIE_MAX_FRAMES    if set, caps the number of OUTPUT frames loaded
                             / rendered (smoke-test mode; also switches the
                             output filename to a "_smoketest" suffix)
    SH2_MOVIE_FPS            output movie frame rate (default 10)
    SH2_MOVIE_KEEP_FRAMES    "1" to keep the temp PNG frame directory after
                             a successful ffmpeg run (default "0", deletes)
    FFMPEG_BIN               path to the ffmpeg binary (default
                              /opt/homebrew/bin/ffmpeg)
"""

import os
import shutil
import subprocess
import tempfile

from pymol import cmd, util

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
REPO_ROOT = "/Users/ivanatang/Developer/SH2"
OUT_DIR = os.path.join(REPO_ROOT, "md_prep", "pymol_renders", "output")

DATA_DIR = (
    "/Users/ivanatang/Library/CloudStorage/OneDrive-UCB-O365/Shirts Lab/"
    "SH2_enhanced_sampling/analysis"
)
PDB = os.path.join(DATA_DIR, "PL_only_frame0.pdb")
XTC = os.path.join(DATA_DIR, "PL_only.xtc")

# 1001 raw frames at 50 ps spacing, spanning simulation time 0-50 ns (raw
# frame index 0 = t = 0.0 ns), confirmed via postprocess_SH2.sh's gmx energy /
# gmx distance output ("Analyzed 1001 frames, last time 50000.000").
RAW_N_FRAMES_TOTAL = 1001
TIME_SPACING_NS = 0.05  # 50 ps
T0_NS = 0.0

PROTEIN_CHAIN = "A"
PEPTIDE_CHAIN = "B"
PEPTIDE_RESN = "LIG"

# Superbinder mutation sites in THIS system's residue numbering (see module
# docstring). Confirmed against sh2_superbinder_pTyr_dodecahedron_protein.itp:
# resi 40=VAL, 45=ALA, 63=LEU = T183V/C188A/K206L.
POCKET_MUT_RESI = "40+45+63"
PHOSPHATE_ATOMS = "P+O1P+O2P+O3P"

# Consistent landmark colors (IBM colorblind-safe palette, same convention as
# biosensors/pymol_renders scripts)
PHOSPHATE_COLOR_HEX = "#FE6100"     # orange
POCKET_MUT_COLOR_HEX = "#785EF0"    # purple
PEPTIDE_COLOR_HEX = "#648FFF"       # blue
PROTEIN_GRAY = "gray80"

VIEWPORT_W = 1600
VIEWPORT_H = 900

FFMPEG_BIN = os.environ.get("FFMPEG_BIN", "/opt/homebrew/bin/ffmpeg")

# --------------------------------------------------------------------------
# Env-var-driven run parameters
# --------------------------------------------------------------------------
STRIDE = int(os.environ.get("SH2_MOVIE_STRIDE", "3"))
FPS = int(os.environ.get("SH2_MOVIE_FPS", "10"))
KEEP_FRAMES = os.environ.get("SH2_MOVIE_KEEP_FRAMES", "0") == "1"
_max_frames_env = os.environ.get("SH2_MOVIE_MAX_FRAMES", "").strip()

FULL_N_OUTPUT_FRAMES = (RAW_N_FRAMES_TOTAL - 1) // STRIDE + 1

if _max_frames_env:
    IS_SMOKE_TEST = True
    N_OUTPUT_FRAMES = min(int(_max_frames_env), FULL_N_OUTPUT_FRAMES)
else:
    IS_SMOKE_TEST = False
    N_OUTPUT_FRAMES = FULL_N_OUTPUT_FRAMES

RAW_STOP = min(RAW_N_FRAMES_TOTAL, (N_OUTPUT_FRAMES - 1) * STRIDE + 1)

OUT_NAME = "sh2_superbinder_pTyr_peptide_0-50ns"
if IS_SMOKE_TEST:
    OUT_NAME += "_smoketest"
OUT_MP4 = os.path.join(OUT_DIR, OUT_NAME + ".mp4")


# --------------------------------------------------------------------------
# Small color helper
# --------------------------------------------------------------------------
def hex_to_rgb01(hex_code):
    hex_code = hex_code.lstrip("#")
    return tuple(int(hex_code[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _sub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def _add(a, b, scale=1.0):
    return tuple(a[i] + b[i] * scale for i in range(3))


def _norm(v):
    length = (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) ** 0.5
    if length < 1e-6:
        return (0.0, 0.0, 1.0)
    return tuple(c / length for c in v)


def _centroid(model):
    n = len(model.atom)
    sx = sum(a.coord[0] for a in model.atom) / n
    sy = sum(a.coord[1] for a in model.atom) / n
    sz = sum(a.coord[2] for a in model.atom) / n
    return (sx, sy, sz)


# --------------------------------------------------------------------------
# Scene construction
# --------------------------------------------------------------------------
def setup_global_render_settings():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("ray_trace_fog", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.2)
    cmd.set("ambient", 0.4)
    cmd.set("defer_builds_mode", 3)


def load_topology_and_trajectory():
    cmd.load(PDB, "mol")
    n_atoms_pdb = cmd.count_atoms("mol")
    print(f"[sh2_movie] loaded topology: {n_atoms_pdb} atoms from {PDB}")

    cmd.load_traj(
        XTC,
        "mol",
        state=1,
        start=1,
        stop=RAW_STOP,
        interval=STRIDE,
        format="xtc",
    )
    n_states = cmd.count_states("mol")
    print(
        f"[sh2_movie] loaded {n_states} states from trajectory "
        f"(stride={STRIDE}, raw stop frame={RAW_STOP} of {RAW_N_FRAMES_TOTAL})"
    )

    cmd.remove("mol and hydro")
    print(f"[sh2_movie] after removing hydrogens: {cmd.count_atoms('mol')} atoms")

    return n_states


def define_selections():
    protein_sel = f"mol and chain {PROTEIN_CHAIN} and polymer"
    peptide_sel = f"mol and chain {PEPTIDE_CHAIN} and resn {PEPTIDE_RESN}"
    phosphate_sel = f"{peptide_sel} and name {PHOSPHATE_ATOMS}"
    pocket_mut_sel = f"{protein_sel} and resi {POCKET_MUT_RESI}"
    core_fit_sel = f"{protein_sel} and name CA"

    cmd.select("protein_sel", protein_sel)
    cmd.select("peptide_sel", peptide_sel)
    cmd.select("phosphate_sel", phosphate_sel)
    cmd.select("pocket_mut_sel", pocket_mut_sel)
    cmd.select("landmark_sel", "peptide_sel or pocket_mut_sel")
    cmd.select("core_fit_sel", core_fit_sel)

    n_core = cmd.count_atoms("core_fit_sel")
    n_pep = cmd.count_atoms("peptide_sel")
    n_phos = cmd.count_atoms("phosphate_sel")
    n_mut = cmd.count_atoms("pocket_mut_sel")
    print(
        f"[sh2_movie] selections: core_fit={n_core} atoms, peptide={n_pep} atoms, "
        f"phosphate={n_phos} atoms, pocket_mut={n_mut} atoms"
    )
    if n_core == 0 or n_pep == 0 or n_phos == 0 or n_mut == 0:
        raise RuntimeError(
            "[sh2_movie] one or more required selections is empty; check "
            "chain/resi/resn conventions before rendering."
        )
    if n_phos != 4:
        raise RuntimeError(
            f"[sh2_movie] expected 4 phosphate atoms (P,O1P,O2P,O3P), found {n_phos}."
        )


def remove_rigid_body_drift():
    """Fit every loaded state onto state 1 using ONLY protein CA atoms (all
    of them -- SH2 has no gate/latch-style loop to exclude the way
    gate_latch_movie.py does), so apparent peptide motion in the movie is
    real, not whole-complex tumbling. The trajectory already had a coarser
    gmx trjconv -fit rot+trans applied (on Protein_Peptide, i.e. including
    the peptide) in postprocess_SH2.sh; this additional protein-only fit
    removes any residual wobble from that broader fit so the SH2 domain
    itself is rock-steady and the peptide's true motion relative to it is
    what's visible.
    """
    rms_list = cmd.intra_fit("core_fit_sel", state=1)
    print(
        f"[sh2_movie] intra_fit on core_fit_sel done; "
        f"max per-state RMSD to state 1 = {max(rms_list):.3f} A"
    )


def style_scene():
    cmd.hide("everything", "mol")

    cmd.show("cartoon", "protein_sel")
    cmd.color(PROTEIN_GRAY, "protein_sel")
    cmd.set("cartoon_transparency", 0.15, "protein_sel")

    cmd.set_color("phosphate_color", list(hex_to_rgb01(PHOSPHATE_COLOR_HEX)))
    cmd.set_color("pocket_mut_color", list(hex_to_rgb01(POCKET_MUT_COLOR_HEX)))
    cmd.set_color("peptide_color", list(hex_to_rgb01(PEPTIDE_COLOR_HEX)))

    # Peptide as sticks (not cartoon -- see module docstring point 1: the
    # whole 11-residue peptide is a single PDB residue, so there's no
    # per-residue backbone for PyMOL's cartoon tracer to follow).
    cmd.show("sticks", "peptide_sel")
    cmd.set("stick_radius", 0.16, "peptide_sel")
    cmd.color("peptide_color", "peptide_sel and elem C")
    util.cnc("peptide_sel")

    # Phosphate group: bigger sticks, distinct color, fully opaque.
    cmd.show("sticks", "phosphate_sel")
    cmd.set("stick_radius", 0.28, "phosphate_sel")
    cmd.color("phosphate_color", "phosphate_sel and elem P")

    # Superbinder mutation side chains: sticks, distinct color.
    cmd.show("sticks", "pocket_mut_sel and not name C+N+O")
    cmd.set("stick_radius", 0.18, "pocket_mut_sel")
    cmd.color("pocket_mut_color", "pocket_mut_sel and elem C")
    cmd.set("cartoon_transparency", 0.0, "pocket_mut_sel")

    cmd.deselect()


def add_time_label_anchor():
    """Fixed-in-model-space label anchor (see gate_latch_movie.py for the
    reasoning -- same trick, no view-matrix decoding needed since the frame
    is drift-corrected and the camera never moves after set_fixed_camera).
    """
    core_model = cmd.get_model("core_fit_sel", state=1)
    peptide_model = cmd.get_model("peptide_sel", state=1)

    core_c = _centroid(core_model)
    peptide_c = _centroid(peptide_model)
    outward = _norm(_sub(peptide_c, core_c))

    peptide_radius = max(
        (
            (a.coord[0] - peptide_c[0]) ** 2
            + (a.coord[1] - peptide_c[1]) ** 2
            + (a.coord[2] - peptide_c[2]) ** 2
        )
        ** 0.5
        for a in peptide_model.atom
    )
    label_pos = _add(peptide_c, outward, peptide_radius + 14.0)

    anchor = "time_label_anchor"
    cmd.pseudoatom(anchor, pos=list(label_pos))
    cmd.hide("everything", anchor)
    cmd.set("label_size", 30, anchor)
    cmd.set("label_color", "black", anchor)
    cmd.set("label_font_id", 7, anchor)
    cmd.set("label_outline_color", "white", anchor)
    return anchor


def set_fixed_camera(anchor):
    """Orient + zoom exactly once, on state 1. No camera-moving command is
    called anywhere else in this script.
    """
    cmd.frame(1)
    cmd.orient("landmark_sel")
    cmd.zoom(f"landmark_sel or {anchor}", buffer=8)


def render_frames(n_states, frame_dir):
    cmd.viewport(VIEWPORT_W, VIEWPORT_H)

    for s in range(1, n_states + 1):
        cmd.frame(s)
        sim_ns = T0_NS + (s - 1) * STRIDE * TIME_SPACING_NS
        cmd.label("time_label_anchor", '"%.1f ns"' % sim_ns)

        png_path = os.path.join(frame_dir, f"frame_{s - 1:05d}.png")
        cmd.png(png_path, ray=0, quiet=1)

        if s == 1 or s % 25 == 0 or s == n_states:
            print(f"[sh2_movie] rendered frame {s}/{n_states} (t={sim_ns:.2f} ns)")


def stitch_movie(frame_dir, n_states):
    os.makedirs(OUT_DIR, exist_ok=True)
    ffmpeg_cmd = [
        FFMPEG_BIN,
        "-y",
        "-framerate", str(FPS),
        "-i", os.path.join(frame_dir, "frame_%05d.png"),
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", "18",
        "-movflags", "+faststart",
        OUT_MP4,
    ]
    print("[sh2_movie] running:", " ".join(ffmpeg_cmd))
    result = subprocess.run(ffmpeg_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(
            f"[sh2_movie] ffmpeg failed (exit {result.returncode}); "
            f"frames left in {frame_dir} for debugging."
        )

    size_bytes = os.path.getsize(OUT_MP4)
    duration_s = n_states / FPS
    print(
        f"[sh2_movie] wrote {OUT_MP4} "
        f"({size_bytes / 1e6:.1f} MB, {n_states} frames, "
        f"~{duration_s:.1f} s @ {FPS} fps)"
    )

    if not KEEP_FRAMES:
        shutil.rmtree(frame_dir, ignore_errors=True)
        print(f"[sh2_movie] cleaned up temp frame dir {frame_dir}")
    else:
        print(f"[sh2_movie] kept temp frame dir {frame_dir} (SH2_MOVIE_KEEP_FRAMES=1)")


def main():
    print(
        f"[sh2_movie] mode={'SMOKE TEST' if IS_SMOKE_TEST else 'FULL RUN'}, "
        f"stride={STRIDE}, n_output_frames={N_OUTPUT_FRAMES}, fps={FPS}, "
        f"out={OUT_MP4}"
    )

    setup_global_render_settings()
    n_states = load_topology_and_trajectory()
    define_selections()
    remove_rigid_body_drift()
    style_scene()
    anchor = add_time_label_anchor()
    set_fixed_camera(anchor)

    frame_dir = tempfile.mkdtemp(prefix="sh2_movie_frames_")
    print(f"[sh2_movie] rendering {n_states} frames into {frame_dir}")
    render_frames(n_states, frame_dir)
    stitch_movie(frame_dir, n_states)


main()
