#!/usr/bin/env python3
"""Coordinates are not in SC_dat.csv. They are in the segmentation masks.

An IMC mask is a LABEL IMAGE: the pixel value IS the CellId for that core.
regionprops gives (label, centroid) = (CellId, y, x) directly. No inference.

OMEandSingleCellMasks.zip is 11.4 GB, but that is dominated by multi-channel
OME-TIFFs we never open. The masks live in a NESTED zip member
(OMEnMasks/Basel_Zuri_masks.zip, ~33 MB). Range-fetch only that nested zip and
open it in memory; ome.zip is never touched.

Filename -> core mapping (verified against Basel_PatientMetadata.csv):
  mask stem ...20170905_121_100_X15Y5_252_a0_full_maks  ->  core BaselTMA_SP41_100_X15Y5
The core id is {TMA}_{SP}_{N}_{XkYk} where N is the number IMMEDIATELY BEFORE the
XkYk token, NOT the trailing number. About 24 cores carry a trailing disambiguator
(e.g. ..._X8Y4_175) when (N, XkYk) is not unique; those append the mask's trailing
number. We snap each candidate to the real core set so the join can never silently
drop a core.
"""
import io, re, struct, sys, zipfile, zlib
import numpy as np, pandas as pd, requests, tifffile
from remotezip import RemoteZip
from skimage.measure import regionprops_table

URL = ("https://zenodo.org/records/3518284/files/"
       "OMEandSingleCellMasks.zip?download=1")
OUT = "data/coords.csv"
META = "data/Basel_PatientMetadata.csv"  # core universe, to disambiguate suffixed cores


def _toc(url):
    with RemoteZip(url) as z:
        return z.infolist()


def _read_member(url, info, sess):
    """Exact compressed byte span. Reading the local header is not optional:
    the central directory's extra-field length can differ from the local one."""
    off = info.header_offset
    lh = sess.get(url, headers={"Range": f"bytes={off}-{off+29}"}).content
    if lh[:4] != b"PK\x03\x04":
        raise RuntimeError("bad local header; server ignoring Range?")
    n, m = struct.unpack("<HH", lh[26:30])
    start = off + 30 + n + m
    body = sess.get(url, headers={"Range": f"bytes={start}-{start+info.compress_size-1}"}).content
    return body if info.compress_type == 0 else zlib.decompress(body, -15)


def _mask_zip(sess):
    """Open the nested Basel_Zuri_masks.zip in memory (range-fetched, ~33 MB)."""
    outer = {i.filename: i for i in _toc(URL)}
    nested = next(n for n in outer if n.lower().endswith(".zip") and "mask" in n.lower())
    zf = zipfile.ZipFile(io.BytesIO(_read_member(URL, outer[nested], sess)))
    return zf, outer[nested].compress_size


def _mask_members(zf, cohort="basel"):
    prefix = {"basel": "BaselTMA", "zurich": "ZTMA"}[cohort]
    return [m for m in zf.infolist()
            if m.filename.lower().endswith((".tif", ".tiff"))
            and m.filename.split("/")[-1].startswith(prefix)]


def _candidates(name):
    """(v1, v2) core-id candidates. v1 = {TMA}_{SP}_{N}_{XkYk} with N the number
    before the XkYk token; v2 appends the mask's trailing number (the disambiguator)."""
    stem = re.sub(r"_a0_full_ma[sk]{2}$", "", name.split("/")[-1].rsplit(".", 1)[0])
    toks = stem.split("_")
    xy_i = next((i for i, t in enumerate(toks) if re.fullmatch(r"X\d+Y\d+", t)), None)
    if xy_i is None or xy_i < 2 or not re.fullmatch(r"\d+", toks[xy_i - 1]):
        return None, None  # e.g. Liver control cores have no numeric N before XY
    v1 = f"{toks[0]}_{toks[1]}_{toks[xy_i - 1]}_{toks[xy_i]}"
    trailing = toks[xy_i + 1] if xy_i + 1 < len(toks) and re.fullmatch(r"\d+", toks[xy_i + 1]) else None
    return v1, (f"{v1}_{trailing}" if trailing else v1)


def core_from_filename(name, valid=None):
    """Map a mask filename to the SC_dat `core`. VERIFY against inspect output.
    Wrong mapping = the join silently drops every cell."""
    v1, v2 = _candidates(name)
    if v1 is None:
        return None
    if valid is None:
        return v1
    if v1 in valid:
        return v1
    if v2 in valid:
        return v2
    return None


def _valid_cores():
    try:
        return set(pd.read_csv(META)["core"].astype(str))
    except Exception:
        return None


def cmd_inspect():
    sess = requests.Session()
    zf, wire = _mask_zip(sess)
    basel, zurich = _mask_members(zf, "basel"), _mask_members(zf, "zurich")
    valid = _valid_cores()
    print(f"nested masks zip: {wire/1e6:.0f} MB on the wire")
    print(f"{len(basel) + len(zurich)} masks  (Basel {len(basel)}, Zurich {len(zurich)})\n")
    for m in basel[:5]:
        print(f"  {m.filename.split('/')[-1]}")
        print(f"      -> core = {core_from_filename(m.filename, valid)!r}")
    if valid is not None:
        hit = sum(1 for m in basel if core_from_filename(m.filename, valid) in valid)
        covered = len({core_from_filename(m.filename, valid) for m in basel} & valid)
        print(f"\nBasel masks -> a real core: {hit}/{len(basel)}  "
              f"(cores covered: {covered}/{len(valid)})")
    print("\nCHECK: derived cores must match SC_dat's `core` (e.g. BaselTMA_SP41_257_X3Y1).")


def cmd_run():
    sess = requests.Session()
    zf, _ = _mask_zip(sess)
    valid = _valid_cores()
    basel = _mask_members(zf, "basel")
    print(f"pulling {len(basel)} Basel masks from the in-memory nested zip")
    frames, dropped, failed = [], [], []
    for k, info in enumerate(basel, 1):
        core = core_from_filename(info.filename, valid)
        if core is None:
            dropped.append(info.filename.split("/")[-1])
            continue
        try:
            img = np.squeeze(tifffile.imread(io.BytesIO(zf.read(info.filename))))
            t = regionprops_table(img, properties=("label", "centroid"))
        except Exception as exc:
            failed.append((info.filename.split("/")[-1], str(exc)))
            continue
        frames.append(pd.DataFrame({
            "core": core,
            "CellId": t["label"].astype(np.int64),
            # regionprops centroid is (row, col) == (y, x). Swapping these silently
            # transposes every tissue map. Do not "fix" this.
            "x": t["centroid-1"],
            "y": t["centroid-0"],
        }))
        print(f"\r  {k}/{len(basel)}  {frames[-1].shape[0]} cells", end="", flush=True)
    coords = pd.concat(frames, ignore_index=True)
    coords["id"] = coords.core + "_" + coords.CellId.astype(str)
    coords.to_csv(OUT, index=False)
    print(f"\n\n{OUT}: {len(coords):,} cells, {coords.core.nunique()} cores")
    if valid is not None:
        print(f"cores in coords & metadata = {len(set(coords.core) & valid)}/{len(valid)}")
    if dropped:
        print(f"dropped {len(dropped)} non-mappable masks (liver controls etc.): {dropped[:3]}")
    if failed:
        print(f"FAILED to read {len(failed)} masks: {failed[:3]}")
    print(coords.head())
    print("\nNEXT: join to SC_dat.csv on `id`. Expect ~1.2M cells for Basel.")
    print("If the join drops most rows, core_from_filename is wrong.")


if __name__ == "__main__":
    {"inspect": cmd_inspect, "run": cmd_run}[sys.argv[1] if len(sys.argv) > 1 else "inspect"]()
