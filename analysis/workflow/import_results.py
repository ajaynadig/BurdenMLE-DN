#!/usr/bin/env python3
"""Explicitly import approved results for a reporting-only rebuild, never implicitly."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import run as w


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",required=True)
    p.add_argument("--source",required=True)
    args=p.parse_args()
    c=w.load_config(args.config); source=Path(args.source).resolve(); root=Path(c["root"])
    if source==root.resolve(): raise ValueError("Source and destination must differ")
    w.initialize(c)
    stages=("models","simulations","upper","korean","public_example")
    if any(w.record_path(c,s).exists() for s in stages): raise RuntimeError("Imported stages already exist")
    # Reporting still reads raw inputs. Pin those inputs and the installed
    # environment just as for a fresh run, before copying expensive fit files.
    subprocess.run([c["rscript"],"--vanilla",str(w.HERE/"preflight.R")],env=w.environment(c),check=True)
    # Copy data products, not old reporting summaries. Model manifests are
    # relocated and validated below; no absolute reference to source fits remains.
    patterns=[f"outputs/data/models*_{c['run_id']}.Rdata",f"outputs/data/model_manifest_*_{c['run_id']}.rds",
              "outputs/derived/upper_rr_no_CES/*", "outputs/derived/korean_wgs/*",
              "outputs/tables/diagnostics/simulations/*", "outputs/tables/supplementary/SupplementaryTable1_*",
              "outputs/figures/main/Figure1.pdf", "outputs/figures/supplementary/SupplementaryFigure1_*",
              "outputs/figures/supplementary/SupplementaryFigure7_*", "outputs/figures/supplementary/SupplementaryFigure16_*",
              "qa/public_example/**/*"]
    copied={}
    for pattern in patterns:
        files=[f for f in source.glob(pattern) if f.is_file()]
        if not files: raise RuntimeError("Missing import input: "+pattern)
        for f in files:
            relative=f.relative_to(source); dest=root/relative
            dest.parent.mkdir(parents=True,exist_ok=True)
            if dest.exists(): raise RuntimeError("Refusing overwrite: "+str(dest))
            shutil.copy2(f,dest)
            if w.sha(dest)!=w.sha(f): raise RuntimeError("Copy checksum mismatch")
            copied[str(relative)]=w.sha(f)
    subprocess.run([c["rscript"],"--vanilla",str(w.HERE/"relocate_manifests.R")],env=w.environment(c),check=True)
    w.atomic_json(root/"provenance/imported_results.json",{"source":str(source),"files":copied,
                  "note":"Explicit result import; reporting rebuild, not fresh model fitting."})
    for stage in stages:
        outputs={str(f.relative_to(root)):w.sha(f) for f in w.products(stage,c)}
        w.atomic_json(w.record_path(c,stage),{"identity":w.sha(root/"workflow/identity.json"),"outputs":outputs,"imported":True})
    print("Imported verified products. Rebuild with --reuse "+",".join(stages))


if __name__=="__main__": main()
