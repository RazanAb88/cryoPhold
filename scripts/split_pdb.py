#!/usr/bin/env python3
import argparse
import os

def split_models(pdb_path, outdir):
    """
    Splits a multi‐MODEL PDB into separate files per MODEL.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    model_lines = []
    current_model = None
    basename = os.path.splitext(os.path.basename(pdb_path))[0]

    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                # If we were already collecting one model, write it out
                if current_model is not None:
                    write_model_file(basename, current_model, model_lines, outdir)
                    model_lines = []
                # Get the model index from the record (e.g. MODEL   1)
                parts = line.split()
                current_model = parts[1] if len(parts) > 1 else str(len(os.listdir(outdir)) + 1)
                model_lines.append(line)
            elif line.startswith("ENDMDL"):
                model_lines.append(line)
                # finish this model
                write_model_file(basename, current_model, model_lines, outdir)
                model_lines = []
                current_model = None
            else:
                # if we're between MODEL/ENDMDL, collect lines
                if current_model is not None:
                    model_lines.append(line)

    # In case the file has no MODEL/ENDMDL but only ATOM records
    if current_model is None and not os.listdir(outdir):
        # just copy the whole thing as model1
        write_model_file(basename, "1", open(pdb_path).read().splitlines(True), outdir)

def write_model_file(basename, model_id, lines, outdir):
    fname = f"{basename}_model{model_id}.pdb"
    outpath = os.path.join(outdir, fname)
    with open(outpath, "w") as out:
        out.writelines(lines)
    print(f"Wrote {outpath}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Split a multi‐MODEL PDB into one file per MODEL block."
    )
    p.add_argument("-i", "--input", required=True,
                   help="Path to input multi‐model PDB")
    p.add_argument("-o", "--outdir", default="models_split",
                   help="Directory to write individual PDBs")
    args = p.parse_args()

    split_models(args.input, args.outdir)
