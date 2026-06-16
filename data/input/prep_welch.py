"""
Extract Welch 2009 data from .xls → clean CSV + combined FASTA.
Run once: python prep_welch.py
"""
import xlrd, os, sys

sys.stdout.reconfigure(encoding="utf-8")
base = os.path.dirname(os.path.abspath(__file__))
xls_path = os.path.join(base, "..", "sequences", "welch_data.xls")
csv_path = os.path.join(base, "welch_variants.csv")
fasta_out = os.path.join(base, "..", "sequences", "Ecoli_K12_MG1655_geneset_with_welch_2009_welch.fa")
bg_fasta  = os.path.join(base, "..", "sequences", "Ecoli_K12_MG1655_geneset_with_welch_2009.fa")

wb = xlrd.open_workbook(xls_path)
ws = wb.sheet_by_index(0)

# MW estimates (Da) from sequence length
MW_SCFV = 281 * 110.0       # 843 bp / 3 = 281 aa
MW_POLY = 575 * 110.0       # 1725 bp / 3 = 575 aa

rows = []
for r in range(1, ws.nrows):
    variant = str(ws.cell_value(r, 0)).strip()
    seq = str(ws.cell_value(r, 16)).strip().upper()
    if len(seq) < 100:
        continue
    expr_raw = ws.cell_value(r, 2)
    if isinstance(expr_raw, float):
        expr_mg_mL = expr_raw
    elif isinstance(expr_raw, str) and expr_raw.strip().lower() == "nd":
        expr_mg_mL = None
    else:
        expr_mg_mL = None

    cai_raw = ws.cell_value(r, 8)
    cai = float(cai_raw) if isinstance(cai_raw, (int, float)) else None

    ptype = "scFv" if len(seq) < 1000 else "Polymerase"
    mw = MW_SCFV if ptype == "scFv" else MW_POLY
    # Welch 2009 "Absolute Expression, mg/ml" is almost certainly µg/mL (mg/L).
    # 206 mg/mL would exceed total E. coli protein; 206 µg/mL at 30% expression is realistic.
    expr_ug_mL = expr_mg_mL  # keep raw column name for traceability
    expr_uM = (expr_ug_mL * 1e-3 * 1e6 / mw) if expr_ug_mL is not None else None  # µg/mL → mg/mL → µM

    rows.append({
        "variant": variant,
        "protein": ptype,
        "seq_len_bp": len(seq),
        "expr_mg_mL": expr_mg_mL,
        "expr_uM": expr_uM,
        "CAI": cai,
        "sequence": seq,
    })

# Write CSV
with open(csv_path, "w") as f:
    f.write("variant,protein,seq_len_bp,expr_mg_mL,expr_uM,CAI,sequence\n")
    for row in rows:
        e1 = f"{row['expr_mg_mL']:.4f}" if row["expr_mg_mL"] is not None else ""
        e2 = f"{row['expr_uM']:.4f}" if row["expr_uM"] is not None else ""
        c = f"{row['CAI']:.6f}" if row["CAI"] is not None else ""
        f.write(f"{row['variant']},{row['protein']},{row['seq_len_bp']},{e1},{e2},{c},{row['sequence']}\n")

print(f"Wrote {len(rows)} variants to {csv_path}")
n_expr = sum(1 for r in rows if r["expr_mg_mL"] is not None)
print(f"  {n_expr} with measured expression, {len(rows) - n_expr} without (nd)")

# Create combined FASTA: copy background genes + append Welch variants
with open(bg_fasta, "r") as f:
    bg_content = f.read()

with open(fasta_out, "w") as f:
    f.write(bg_content)
    if not bg_content.endswith("\n"):
        f.write("\n")
    for row in rows:
        f.write(f">{row['variant']}\n")
        seq = row["sequence"]
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

print(f"Wrote combined FASTA to {fasta_out}")
print(f"  Background genes + {len(rows)} Welch variants")
