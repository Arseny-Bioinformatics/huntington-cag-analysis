"""
huntingtons_cag_analysis.py
────────────────────────────────────────────────────────────────────
Bioinformatics script to detect and categorize CAG repeat expansions
in the HTT gene — a key diagnostic marker for Huntington's Disease.

Author  : Biology Student
Purpose : Learn Python + Bioinformatics fundamentals
Concepts: String parsing, regex, matplotlib visualisation, pandas scaling
────────────────────────────────────────────────────────────────────
"""

import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ─── 1. CONSTANTS ────────────────────────────────────────────────────────────
# Clinical thresholds are defined by the Huntington's Disease Society of America.
# These numbers represent the number of CAG triplet repeats in the HTT gene.

THRESHOLD_NORMAL       = 35   # ≤35  → Normal (unaffected)
THRESHOLD_INTERMEDIATE = 39   # 36–39 → Intermediate (reduced penetrance)
# >40 repeats           → Affected (full penetrance / HD diagnosis)

TRIPLET = "CAG"               # The nucleotide motif we are searching for


# ─── 2. CORE SEQUENCE ANALYSIS ───────────────────────────────────────────────

def find_longest_cag_run(dna_sequence: str) -> int:
    """
    Find the longest *consecutive* run of CAG triplets in a DNA string.

    How it works
    ────────────
    Regular expression breakdown:
        (?:{TRIPLET})+    →  one or more non-capturing groups of "CAG"

    re.findall() returns ALL non-overlapping matches of this pattern.
    We then count how many times "CAG" appears inside the longest match,
    which gives us the number of consecutive repeats.

    Parameters
    ──────────
    dna_sequence : str
        Raw DNA string (A, T, C, G). Case-insensitive.

    Returns
    ───────
    int : Length of the longest CAG repeat run (0 if none found).
    """
    # Normalise to uppercase so "cag" and "CAG" are treated equally
    dna_sequence = dna_sequence.upper().strip()

    # Build regex pattern: one or more consecutive CAG blocks
    pattern = f"(?:{TRIPLET})+"

    # Find all consecutive CAG runs in the sequence
    matches = re.findall(pattern, dna_sequence)

    if not matches:
        return 0

    # The longest match string divided by 3 gives the repeat count
    longest_match = max(matches, key=len)
    repeat_count  = len(longest_match) // len(TRIPLET)

    return repeat_count


def categorize_result(repeat_count: int) -> tuple[str, str]:
    """
    Translate a CAG repeat count into a clinical category.

    Parameters
    ──────────
    repeat_count : int  — number of consecutive CAG repeats

    Returns
    ───────
    (category_label, hex_color) : tuple for display and visualisation
    """
    if repeat_count <= THRESHOLD_NORMAL:
        return "Normal",       "#2ecc71"   # Green
    elif repeat_count <= THRESHOLD_INTERMEDIATE:
        return "Intermediate", "#f39c12"   # Amber
    else:
        return "Affected",     "#e74c3c"   # Red


# ─── 3. VISUALISATION ────────────────────────────────────────────────────────

def plot_cag_chart(repeat_count: int, category: str, color: str) -> None:
    """
    Draw a bar chart comparing the patient's CAG count against clinical
    threshold lines for Normal, Intermediate, and Affected ranges.

    Parameters
    ──────────
    repeat_count : int   — patient's longest CAG run
    category     : str   — clinical label ("Normal" / "Intermediate" / "Affected")
    color        : str   — hex color corresponding to the category
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("#1a1a2e")       # Dark background
    ax.set_facecolor("#16213e")

    # ── Bar: patient result ──────────────────────────────────────────────────
    bar = ax.bar(
        x       = ["Patient CAG Repeats"],
        height  = [repeat_count],
        color   = color,
        width   = 0.4,
        zorder  = 3,
        edgecolor = "white",
        linewidth = 0.8,
    )

    # Annotate the bar with the numeric count
    ax.text(
        0, repeat_count + 0.5,
        str(repeat_count),
        ha="center", va="bottom",
        fontsize=18, fontweight="bold", color="white"
    )

    # ── Threshold lines ──────────────────────────────────────────────────────
    thresholds = [
        (THRESHOLD_NORMAL,       "#2ecc71", "Normal / Intermediate boundary (35)"),
        (THRESHOLD_INTERMEDIATE, "#f39c12", "Intermediate / Affected boundary (39)"),
    ]
    for y_val, line_color, label in thresholds:
        ax.axhline(y=y_val, color=line_color, linewidth=1.5,
                   linestyle="--", alpha=0.8, label=label, zorder=2)

    # ── Shaded regions ───────────────────────────────────────────────────────
    y_max = max(repeat_count + 10, 55)
    ax.axhspan(0,                        THRESHOLD_NORMAL,       alpha=0.06, color="#2ecc71")
    ax.axhspan(THRESHOLD_NORMAL + 1,     THRESHOLD_INTERMEDIATE, alpha=0.06, color="#f39c12")
    ax.axhspan(THRESHOLD_INTERMEDIATE + 1, y_max,                alpha=0.06, color="#e74c3c")

    # ── Axis styling ─────────────────────────────────────────────────────────
    ax.set_ylim(0, y_max)
    ax.set_ylabel("CAG Repeat Count", color="white", fontsize=12)
    ax.set_title(
        f"HTT Gene — CAG Repeat Analysis\nResult: {category}  ({repeat_count} repeats)",
        color="white", fontsize=14, fontweight="bold", pad=15
    )
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#444")

    # ── Legend ───────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color="#2ecc71", label="Normal (≤35)"),
        mpatches.Patch(color="#f39c12", label="Intermediate (36–39)"),
        mpatches.Patch(color="#e74c3c", label="Affected (>40)"),
    ]
    ax.legend(
        handles=legend_patches, loc="upper left",
        facecolor="#1a1a2e", edgecolor="#444",
        labelcolor="white", fontsize=9
    )

    plt.tight_layout()
    plt.savefig("cag_analysis_chart.png", dpi=150, bbox_inches="tight")
    print("  → Chart saved as 'cag_analysis_chart.png'")
    plt.show()


# ─── 4. MAIN PIPELINE ────────────────────────────────────────────────────────

def analyze_patient(dna_sequence: str, patient_id: str = "Patient") -> dict:
    """
    Full analysis pipeline for a single DNA sequence.

    Steps
    ─────
    1. Find the longest CAG run  (find_longest_cag_run)
    2. Categorize clinically     (categorize_result)
    3. Print a summary report
    4. Plot the result           (plot_cag_chart)

    Returns
    ───────
    dict with keys: patient_id, repeat_count, category
    """
    print(f"\n{'═' * 50}")
    print("  Huntington's Disease — CAG Repeat Analysis")
    print(f"  Patient : {patient_id}")
    print(f"{'═' * 50}")

    # Step 1 – Count
    count = find_longest_cag_run(dna_sequence)
    print(f"  Longest consecutive CAG run : {count} repeats")

    # Step 2 – Categorise
    category, color = categorize_result(count)
    print(f"  Clinical category           : {category}")

    # Step 3 – Clinical note
    notes = {
        "Normal":       "No elevated risk detected.",
        "Intermediate": "Reduced penetrance — genetic counselling recommended.",
        "Affected":     "Full penetrance — consistent with HD diagnosis.",
    }
    print(f"  Clinical note               : {notes[category]}")
    print(f"{'─' * 50}\n")

    # Step 4 – Chart
    plot_cag_chart(count, category, color)

    return {"patient_id": patient_id, "repeat_count": count, "category": category}


# ─── 5. SCALING WITH PANDAS (demonstration) ──────────────────────────────────

def analyze_cohort_with_pandas(csv_filepath: str) -> None:
    """
    ─────────────────────────────────────────────────────────────
    HOW TO SCALE TO 1,000 PATIENTS USING PANDAS
    ─────────────────────────────────────────────────────────────

    Expected CSV format:
        patient_id,dna_sequence
        P001,ATGATGCAGCAGCAGCAGCAGCAGCAG...
        P002,GCATCAGCAGCAGCAGCAGCAGCAGCAG...

    This function shows the pattern — it will work once you have
    a real CSV file at `csv_filepath`.
    ─────────────────────────────────────────────────────────────
    """
    import pandas as pd

    # ── Load ──────────────────────────────────────────────────────────────────
    df = pd.read_csv(csv_filepath)                          # Read all rows at once

    # ── Vectorised analysis ───────────────────────────────────────────────────
    # .apply() runs our function on every row — no manual loop needed
    df["repeat_count"] = df["dna_sequence"].apply(find_longest_cag_run)

    # Map the count to a category label using the same logic
    df["category"] = df["repeat_count"].apply(
        lambda n: categorize_result(n)[0]   # [0] = label, [1] = color
    )

    # ── Summary statistics ────────────────────────────────────────────────────
    print("\n── Cohort Summary ──────────────────────────────")
    print(df["category"].value_counts().to_string())
    print(f"\nMean CAG repeats : {df['repeat_count'].mean():.1f}")
    print(f"Max  CAG repeats : {df['repeat_count'].max()}")

    # ── Export enriched results ───────────────────────────────────────────────
    df.to_csv("cohort_results.csv", index=False)
    print("  → Results saved to 'cohort_results.csv'")
        # ── Population-level bar chart ────────────────────────────────────────────
    

    # ── Population-level bar chart ────────────────────────────────────────────
# ── Population-level bar chart ────────────────────────────────────────────

# ─── 6. ENTRY POINT ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── Example sequences ─────────────────────────────────────────────────────
    # Modify these or replace with real FASTA data.

    # Simulated HTT exon-1 region with 42 CAG repeats  → Affected
    sequence_affected = (
        "ATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAG"
        "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
        "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
        "CAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCG"
    )

    # Simulated sequence with 20 CAG repeats  → Normal
    sequence_normal = (
        "ATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAG"
        "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCG"
    )

    # Run analysis for the affected sample
    #result = analyze_patient(sequence_affected, patient_id="Sample_HD_001")

    # Uncomment the next line to test the normal sample:
    #result = analyze_patient(sequence_normal, patient_id="Sample_CTL_001")

    # Uncomment to run the Pandas cohort analysis (requires a CSV file):
    analyze_cohort_with_pandas("huntington_demo_cohort_30.csv")