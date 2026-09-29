# Step 1 
```
#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --output=output_blast.log
#SBATCH --account=acc_jfierst
#SBATCH --partition=highmem1-sapphirerapids
#SBATCH --qos=highmem1
#SBATCH --mem=384G
#SBATCH --cpus-per-task=32

# Load BLAST+ and set database
module load blast-plus
export BLASTDB='/home/data/jfierst/blastPractice/nt_db/'

# FASTA file to process
FASTA_FILE="JU3778_p_purged_new.fa"

# Check that the file exists
if [ ! -f "$FASTA_FILE" ]; then
    echo "FASTA file not found: $FASTA_FILE"
    exit 1
fi

# Output filename
OUTPUT_FILE="${FASTA_FILE%.fa}.out"

# Run BLAST
echo "Starting BLAST for file: $FASTA_FILE"

blastn \
    -query "$FASTA_FILE" \
    -db nt \
    -culling_limit 5 \
    -evalue 1e-25 \
    -num_threads 32 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -out "$OUTPUT_FILE"

echo "BLAST completed. Output saved to: $OUTPUT_FILE"
```

# Step 2
```
#!/bin/bash
#SBATCH --job-name=trim_contaminants
#SBATCH --output=trim_contaminants.log
#SBATCH --account=acc_jfierst
#SBATCH --partition=default-part
#SBATCH --qos=standard
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1


# ------------------------------------------------------------
# Step 2
#
# Conservative classification of BLAST hits.
#
# A BLAST hit is classified as contaminant only when:
#
#   1. Subject genus is in contaminants.txt
#   2. Percent identity >= MIN_PIDENT
#   3. Alignment length >= MIN_LENGTH
#   4. E-value <= MAX_EVALUE
#
# A PTG is placed in remove.txt only if it has at least one
# qualifying contaminant hit.
#
# IMPORTANT:
#   remove.txt identifies PTGs containing contaminant regions.
#   Step 3 uses the individual BLAST qstart/qend coordinates
#   to remove only those regions.
#
# Everything that is not confidently classified as contaminant
# goes into check.txt rather than remove.txt.
# ------------------------------------------------------------

SP="JU4645"
TYPE="p"

wd="/home/data/jfierst/Karolina/${SP}"
file="${wd}/${SP}.${TYPE}.fa.out"

keep="${wd}/${SP}_${TYPE}_keep.txt"
remove="${wd}/${SP}_${TYPE}_remove.txt"
check="${wd}/${SP}_${TYPE}_check.txt"

# Audit file containing the actual contaminant BLAST hits
contaminant_hits="${wd}/${SP}_${TYPE}_contaminant_hits.txt"

# ------------------------------------------------------------
# Conservative thresholds
# ------------------------------------------------------------

MIN_PIDENT=95
MIN_LENGTH=100
MAX_EVALUE="1e-25"

# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if [[ ! -s "$file" ]]; then
    echo "ERROR: Input BLAST file is missing or empty:"
    echo "$file" >&2
    exit 1
fi

# ------------------------------------------------------------
# Empty output files
# ------------------------------------------------------------

: > "$keep"
: > "$remove"
: > "$check"
: > "$contaminant_hits"

# ------------------------------------------------------------
# Classification
# ------------------------------------------------------------

awk -F '\t' \
    -v keep="$keep" \
    -v remove="$remove" \
    -v check="$check" \
    -v contaminant_hits="$contaminant_hits" \
    -v min_pident="$MIN_PIDENT" \
    -v min_length="$MIN_LENGTH" \
    -v max_evalue="$MAX_EVALUE" '

BEGIN {

    # --------------------------------------------------------
    # Read nematode genera
    # --------------------------------------------------------

    while ((getline line < "nematodes.txt") > 0) {

        gsub(/\r/, "", line)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)

        if (line != "")
            nem[line] = 1
    }

    close("nematodes.txt")


    # --------------------------------------------------------
    # Read contaminant genera
    # --------------------------------------------------------

    while ((getline line < "contaminants.txt") > 0) {

        gsub(/\r/, "", line)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)

        if (line != "")
            bac[line] = 1
    }

    close("contaminants.txt")
}


{

    # --------------------------------------------------------
    # BLAST columns
    #
    # $1  qseqid
    # $2  sseqid
    # $3  pident
    # $4  length
    # $5  mismatch
    # $6  gapopen
    # $7  qstart
    # $8  qend
    # $9  sstart
    # $10 send
    # $11 evalue
    # $12 bitscore
    # $13 stitle
    # --------------------------------------------------------

    ptg = $1
    pident = $3 + 0
    aln_length = $4 + 0
    evalue = $11 + 0
    description = $13


    # --------------------------------------------------------
    # Extract genus from subject description.
    #
    # Expected format:
    #
    # accession Genus species ...
    #
    # Example:
    #
    # CP115866.1 Caenorhabditis nigoni ...
    #
    # fields[2] = Caenorhabditis
    # --------------------------------------------------------

    split(description, fields, /[[:space:]]+/)

    genus = fields[2]


    # --------------------------------------------------------
    # Create unique key for this particular BLAST hit.
    #
    # Coordinates are included so that separate hits are
    # retained in the audit file.
    # --------------------------------------------------------

    hitkey = ptg "\t" \
             $2 "\t" \
             $7 "\t" \
             $8 "\t" \
             $11


    # --------------------------------------------------------
    # Strong contaminant hit
    # --------------------------------------------------------

    if ((genus in bac) &&
        pident >= min_pident &&
        aln_length >= min_length &&
        evalue <= max_evalue) {

        # ----------------------------------------------------
        # Record PTG as containing contaminant sequence.
        # ----------------------------------------------------

        contaminant_ptg[ptg] = 1

        # ----------------------------------------------------
        # Record individual contaminant hit.
        #
        # This file is useful for auditing exactly what Step 3
        # will remove.
        # ----------------------------------------------------

        if (!(hitkey in seen_contaminant)) {

            print ptg "\t" \
                  $2 "\t" \
                  pident "\t" \
                  aln_length "\t" \
                  $7 "\t" \
                  $8 "\t" \
                  evalue "\t" \
                  $12 "\t" \
                  genus "\t" \
                  description \
                  >> contaminant_hits

            seen_contaminant[hitkey] = 1
        }

        next
    }


    # --------------------------------------------------------
    # Nematode hit
    #
    # A strong nematode hit is retained in keep.txt.
    #
    # We do NOT use this to cancel a contaminant hit. A PTG can
    # legitimately contain both nematode and contaminant
    # sequence. Step 3 will use coordinates to retain the
    # nematode portion.
    # --------------------------------------------------------

    if (genus in nem) {

        key = ptg "\t" description

        if (!(key in seen_keep)) {

            print ptg "\t" description >> keep

            seen_keep[key] = 1
        }

        next
    }


    # --------------------------------------------------------
    # Everything else goes to check.txt.
    #
    # This includes:
    #
    #   - weak contaminant hits
    #   - short contaminant hits
    #   - lower identity contaminant hits
    #   - unknown genera
    #   - poorly formatted descriptions
    # --------------------------------------------------------

    key = ptg "\t" description

    if (!(key in seen_check)) {

        print ptg "\t" description >> check

        seen_check[key] = 1
    }
}


END {

    # --------------------------------------------------------
    # Write remove.txt only after all BLAST records have been
    # examined.
    #
    # This guarantees that a PTG is listed once.
    # --------------------------------------------------------

    for (ptg in contaminant_ptg) {

        print ptg "\tstrong contaminant BLAST hit" >> remove
    }
}

' "$file"


# ------------------------------------------------------------
# Sort output files for reproducibility
# ------------------------------------------------------------

sort -u "$keep" > "${keep}.tmp"
mv "${keep}.tmp" "$keep"

sort -u "$remove" > "${remove}.tmp"
mv "${remove}.tmp" "$remove"

sort -u "$check" > "${check}.tmp"
mv "${check}.tmp" "$check"


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

N_KEEP=$(cut -f1 "$keep" | sort -u | wc -l)
N_REMOVE=$(cut -f1 "$remove" | sort -u | wc -l)
N_CHECK=$(cut -f1 "$check" | sort -u | wc -l)
N_CONTAMINANT_HITS=$(wc -l < "$contaminant_hits")

echo
echo "=============================================="
echo "Conservative BLAST classification complete"
echo "=============================================="
echo
echo "Input:"
echo "  $file"
echo
echo "Thresholds:"
echo "  Minimum identity : ${MIN_PIDENT}%"
echo "  Minimum length   : ${MIN_LENGTH} bp"
echo "  Maximum e-value  : ${MAX_EVALUE}"
echo
echo "PTGs:"
echo "  Keep              : $N_KEEP"
echo "  Contaminant       : $N_REMOVE"
echo "  Check             : $N_CHECK"
echo
echo "Strong contaminant BLAST hits:"
echo "  $N_CONTAMINANT_HITS"
echo
echo "Output:"
echo "  Keep              : $keep"
echo "  Remove            : $remove"
echo "  Check             : $check"
echo "  Contaminant hits  : $contaminant_hits"
echo
```

# Step 3

```
#!/bin/bash
#SBATCH --job-name=trim_contaminants
#SBATCH --output=trim_contaminants.log
#SBATCH --account=acc_jfierst
#SBATCH --partition=default-part
#SBATCH --qos=standard
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

set -euo pipefail

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

SPECIES="JU4645"
TYPE="p"

WD="/home/data/jfierst/Karolina/${SPECIES}"

INPUT_FASTA="${WD}/${SPECIES}.${TYPE}.fa"
BLAST_FILE="${WD}/${SPECIES}.${TYPE}.fa.out"
REMOVE_FILE="${WD}/${SPECIES}_${TYPE}_remove.txt"

OUT_FILE="${WD}/${SPECIES}_${TYPE}_decontaminated.fa"

INTERVAL_FILE="${WD}/${SPECIES}_${TYPE}_contaminant_intervals.tmp"
MERGED_FILE="${WD}/${SPECIES}_${TYPE}_merged_intervals.tmp"
CONTAMINANT_PTG="${WD}/${SPECIES}_${TYPE}_contaminant_ptgs.tmp"

# ------------------------------------------------------------
# Check required files
# ------------------------------------------------------------

echo "Running contaminant trimming"
echo
echo "  Input FASTA      : $INPUT_FASTA"
echo "  BLAST output     : $BLAST_FILE"
echo "  Remove file      : $REMOVE_FILE"
echo "  Output FASTA     : $OUT_FILE"
echo

if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA not found or empty:"
    echo "$INPUT_FASTA"
    exit 1
fi
if [[ ! -s "$BLAST_FILE" ]]; then
    echo "ERROR: BLAST output not found or empty:"
    echo "$BLAST_FILE"
    exit 1
fi

if [[ ! -s "$REMOVE_FILE" ]]; then
    echo "ERROR: Remove file not found or empty:"
    echo "$REMOVE_FILE"
    exit 1
fi

# ------------------------------------------------------------
# Step 1
#
# Extract PTG IDs classified as contaminants.
#
# remove.txt format:
#
# PTG_ID    description
#
# Only column 1 is needed.
# ------------------------------------------------------------

cut -f1 "$REMOVE_FILE" |
    sed '/^[[:space:]]*$/d' |
    sort -u > "$CONTAMINANT_PTG"

N_CONTAMINANT_PTG=$(wc -l < "$CONTAMINANT_PTG")

echo "Contaminant PTGs: $N_CONTAMINANT_PTG"

# ------------------------------------------------------------
# Step 2
#
# Extract contaminant BLAST coordinates.
#
# BLAST outfmt columns:
#
#   1  qseqid
#   2  sseqid
#   3  pident
#   4  length
#   5  mismatch
#   6  gapopen
#   7  qstart
#   8  qend
#   9  sstart
#   10 send
#   11 evalue
#   12 bitscore
#   13 stitle
#
# qstart/qend are coordinates on the PTG/query.
#
# These coordinates are 1-based and inclusive.
# ------------------------------------------------------------

awk '
    BEGIN {
        FS = "\t"
        OFS = "\t"

        while ((getline ptg < "'"$CONTAMINANT_PTG"'") > 0) {
            gsub(/\r/, "", ptg)

            if (ptg != "")
                contaminant[ptg] = 1
        }

        close("'"$CONTAMINANT_PTG"'")
    }

    NF >= 8 && ($1 in contaminant) {

        start = $7
        end   = $8

        # Normalize orientation.
        if (start <= end) {
            s = start
            e = end
        }
        else {
            s = end
            e = start
        }
        # Ignore malformed coordinates.
        if (s ~ /^[0-9]+$/ &&
            e ~ /^[0-9]+$/ &&
            s > 0 &&
            e >= s) {

            print $1, s, e
        }
    }
' "$BLAST_FILE" > "$INTERVAL_FILE"

N_INTERVALS=$(wc -l < "$INTERVAL_FILE")

echo "Contaminant BLAST intervals: $N_INTERVALS"

# ------------------------------------------------------------
# Step 3
#
# Merge overlapping or adjacent contaminant intervals.
#
# Example:
#
#   PTG1  100 200
#   PTG1  180 250
#   PTG1  251 300
#
# becomes:
#
#   PTG1  100 300
#
# This prevents overlapping BLAST hits from causing problems
# during trimming.
# ------------------------------------------------------------
sort -k1,1 -k2,2n -k3,3n "$INTERVAL_FILE" |
awk '
    BEGIN {
        FS = "\t"
        OFS = "\t"
    }

    {
        ptg = $1
        start = $2
        end = $3

        if (!have) {
            current_ptg = ptg
            current_start = start
            current_end = end
            have = 1
            next
        }

        # Same PTG and overlapping/adjacent interval.
        if (ptg == current_ptg &&
            start <= current_end + 1) {

            if (end > current_end)
                current_end = end
        }
        else {

            print current_ptg, current_start, current_end

            current_ptg = ptg
            current_start = start
            current_end = end
        }
    }

    END {
        if (have)
            print current_ptg, current_start, current_end
    }
' > "$MERGED_FILE"

N_MERGED=$(wc -l < "$MERGED_FILE")
echo "Merged contaminant intervals: $N_MERGED"

# ------------------------------------------------------------
# Step 4
#
# Read the original FASTA and remove contaminant intervals.
#
# IMPORTANT:
#
# The entire PTG is NOT removed.
#
# Only the sequence between qstart and qend is removed.
#
# If a PTG has:
#
#   1-1000       contaminant
#   1001-5000    nematode
#
# only bases 1-1000 are removed.
#
# If the contaminant covers the entire PTG, the PTG is omitted
# from the output because no sequence remains.
# ------------------------------------------------------------

: > "$OUT_FILE"

awk '
    BEGIN {
        FS = "\t"

        # ----------------------------------------------------
        # Read merged contaminant intervals.
        # ----------------------------------------------------

        while ((getline line < "'"$MERGED_FILE"'") > 0) {

            split(line, a, "\t")

            ptg = a[1]
            start = a[2] + 0
            end = a[3] + 0
            n[ptg]++

            interval_start[ptg, n[ptg]] = start
            interval_end[ptg, n[ptg]] = end
        }

        close("'"$MERGED_FILE"'")

        current_id = ""
        sequence = ""

        processed = 0
        trimmed = 0
        unchanged = 0
        removed = 0
    }

    # --------------------------------------------------------
    # FASTA header
    # --------------------------------------------------------

    /^>/ {

        # Process previous sequence.
        if (current_id != "") {
            process_sequence(current_id, sequence)
        }

        current_id = $1
        sub(/^>/, "", current_id)

        sequence = ""

        next
    }

    # --------------------------------------------------------
    # FASTA sequence
    # --------------------------------------------------------
 {
        gsub(/[[:space:]]/, "", $0)
        sequence = sequence $0
    }

    # --------------------------------------------------------
    # End of FASTA
    # --------------------------------------------------------

    END {

        if (current_id != "") {
            process_sequence(current_id, sequence)
        }

        print ""
        print "FASTA records processed:", processed > "/dev/stderr"
        print "PTGs trimmed:", trimmed > "/dev/stderr"
        print "PTGs unchanged:", unchanged > "/dev/stderr"
        print "PTGs completely removed:", removed > "/dev/stderr"
    }

    # --------------------------------------------------------
    # Remove contaminant intervals from one PTG.
    # --------------------------------------------------------

    function process_sequence(id, seq,
                              i, start, end,
                              pos, result,
                              seq_length) {

        processed++

        seq_length = length(seq)
        # ----------------------------------------------------
        # No contaminant hit for this PTG.
        #
        # Keep the original sequence unchanged.
        # ----------------------------------------------------

        if (!(id in n)) {

            unchanged++

            print ">" id
            print_wrapped(seq)

            return
        }

        result = ""
        pos = 1

        # ----------------------------------------------------
        # Remove each contaminant interval.
        # ----------------------------------------------------

        for (i = 1; i <= n[id]; i++) {

            start = interval_start[id, i]
            end   = interval_end[id, i]

            # Keep coordinates inside the actual sequence.
            if (start < 1)
                start = 1

            if (end > seq_length)
                end = seq_length

            if (start > seq_length)
                continue

            # Keep sequence before contaminant region.
            if (start > pos) {
                result = result substr(seq, pos, start - pos)
            }

            # Move past contaminant region.
            pos = end + 1
        }

        # ----------------------------------------------------
        # Keep sequence after final contaminant region.
        # ----------------------------------------------------
           if (pos <= seq_length) {
            result = result substr(seq, pos)
        }

        # ----------------------------------------------------
        # Write remaining sequence.
        # ----------------------------------------------------

        if (length(result) == 0) {

            # Entire PTG was contaminant.
            removed++

        }
        else {

            trimmed++

            print ">" id
            print_wrapped(result)
        }
    }

    # --------------------------------------------------------
    # Write FASTA sequence at 60 bp per line.
    # --------------------------------------------------------

    function print_wrapped(seq, i, width) {

        width = 60

        for (i = 1; i <= length(seq); i += width) {
            print substr(seq, i, width)
        }
    }

' "$INPUT_FASTA" > "$OUT_FILE"

# ------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------
rm -f "$CONTAMINANT_PTG"
rm -f "$INTERVAL_FILE"
rm -f "$MERGED_FILE"

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "============================================"
echo "Contaminant trimming complete"
echo "============================================"
echo
echo "Input FASTA : $INPUT_FASTA"
echo "BLAST file  : $BLAST_FILE"
echo "Remove file : $REMOVE_FILE"
echo "Output FASTA: $OUT_FILE"
echo
echo "Contaminant PTGs       : $N_CONTAMINANT_PTG"
echo "BLAST intervals        : $N_INTERVALS"
echo "Merged intervals       : $N_MERGED"
echo
echo "Only contaminant regions"
echo "identified by BLAST qstart/qend"
echo "were removed."
echo
```
