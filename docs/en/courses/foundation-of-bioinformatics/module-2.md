---
draft: true
---

# Module 2: Programming Fundamentals
**Python for Bioinformatics**

> **Core Concept**: Bioinformatics requires computational skills to manipulate, analyze, and extract insights from biological data. Python provides the tools—you provide the biological thinking.

---

## Module Overview

**Duration**: 3 weeks
**Prerequisites**: Module 1 (understanding biological data types)
**Programming Level**: Beginner-friendly (no prior experience required)

### Learning Objectives

By the end of this module, you will be able to:

1. ✓ Write Python programs to manipulate biological sequences
2. ✓ Use Biopython library for common bioinformatics tasks
3. ✓ Parse and analyze FASTA, FASTQ, and GenBank file formats
4. ✓ Access NCBI databases programmatically using Entrez
5. ✓ Write reusable, well-documented bioinformatics scripts
6. ✓ Apply best practices for reproducible computational biology

---

## Topics

### 1. Python Fundamentals

#### 🐍 Why Python for Bioinformatics?

**Advantages:**
- ✓ **Readable syntax** - Code looks like pseudocode
- ✓ **Rich ecosystem** - Biopython, NumPy, Pandas, Matplotlib
- ✓ **Interactive** - Test ideas quickly in Jupyter notebooks
- ✓ **Community** - Extensive bioinformatics resources

!!! info "Python vs. Other Languages"
    - **R**: Better for statistics/visualization
    - **Perl**: Legacy bioinformatics scripts (being replaced)
    - **Python**: Best balance for bioinformatics workflows

---

#### 📦 Data Types for Biological Data

=== "Strings (Sequences)"
    ```python
    # DNA sequence as string
    dna_seq = "ATGCGATCGTAGCTAGCT"

    # String operations
    length = len(dna_seq)  # 18
    first_codon = dna_seq[0:3]  # "ATG"
    gc_count = dna_seq.count('G') + dna_seq.count('C')  # 10

    # String methods
    rna_seq = dna_seq.replace('T', 'U')  # "AUGCGAUCGUAGCUAGCU"
    ```

    **Why strings?** Biological sequences are naturally text data

=== "Lists (Collections)"
    ```python
    # List of gene names
    genes = ["BRCA1", "TP53", "EGFR", "MYC"]

    # List operations
    genes.append("KRAS")  # Add element
    genes.sort()  # Sort alphabetically
    first_gene = genes[0]  # Access by index

    # List of expression values
    expression = [145.3, 523.8, 189.2, 856.1]
    mean_expr = sum(expression) / len(expression)
    ```

    **Why lists?** Store multiple values (gene names, counts, coordinates)

=== "Dictionaries (Mappings)"
    ```python
    # Map gene names to expression values
    gene_expression = {
        "BRCA1": 145.3,
        "TP53": 523.8,
        "EGFR": 189.2,
        "MYC": 856.1
    }

    # Dictionary operations
    brca1_expr = gene_expression["BRCA1"]  # 145.3
    gene_expression["KRAS"] = 234.5  # Add new entry

    # Check if gene exists
    if "TP53" in gene_expression:
        print(f"TP53 expression: {gene_expression['TP53']}")
    ```

    **Why dictionaries?** Natural for key-value relationships (gene→expression, codon→amino acid)

---

#### 🔁 Control Flow

**Making Decisions with Biological Data:**

```python
def classify_gc_content(sequence):
    """Classify sequence by GC content."""
    gc_count = sequence.count('G') + sequence.count('C')
    gc_percent = (gc_count / len(sequence)) * 100

    if gc_percent < 40:
        return "AT-rich"
    elif gc_percent < 60:
        return "Balanced"
    else:
        return "GC-rich"

# Example
seq = "ATGCGATCGTAGCTAGCT"
classification = classify_gc_content(seq)
print(f"Sequence is {classification}")  # "GC-rich"
```

**Iterating Over Biological Data:**

```python
# Process multiple sequences
sequences = ["ATGCGT", "GCGCGC", "ATATAT"]

for seq in sequences:
    gc = classify_gc_content(seq)
    print(f"{seq}: {gc}")

# Output:
# ATGCGT: Balanced
# GCGCGC: GC-rich
# ATATAT: AT-rich
```

---

#### 🔧 Functions: Reusable Bioinformatics Tools

```python
def reverse_complement(dna_seq):
    """
    Return the reverse complement of a DNA sequence.

    Args:
        dna_seq (str): DNA sequence (A, T, G, C)

    Returns:
        str: Reverse complement sequence

    Example:
        >>> reverse_complement("ATGC")
        'GCAT'
    """
    # Complement mapping
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    # Build complement
    comp_seq = ''.join([complement[base] for base in dna_seq])

    # Reverse
    return comp_seq[::-1]

# Test
original = "ATGCGATCG"
rev_comp = reverse_complement(original)
print(f"Original: {original}")
print(f"RevComp:  {rev_comp}")
```

!!! tip "Function Design Principles"
    1. **Single purpose** - One function, one task
    2. **Descriptive names** - `reverse_complement` not `rc`
    3. **Docstrings** - Explain what, why, and how
    4. **Type hints** - Document expected input/output types

---

### 2. Biopython: Bioinformatics Library

#### 📚 Introduction to Biopython

**Biopython** provides data structures and tools for:
- Sequence manipulation
- File format parsing (FASTA, GenBank, PDB)
- Database access (NCBI, UniProt)
- Sequence alignment
- Phylogenetics

**Installation:**
```bash
pip install biopython
```

---

#### 🧬 Seq Objects: Better Than Strings

=== "Basic Seq Usage"
    ```python
    from Bio.Seq import Seq

    # Create Seq object
    dna_seq = Seq("ATGCGATCGTAGCTAGCT")

    # Biological operations
    rna_seq = dna_seq.transcribe()
    print(rna_seq)  # AUGCGAUCGUAGCUAGCU

    # Reverse complement
    rev_comp = dna_seq.reverse_complement()
    print(rev_comp)  # AGCTAGCTACGATCGCAT

    # Translation
    protein = dna_seq.translate()
    print(protein)  # MRSSS*
    ```

=== "Why Seq vs String?"
    ```python
    # String limitations
    dna_string = "ATGCGT"
    # No biological methods
    # dna_string.transcribe()  # ❌ AttributeError

    # Seq advantages
    dna_seq = Seq("ATGCGT")
    rna_seq = dna_seq.transcribe()  # ✓ Works
    protein = dna_seq.translate()    # ✓ Works

    # Seq validates biological operations
    protein_seq = Seq("MKTAYIAK")
    # protein_seq.transcribe()  # ❌ Error: can't transcribe protein
    ```

=== "Seq Operations"
    ```python
    from Bio.Seq import Seq

    seq = Seq("ATGCGATCGTAGCT")

    # Count nucleotides
    print(f"A: {seq.count('A')}")  # 3
    print(f"G: {seq.count('G')}")  # 4

    # GC content
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    print(f"GC%: {gc_content * 100:.1f}")  # 57.1%

    # Find motifs
    position = seq.find("TCG")
    print(f"TCG found at position: {position}")  # 6

    # Slicing
    first_codon = seq[0:3]  # ATG
    second_codon = seq[3:6]  # CGA
    ```

---

### 3. Working with Biological File Formats

#### 📄 FASTA Format

**Structure:**
```
>seq_id description
ATGCGATCGTAGCTAGCTGATCGATCG
TCGATCGATCGTACGATCGATCGATCG
>another_seq more info
GCGCGCGCGCGCGCGCGC
```

**Parsing FASTA:**

```python
from Bio import SeqIO

# Read single sequence
for record in SeqIO.parse("sequence.fasta", "fasta"):
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Sequence: {record.seq}")
    print(f"Length: {len(record)}")
```

**Practical Example: Filter by Length**

```python
from Bio import SeqIO

def filter_by_length(input_file, output_file, min_length):
    """
    Filter sequences by minimum length.

    Args:
        input_file: Input FASTA file
        output_file: Output FASTA file
        min_length: Minimum sequence length
    """
    sequences = []

    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq) >= min_length:
            sequences.append(record)

    # Write filtered sequences
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Kept {len(sequences)} sequences >= {min_length} bp")

# Example usage
filter_by_length("all_seqs.fasta", "long_seqs.fasta", min_length=500)
```

---

#### 🧬 FASTQ Format (With Quality Scores)

**Structure:**
```
@seq_id
ATGCGATCGTAGCT
+
IIIHHGGGFFFEEE
```

**Quality Scores:**
- ASCII characters encode quality (Phred scores)
- `I` = high quality (Q=40, 99.99% accurate)
- `E` = lower quality (Q=36, 99.97% accurate)

**Parsing FASTQ:**

```python
from Bio import SeqIO

for record in SeqIO.parse("reads.fastq", "fastq"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")

    # Quality scores (Phred scores)
    qualities = record.letter_annotations["phred_quality"]
    mean_quality = sum(qualities) / len(qualities)
    print(f"Mean quality: {mean_quality:.1f}")
```

**Quality Filtering:**

```python
def filter_by_quality(input_fastq, output_fastq, min_quality=30):
    """Keep only high-quality reads."""
    high_quality = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        qualities = record.letter_annotations["phred_quality"]
        mean_qual = sum(qualities) / len(qualities)

        if mean_qual >= min_quality:
            high_quality.append(record)

    SeqIO.write(high_quality, output_fastq, "fastq")
    print(f"Kept {len(high_quality)} high-quality reads")
```

---

#### 🧬 GenBank Format (Rich Annotations)

**GenBank contains:**
- Sequence
- Features (genes, CDS, promoters)
- References
- Organism information

**Parsing GenBank:**

```python
from Bio import SeqIO

# Read GenBank file
record = SeqIO.read("NC_000913.gb", "genbank")

print(f"ID: {record.id}")
print(f"Description: {record.description}")
print(f"Organism: {record.annotations['organism']}")
print(f"Sequence length: {len(record.seq)}")

# Extract features
for feature in record.features:
    if feature.type == "CDS":  # Coding sequence
        gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
        location = feature.location
        print(f"Gene {gene_name} at {location}")
```

---

### 4. NCBI Programmatic Access

#### 🌐 Entrez: The NCBI API

**What is Entrez?**
- Unified API for all NCBI databases
- Programmatic access to GenBank, PubMed, SRA, etc.
- Free but requires email registration

**Setup:**

```python
from Bio import Entrez

# ALWAYS set your email (required by NCBI)
Entrez.email = "your.email@example.com"
```

!!! warning "NCBI Usage Policy"
    - **Provide your email** - Required by NCBI
    - **Limit requests** - Max 3 per second (10/sec with API key)
    - **Use Entrez.read()** - Parse XML responses
    - **Cache results** - Don't re-download unnecessarily

---

#### 🔍 Searching NCBI Databases

**Example: Search PubMed**

```python
from Bio import Entrez

Entrez.email = "your.email@example.com"

# Search PubMed
handle = Entrez.esearch(db="pubmed",
                        term="CRISPR AND 2023[PDAT]",
                        retmax=10)
record = Entrez.read(handle)
handle.close()

print(f"Found {record['Count']} articles")
print(f"First 10 PMIDs: {record['IdList']}")
```

**Example: Search Nucleotide Database**

```python
# Search for BRCA1 sequences
handle = Entrez.esearch(db="nucleotide",
                        term="BRCA1[Gene] AND Homo sapiens[Organism]",
                        retmax=5)
record = Entrez.read(handle)
handle.close()

print(f"Found {record['Count']} sequences")
for seq_id in record['IdList']:
    print(f"  {seq_id}")
```

---

#### 📥 Fetching Records from NCBI

**Fetch GenBank Record:**

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

# Fetch GenBank record
handle = Entrez.efetch(db="nucleotide",
                       id="NM_007294",  # BRCA1 mRNA
                       rettype="gb",
                       retmode="text")

# Parse as GenBank
record = SeqIO.read(handle, "genbank")
handle.close()

print(f"ID: {record.id}")
print(f"Description: {record.description}")
print(f"Length: {len(record.seq)} bp")
print(f"Organism: {record.annotations['organism']}")

# Extract CDS features
for feature in record.features:
    if feature.type == "CDS":
        print(f"Coding sequence: {feature.location}")
```

**Batch Download:**

```python
def download_sequences(id_list, output_file):
    """
    Download multiple sequences from NCBI.

    Args:
        id_list: List of GenBank IDs
        output_file: Output FASTA file
    """
    Entrez.email = "your.email@example.com"

    # Fetch all sequences at once (efficient)
    handle = Entrez.efetch(db="nucleotide",
                           id=id_list,
                           rettype="fasta",
                           retmode="text")

    # Parse and save
    records = SeqIO.parse(handle, "fasta")
    count = SeqIO.write(records, output_file, "fasta")
    handle.close()

    print(f"Downloaded {count} sequences")

# Example
ids = ["NM_007294", "NM_000546", "NM_001904"]
download_sequences(ids, "genes.fasta")
```

---

### 5. Best Practices for Bioinformatics Code

#### ✅ Code Organization

**Good Script Structure:**

```python
"""
analyze_sequences.py - Analyze GC content of FASTA sequences

Usage:
    python analyze_sequences.py input.fasta output.csv
"""

from Bio import SeqIO
import sys

def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def analyze_fasta(input_file, output_file):
    """
    Analyze GC content for all sequences in FASTA.

    Args:
        input_file: Input FASTA file path
        output_file: Output CSV file path
    """
    results = []

    for record in SeqIO.parse(input_file, "fasta"):
        gc = calculate_gc_content(str(record.seq))
        results.append(f"{record.id},{len(record.seq)},{gc:.2f}")

    # Write results
    with open(output_file, 'w') as f:
        f.write("seq_id,length,gc_content\n")
        for line in results:
            f.write(line + "\n")

    print(f"Analyzed {len(results)} sequences")

def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: python analyze_sequences.py input.fasta output.csv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    analyze_fasta(input_file, output_file)

if __name__ == "__main__":
    main()
```

---

#### 📝 Documentation

!!! tip "Documentation Levels"
    1. **Module docstring** - What the script does
    2. **Function docstrings** - What each function does
    3. **Inline comments** - Why (not what) for complex logic
    4. **README** - How to use the script

**Example Docstring:**

```python
def align_sequences(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Perform pairwise sequence alignment using Needleman-Wunsch.

    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        match (int): Score for matching residues (default: 1)
        mismatch (int): Penalty for mismatches (default: -1)
        gap (int): Penalty for gaps (default: -2)

    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)

    Example:
        >>> align_sequences("ATGC", "ATGC")
        ('ATGC', 'ATGC', 4)

    Raises:
        ValueError: If sequences are empty
    """
    # Implementation here...
```

---

#### 🔬 Testing Your Code

**Simple Testing:**

```python
def test_gc_content():
    """Test GC content calculation."""
    # Test case 1: All GC
    assert calculate_gc_content("GCGCGC") == 100.0

    # Test case 2: No GC
    assert calculate_gc_content("ATATAT") == 0.0

    # Test case 3: Mixed
    result = calculate_gc_content("ATGC")
    assert 49.0 < result < 51.0  # ~50%

    print("All tests passed!")

test_gc_content()
```

---

#### 🚀 Performance Considerations

??? note "Efficient vs. Inefficient Code"
    **Inefficient:**
    ```python
    # Reading entire file into memory
    sequences = []
    for record in SeqIO.parse("huge_file.fasta", "fasta"):
        sequences.append(record)  # ❌ Loads everything

    # Process sequences
    for seq in sequences:
        print(len(seq))
    ```

    **Efficient:**
    ```python
    # Process one at a time (streaming)
    for record in SeqIO.parse("huge_file.fasta", "fasta"):
        print(len(record))  # ✓ Process immediately, don't store
    ```

---

## 📝 Exercises

### Exercise 1: Sequence Statistics

Write a function that calculates comprehensive statistics for a DNA sequence:

```python
def sequence_stats(dna_seq):
    """
    Calculate statistics for a DNA sequence.

    Should return a dictionary with:
    - length
    - gc_content (%)
    - at_content (%)
    - nucleotide_counts (dict)
    """
    # Your code here
    pass

# Test
seq = "ATGCGATCGTAGCTAGCT"
stats = sequence_stats(seq)
print(stats)
# Expected output:
# {
#   'length': 18,
#   'gc_content': 55.56,
#   'at_content': 44.44,
#   'nucleotide_counts': {'A': 4, 'T': 4, 'G': 6, 'C': 4}
# }
```

### Exercise 2: FASTA Parser

Write a script that:
1. Reads a FASTA file
2. Calculates GC content for each sequence
3. Writes results to a CSV file with columns: `id,length,gc_content`

### Exercise 3: NCBI Download

Write a function that:
1. Searches NCBI for a gene name
2. Downloads the top 5 results
3. Saves them to a FASTA file

---

## 📚 Readings

### Required

1. **Python for Biologists** - Martin Jones (Chapters 1-5)
   *Focus*: Python fundamentals with biological examples

2. **Biopython Tutorial and Cookbook**
   *Focus*: SeqIO, Seq objects, Entrez

### Supplementary

3. **Python Documentation** - Built-in types and functions
4. **NCBI E-utilities Documentation** - Entrez API reference

---

## ✅ Self-Assessment

After completing this module, you should be able to:

- [ ] Write Python scripts using variables, loops, and functions
- [ ] Use Biopython's Seq object for sequence manipulation
- [ ] Parse FASTA, FASTQ, and GenBank files
- [ ] Access NCBI databases programmatically using Entrez
- [ ] Calculate basic sequence statistics (GC content, length, composition)
- [ ] Filter sequences based on criteria (length, quality, content)
- [ ] Write well-documented, reusable bioinformatics functions
- [ ] Test your code with simple assertions
- [ ] Follow best practices for code organization

!!! tip "Practice Makes Perfect"
    The only way to learn programming is by writing code. Complete all exercises and experiment with your own biological data.

---

## 🔗 Connection to Future Modules

!!! info "Why Programming Matters"
    **Module 3** (Databases) requires:
    - Entrez skills for querying NCBI
    - File parsing for local database files
    - Data structures for storing results

    **Module 4** (Sequence Alignment) requires:
    - Sequence manipulation
    - File I/O for reading/writing alignments
    - Functions for implementing algorithms

    **Module 5** (Statistics) requires:
    - Data structures (lists, dictionaries)
    - File parsing for count matrices
    - Python for statistical computing

**Next Module**: [Module 3: Biological Databases](module-3.md) - Now that you can program, you'll learn to access and integrate data from public repositories.

---

[↑ Course Index](index.md) | [← Module 1](module-1.md) | [Next: Module 3 →](module-3.md) | [🌐 View in Arabic](/ar/courses/foundation-of-bioinformatics/module-2/)
