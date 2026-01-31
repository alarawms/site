# MkDocs + Emacs Workflow Cheatsheet

Quick reference for creating rich content courses using your Emacs MkDocs workflow.

---

## 🚀 Quick Start

### Create a New Course (from Emacs)

1. Press `SPC w C` (your keybinding for `mkdocs-new-course`)
2. Enter course title: "Foundation in Bioinformatics"
3. Enter number of modules: 8
4. Files auto-generated with navigation!

**Result:**
```
docs/courses/foundation-in-bioinformatics/
├── index.md
├── module-1.md
├── module-2.md
├── module-3.md
├── module-4.md
├── module-5.md
├── module-6.md
├── module-7.md
└── module-8.md
```

---

## 📝 Insert Rich Content (Emacs Commands)

### Your Custom Keybindings

| Command | Keybinding | Description |
|---------|------------|-------------|
| `mkdocs-insert-admonition` | `SPC w i a` | Insert info box |
| `mkdocs-insert-code-block` | `SPC w i c` | Insert code block |
| `mkdocs-insert-mermaid` | `SPC w i m` | Insert diagram |
| `mkdocs-insert-table` | `SPC w i t` | Insert table |

### Quick Snippets (Type Manually)

#### Image with Lightbox
```markdown
![DNA Structure](../images/dna-structure.png)
```

#### External Link
```markdown
Visit [NCBI](https://www.ncbi.nlm.nih.gov/) database.
```

#### Internal Link
```markdown
See [Module 2](module-2.md) for alignment algorithms.
```

#### Admonition
Use `SPC w i a` or type:
```markdown
!!! tip
    Pro tip: Always validate your sequence data!
```

#### Code Block
Use `SPC w i c` or type:
```markdown
```python
def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100
` ``
```

#### Mermaid Diagram
Use `SPC w i m` or type:
```markdown
` ``mermaid
graph LR
    A[Input] --> B[Process]
    B --> C[Output]
` ``
```

---

## 🎯 Common Patterns for Bioinformatics Courses

### Pattern 1: Workflow Diagram + Code

```markdown
## Sequence Alignment Workflow

` ``mermaid
graph LR
    A[FASTQ Files] --> B[Quality Control]
    B --> C[Alignment]
    C --> D[BAM File]
` ``

### Implementation

` ``python
import pysam

# Load alignment file
bamfile = pysam.AlignmentFile("sample.bam", "rb")
for read in bamfile.fetch():
    print(read.query_name, read.reference_start)
` ``
```

### Pattern 2: Comparison Tabs

```markdown
## Alignment Tools Comparison

=== "BWA"
    ` ``bash
    bwa mem ref.fa reads.fq > aligned.sam
    ` ``

=== "Bowtie2"
    ` ``bash
    bowtie2 -x ref -U reads.fq -S aligned.sam
    ` ``

=== "STAR"
    ` ``bash
    STAR --genomeDir ref_dir --readFilesIn reads.fq
    ` ``
```

### Pattern 3: Interactive Data Visualization

```markdown
## Quality Score Distribution

` ``vegalite
{
  "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
  "data": {
    "values": [
      {"position": 1, "q30": 95},
      {"position": 50, "q30": 92},
      {"position": 100, "q30": 88}
    ]
  },
  "mark": "line",
  "encoding": {
    "x": {"field": "position", "type": "quantitative"},
    "y": {"field": "q30", "type": "quantitative"}
  }
}
` ``
```

### Pattern 4: Reference Tables from CSV

Create `docs/data/enzymes.csv`:
```csv
Enzyme,Function,Recognition Site
EcoRI,Restriction enzyme,GAATTC
BamHI,Restriction enzyme,GGATCC
HindIII,Restriction enzyme,AAGCTT
```

Then in your module (remove spaces between braces):
```markdown
## Common Restriction Enzymes

{ { read_csv('data/enzymes.csv') } }
```

**Live example:**

## Common Restriction Enzymes

{{ read_csv('data/enzymes.csv') }}

---

## 📁 Recommended Directory Structure

```
docs/
├── courses/
│   └── foundation-in-bioinformatics/
│       ├── index.md
│       ├── module-1.md
│       ├── module-2.md
│       ├── images/           # Course-specific images
│       │   ├── dna-structure.png
│       │   └── alignment-workflow.png
│       └── data/             # Course datasets
│           ├── quality-metrics.csv
│           └── gene-expression.csv
├── images/                    # Global images
└── data/                      # Global datasets
```

---

## 🔄 Complete Workflow Example

### 1. Create Course (Emacs)
```
SPC w C
→ "Foundation in Bioinformatics"
→ 8 modules
```

### 2. Edit Module (Emacs)
```
SPC w o          # Open MkDocs root
Navigate to module-1.md
```

### 3. Add Content
```markdown
# Module 1: Introduction to Bioinformatics

## What is Bioinformatics?

!!! info "Definition"
    Bioinformatics combines biology, computer science, and statistics
    to analyze biological data.

## The Central Dogma

` ``mermaid
sequenceDiagram
    DNA->>RNA: Transcription
    RNA->>Protein: Translation
` ``

## Key Databases

- [NCBI](https://www.ncbi.nlm.nih.gov/)
- [UniProt](https://www.uniprot.org/)
- [Ensembl](https://www.ensembl.org/)
```

### 4. Preview (Emacs)
```
SPC w p          # Start preview server
```

### 5. Publish (Emacs)
```
SPC w P          # Publish with custom message
# or
SPC w Q          # Quick publish with timestamp
```

---

## 💡 Pro Tips

### 1. Image Workflow
```bash
# Store images in course directory
docs/courses/foundation-in-bioinformatics/images/

# Reference in markdown
![Diagram](images/workflow.png)
```

### 2. Reusable Code Snippets
Create `docs/snippets/common-imports.py`:
```python
import pandas as pd
import numpy as np
from Bio import SeqIO
```

Include in modules:
```markdown
` ``python
--8<-- "snippets/common-imports.py"

# Your code here
` ``
```

### 3. Cross-Reference Modules
```markdown
See [alignment algorithms](module-3.md#smith-waterman) in Module 3.

Review [quality control](module-2.md) before proceeding.
```

### 4. Math Equations
```markdown
The alignment score is calculated as:

$$
S = \sum_{i=1}^{n} s(a_i, b_i) - \sum_{j=1}^{m} g_j
$$

where $s$ is the substitution score and $g$ is the gap penalty.
```

---

## 🎨 Style Conventions

### Consistent Heading Structure
```markdown
# Module Title (H1)

## Main Section (H2)

### Subsection (H3)

#### Detail (H4)
```

### Admonition Usage
- `!!! note` - General information
- `!!! tip` - Best practices
- `!!! warning` - Common pitfalls
- `!!! example` - Code examples
- `!!! question` - Exercises

### Code Language Tags
Always specify language:
- `python` - Python code
- `bash` - Shell commands
- `r` - R code
- `text` - Plain text output
- `mermaid` - Diagrams
- `vegalite` - Charts

---

## 📚 Quick Reference Links

- [Rich Content Guide](rich-content-guide.md) - Complete reference
- [Example Module](example-module.md) - Working example
- [Plugin Demo](plugin-demo.md) - All plugins demonstrated

---

## ⌨️ Your Complete MkDocs Keybinding Map

```elisp
SPC w n     # New blog post
SPC w C     # New course
SPC w m     # New module (add to existing course)
SPC w p     # Preview site (start server)
SPC w s     # Stop preview
SPC w b     # Build site
SPC w P     # Publish (git commit + push)
SPC w Q     # Quick publish (auto message)
SPC w o     # Open project root
SPC w c     # Open mkdocs.yml config

# Insert helpers
SPC w i a   # Insert admonition
SPC w i c   # Insert code block
SPC w i m   # Insert mermaid diagram
SPC w i t   # Insert table
```

---

## 🚀 Next Steps

1. ✅ Create your first course: `SPC w C`
2. ✅ Check [example-module.md](example-module.md) for inspiration
3. ✅ Use [rich-content-guide.md](rich-content-guide.md) as reference
4. ✅ Preview frequently: `SPC w p`
5. ✅ Publish when ready: `SPC w P`

Happy course creation! 🎓
