# Citation and Reference Guide

Complete guide for citing URLs, papers, and creating bibliographies in your MkDocs courses.

---

## 🔗 Method 1: Inline URL Citations

### Basic Inline Citation

```markdown
According to [NCBI documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/),
the E-utilities API provides programmatic access to databases.
```

**Renders as:**

According to [NCBI documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/),
the E-utilities API provides programmatic access to databases.

### Inline with Access Date

```markdown
The UniProt database [1] contains over 200 million protein sequences
(accessed January 2026).

[1]: https://www.uniprot.org/
```

---

## 📚 Method 2: Reference-Style Links (Recommended)

### At End of Section

```markdown
## Sequence Alignment

The Smith-Waterman algorithm [1] is the standard for local alignment,
while BLAST [2] provides a faster heuristic approach. Modern tools like
Diamond [3] offer even greater speed.

### References

[1]: Smith, T.F., & Waterman, M.S. (1981). [Identification of common
molecular subsequences](https://doi.org/10.1016/0022-2836(81)90087-5).
*Journal of Molecular Biology*, 147(1), 195-197.

[2]: Altschul, S.F., et al. (1990). [Basic local alignment search tool]
(https://www.ncbi.nlm.nih.gov/pubmed/2231712). *Journal of Molecular Biology*,
215(3), 403-410.

[3]: Buchfink, B., et al. (2015). [Fast and sensitive protein alignment
using DIAMOND](https://www.nature.com/articles/nmeth.3176). *Nature Methods*,
12(1), 59-60.
```

### Numbered References

```markdown
## Next-Generation Sequencing

Illumina sequencing uses sequencing-by-synthesis chemistry [1], while
Oxford Nanopore uses nanopore technology for long reads [2]. PacBio offers
HiFi reads combining length and accuracy [3].

---

### References

1. Bentley, D.R., et al. (2008). Accurate whole human genome sequencing
   using reversible terminator chemistry. *Nature*, 456(7218), 53-59.
   [https://doi.org/10.1038/nature07517](https://doi.org/10.1038/nature07517)

2. Oxford Nanopore Technologies. (2023). *Nanopore Sequencing Technology*.
   Retrieved from [https://nanoporetech.com/how-it-works](https://nanoporetech.com/how-it-works)

3. Pacific Biosciences. (2023). *HiFi Sequencing*.
   [https://www.pacb.com/technology/hifi-sequencing/](https://www.pacb.com/technology/hifi-sequencing/)
```

---

## 📖 Method 3: Footnote Citations

Perfect for academic-style references.

```markdown
The ENCODE project[^1] has identified functional elements across the human genome.
RNA-seq analysis[^2] has become the standard for transcriptome profiling.

[^1]: ENCODE Project Consortium. (2012). An integrated encyclopedia of DNA
      elements in the human genome. *Nature*, 489(7414), 57-74.
      DOI: [10.1038/nature11247](https://doi.org/10.1038/nature11247)

[^2]: Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: a revolutionary
      tool for transcriptomics. *Nature Reviews Genetics*, 10(1), 57-63.
      Available at: [https://www.nature.com/articles/nrg2484](https://www.nature.com/articles/nrg2484)
```

**Renders as:**

The ENCODE project[^1] has identified functional elements across the human genome.
RNA-seq analysis[^2] has become the standard for transcriptome profiling.

[^1]: ENCODE Project Consortium. (2012). An integrated encyclopedia of DNA
      elements in the human genome. *Nature*, 489(7414), 57-74.
      DOI: [10.1038/nature11247](https://doi.org/10.1038/nature11247)

[^2]: Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: a revolutionary
      tool for transcriptomics. *Nature Reviews Genetics*, 10(1), 57-63.
      Available at: [https://www.nature.com/articles/nrg2484](https://www.nature.com/articles/nrg2484)

---

## 🎯 Method 4: Admonition References

For important or highlighted citations.

```markdown
!!! cite "Key Reference"
    **Watson, J.D., & Crick, F.H.C. (1953).** Molecular structure of
    nucleic acids: A structure for deoxyribose nucleic acid. *Nature*,
    171(4356), 737-738.

    DOI: [10.1038/171737a0](https://doi.org/10.1038/171737a0)

    This seminal paper proposed the double helix structure of DNA.

!!! info "Online Resource"
    **NCBI Resource Coordinators (2023).** Database resources of the
    National Center for Biotechnology Information. *Nucleic Acids Research*.

    Website: [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)

    Accessed: January 20, 2026
```

**Renders as:**

!!! cite "Key Reference"
    **Watson, J.D., & Crick, F.H.C. (1953).** Molecular structure of
    nucleic acids: A structure for deoxyribose nucleic acid. *Nature*,
    171(4356), 737-738.

    DOI: [10.1038/171737a0](https://doi.org/10.1038/171737a0)

    This seminal paper proposed the double helix structure of DNA.

---

## 📊 Method 5: Reference Table

For courses with many citations.

```markdown
## Bibliography

| # | Citation | Link | Notes |
|---|----------|------|-------|
| 1 | Smith & Waterman (1981) | [DOI](https://doi.org/10.1016/0022-2836(81)90087-5) | Smith-Waterman algorithm |
| 2 | NCBI E-utilities | [Docs](https://www.ncbi.nlm.nih.gov/books/NBK25499/) | API documentation |
| 3 | UniProt Database | [Web](https://www.uniprot.org/) | Protein sequences |
| 4 | Ensembl Genome Browser | [Web](https://www.ensembl.org/) | Genome annotations |
```

**Renders as:**

| # | Citation | Link | Notes |
|---|----------|------|-------|
| 1 | Smith & Waterman (1981) | [DOI](https://doi.org/10.1016/0022-2836(81)90087-5) | Smith-Waterman algorithm |
| 2 | NCBI E-utilities | [Docs](https://www.ncbi.nlm.nih.gov/books/NBK25499/) | API documentation |
| 3 | UniProt Database | [Web](https://www.uniprot.org/) | Protein sequences |
| 4 | Ensembl Genome Browser | [Web](https://www.ensembl.org/) | Genome annotations |

---

## 🌐 Citation Formats for Different Web Resources

### Database Citation

```markdown
**NCBI Resource Coordinators.** (2023). Database resources of the National
Center for Biotechnology Information. *Nucleic Acids Research*, 51(D1), D1-D9.

- DOI: [10.1093/nar/gkac1032](https://doi.org/10.1093/nar/gkac1032)
- Website: [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)
- Accessed: January 20, 2026
```

### Software/Tool Citation

```markdown
**Langmead, B., & Salzberg, S.L.** (2012). Fast gapped-read alignment
with Bowtie 2. *Nature Methods*, 9(4), 357-359.

- Software: Bowtie2
- Version: 2.5.0
- URL: [http://bowtie-bio.sourceforge.net/bowtie2/](http://bowtie-bio.sourceforge.net/bowtie2/)
- GitHub: [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)
```

### Documentation/Tutorial Citation

```markdown
**Galaxy Community.** (2023). *Galaxy Training Materials: RNA-Seq Analysis*.
Retrieved January 20, 2026, from
[https://training.galaxyproject.org/training-material/topics/transcriptomics/](https://training.galaxyproject.org/training-material/topics/transcriptomics/)
```

### Preprint Citation

```markdown
**Author, A.B., & Author, C.D.** (2023). Title of preprint.
*bioRxiv* preprint.

DOI: [10.1101/2023.01.01.123456](https://doi.org/10.1101/2023.01.01.123456)

!!! warning
    Preprint - not yet peer reviewed
```

### GitHub Repository

```markdown
**Developer, A.** (2023). *Repository Name* (Version 1.2.3) [Software].
GitHub. [https://github.com/user/repo](https://github.com/user/repo)
```

---

## 🎓 Complete Module Example

Here's how to structure a complete module with references:

```markdown
# Module 3: Sequence Alignment

## Introduction

Sequence alignment is fundamental to bioinformatics[^align]. The Smith-Waterman
algorithm[^sw] provides optimal local alignment, while BLAST[^blast] offers
a practical heuristic approach.

## Algorithms

### Smith-Waterman

The Smith-Waterman algorithm uses dynamic programming to find optimal local
alignments[^sw]. Implementation details can be found in the original paper
and modern tutorials[^biostar].

!!! cite "Original Publication"
    **Smith, T.F., & Waterman, M.S.** (1981). Identification of common
    molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.

    DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5)

### BLAST

BLAST revolutionized sequence similarity searching by introducing heuristic
methods[^blast]. The NCBI BLAST suite[^ncbi-blast] remains the most widely
used implementation.

## Practical Tools

Common alignment tools include:

| Tool | Purpose | Citation | URL |
|------|---------|----------|-----|
| BLAST | Database searching | Altschul et al. (1990) | [NCBI](https://blast.ncbi.nlm.nih.gov/) |
| Bowtie2 | Short read alignment | Langmead & Salzberg (2012) | [Website](http://bowtie-bio.sourceforge.net/bowtie2/) |
| BWA | Read mapping | Li & Durbin (2009) | [GitHub](https://github.com/lh3/bwa) |

## Further Reading

- [NCBI BLAST Documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [Biostars: Sequence Alignment Tutorials](https://www.biostars.org/t/sequence-alignment/)

---

## References

[^align]: Mount, D.W. (2004). *Bioinformatics: Sequence and Genome Analysis*
          (2nd ed.). Cold Spring Harbor Laboratory Press.

[^sw]: Smith, T.F., & Waterman, M.S. (1981). Identification of common
       molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.
       [https://doi.org/10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5)

[^blast]: Altschul, S.F., Gish, W., Miller, W., Myers, E.W., & Lipman, D.J. (1990).
          Basic local alignment search tool. *Journal of Molecular Biology*, 215(3), 403-410.
          [https://doi.org/10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)

[^ncbi-blast]: NCBI. (2023). *BLAST: Basic Local Alignment Search Tool*.
               National Center for Biotechnology Information.
               [https://blast.ncbi.nlm.nih.gov/](https://blast.ncbi.nlm.nih.gov/)
               Accessed January 20, 2026.

[^biostar]: Biostars Community. (2023). *Sequence Alignment Tutorials*.
            [https://www.biostars.org/](https://www.biostars.org/)

---

**Navigation:** [← Module 2](module-2.md) | [↑ Index](index.md) | [Module 4 →](module-4.md)
```

---

## 💡 Best Practices

### 1. **Consistent Format**

Choose one citation style and stick to it:
- **APA**: Author (Year). Title. *Journal*, volume(issue), pages.
- **IEEE**: [1] Author, "Title," *Journal*, vol. X, no. Y, pp. ZZ-ZZ, Year.
- **Vancouver**: Author AB. Title. Journal. Year;Vol(Issue):Pages.

### 2. **Include Access Dates for Web Resources**

```markdown
Retrieved January 20, 2026, from https://example.com
```

### 3. **Use DOIs When Available**

```markdown
DOI: [10.1038/nature12345](https://doi.org/10.1038/nature12345)
```

### 4. **Distinguish Resource Types**

```markdown
!!! info "Database"
    UniProt - Protein sequences and annotations

!!! example "Software"
    STAR - RNA-seq aligner

!!! quote "Paper"
    Smith & Waterman (1981) - Original algorithm
```

### 5. **Group by Module**

Create a references section at the end of each module rather than one massive bibliography.

---

## 📝 Quick Templates

### Template 1: Research Paper

```markdown
[^ref]: Author, A.B., Author, C.D., & Author, E.F. (Year). Title of paper.
        *Journal Name*, volume(issue), pages.
        DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)
```

### Template 2: Website/Database

```markdown
[^ref]: Organization/Author. (Year). *Resource Name*.
        Retrieved Month Day, Year, from https://example.com
```

### Template 3: Software

```markdown
[^ref]: Developer, A. (Year). *Software Name* (Version X.Y.Z) [Software].
        Platform. https://example.com
        DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx) (if available)
```

### Template 4: Book

```markdown
[^ref]: Author, A.B. (Year). *Book Title* (Edition). Publisher.
        ISBN: XXXX-XXXX-XXXX
```

---

## 🔧 Integration with Your Emacs Workflow

Your Emacs config has **citar** and **org-roam** for bibliography management:

### From Org-Mode to MkDocs

1. **Maintain bibliography** in `~/org/bibliography/references/library.bib`
2. **Export formatted citations** to your markdown modules
3. Use **org-cite** syntax in Org files, then export

Example Org source:
```org
* Module Notes
According to [cite:@smith1981] the algorithm uses dynamic programming.

#+print_bibliography:
```

Export to Markdown, then copy to your MkDocs course.

---

## 📚 Common Bioinformatics Resources to Cite

### Essential Databases

```markdown
- **NCBI GenBank**: https://www.ncbi.nlm.nih.gov/genbank/
- **UniProt**: https://www.uniprot.org/
- **Ensembl**: https://www.ensembl.org/
- **PDB (Protein Data Bank)**: https://www.rcsb.org/
- **GEO (Gene Expression Omnibus)**: https://www.ncbi.nlm.nih.gov/geo/
```

### Essential Tools

```markdown
- **BLAST**: https://blast.ncbi.nlm.nih.gov/
- **Bioconductor**: https://www.bioconductor.org/
- **Galaxy**: https://usegalaxy.org/
- **UCSC Genome Browser**: https://genome.ucsc.edu/
```

### Essential Papers (Foundational)

```markdown
- Watson & Crick (1953) - DNA structure
- Smith & Waterman (1981) - Sequence alignment
- Altschul et al. (1990) - BLAST algorithm
- Lander et al. (2001) - Human genome sequence
```

---

## 🎯 Summary

**For Academic Rigor:** Use footnotes (`[^1]`) with full citations at module end

**For Quick References:** Use inline links with descriptive text

**For Important Sources:** Use admonition blocks with `!!! cite`

**For Many Citations:** Use reference tables or numbered lists

**Best Practice:** Combine methods - inline for casual mentions, footnotes for key references

---

**See Also:**
- [Rich Content Guide](rich-content-guide.md) - For other content types
- [Example Module](example-module.md) - See citations in context
