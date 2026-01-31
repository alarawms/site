# Citation Quick Reference Card

Quick copy-paste templates for citing URLs and resources in your courses.

---

## 🚀 Quick Templates

### 1. Inline URL Citation

```markdown
According to the [NCBI documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/),
E-utilities provides programmatic access to databases.
```

---

### 2. Footnote Citation (Recommended for Academic)

```markdown
The Smith-Waterman algorithm[^1] provides optimal local alignments.

[^1]: Smith, T.F., & Waterman, M.S. (1981). Identification of common
      molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.
      DOI: [10.1016/0022-2836(81)90087-5](https://doi.org/10.1016/0022-2836(81)90087-5)
```

---

### 3. Reference Section at End

```markdown
## Your Content Here

Use reference [1] for background and [2] for methods.

---

## References

1. Author, A.B. (Year). Title. *Journal*, Vol(Issue), pages.
   [https://doi.org/xxxxx](https://doi.org/xxxxx)

2. Organization. (Year). *Resource Name*.
   Retrieved from [https://example.com](https://example.com)
```

---

### 4. Citation in Admonition Box

```markdown
!!! cite "Key Reference"
    **Watson, J.D., & Crick, F.H.C. (1953).** Molecular structure of
    nucleic acids. *Nature*, 171(4356), 737-738.

    DOI: [10.1038/171737a0](https://doi.org/10.1038/171737a0)
```

---

## 📚 Specific Templates

### Research Paper

```markdown
[^ref]: Author, A.B., Author, C.D., & Author, E.F. (Year). Title of paper.
        *Journal Name*, volume(issue), pages.
        DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)
```

**Example:**
```markdown
[^blast]: Altschul, S.F., et al. (1990). Basic local alignment search tool.
          *Journal of Molecular Biology*, 215(3), 403-410.
          DOI: [10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)
```

---

### Database/Website

```markdown
[^ref]: Organization. (Year). *Resource Name*.
        Retrieved Month Day, Year, from https://example.com
```

**Example:**
```markdown
[^ncbi]: NCBI. (2023). *GenBank Database*.
         Retrieved January 20, 2026, from https://www.ncbi.nlm.nih.gov/genbank/
```

---

### Software Tool

```markdown
[^ref]: Developer, A. (Year). *Software Name* (Version X.Y) [Software].
        Platform. https://github.com/user/repo
        DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx) (if available)
```

**Example:**
```markdown
[^diamond]: Buchfink, B., et al. (2015). *DIAMOND* (v2.0.15) [Software].
            GitHub. https://github.com/bbuchfink/diamond
            DOI: [10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)
```

---

### Book

```markdown
[^ref]: Author, A.B. (Year). *Book Title* (Edition). Publisher.
```

**Example:**
```markdown
[^mount]: Mount, D.W. (2004). *Bioinformatics: Sequence and Genome Analysis*
          (2nd ed.). Cold Spring Harbor Laboratory Press.
```

---

## 🎯 Bioinformatics-Specific Examples

### NCBI Database

```markdown
[^genbank]: NCBI Resource Coordinators. (2023). Database resources of the
            National Center for Biotechnology Information.
            *Nucleic Acids Research*, 51(D1), D1-D9.
            Website: https://www.ncbi.nlm.nih.gov/
            DOI: [10.1093/nar/gkac1032](https://doi.org/10.1093/nar/gkac1032)
```

---

### UniProt

```markdown
[^uniprot]: The UniProt Consortium. (2023). UniProt: the Universal Protein
            Knowledgebase in 2023. *Nucleic Acids Research*, 51(D1), D523-D531.
            Website: https://www.uniprot.org/
            DOI: [10.1093/nar/gkac1052](https://doi.org/10.1093/nar/gkac1052)
```

---

### Ensembl

```markdown
[^ensembl]: Martin, F.J., et al. (2023). Ensembl 2023. *Nucleic Acids Research*,
            51(D1), D933-D941.
            Website: https://www.ensembl.org/
            DOI: [10.1093/nar/gkac958](https://doi.org/10.1093/nar/gkac958)
```

---

### BLAST

```markdown
[^blast]: Altschul, S.F., et al. (1990). Basic local alignment search tool.
          *Journal of Molecular Biology*, 215(3), 403-410.
          DOI: [10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)
          Web interface: https://blast.ncbi.nlm.nih.gov/
```

---

### Bioconductor Package

```markdown
[^deseq2]: Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of
           fold change and dispersion for RNA-seq data with DESeq2.
           *Genome Biology*, 15(12), 550.
           Package: https://bioconductor.org/packages/DESeq2/
           DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
```

---

### Galaxy Workflow

```markdown
[^galaxy]: Galaxy Community. (2023). *RNA-Seq Analysis Tutorial*.
           Galaxy Training Network.
           https://training.galaxyproject.org/training-material/topics/transcriptomics/
           Accessed: January 20, 2026
```

---

## 💡 Best Practices

### ✅ DO

- Use footnotes `[^1]` for academic rigor
- Include DOI links when available
- Add access dates for web resources
- Group references at module end
- Use consistent citation style

### ❌ DON'T

- Don't use "click here" as link text
- Don't forget access dates for web resources
- Don't mix citation styles in one course
- Don't cite unreliable sources

---

## 📋 Complete Module Template

```markdown
# Module Title

Your introduction with citation[^intro].

## Section

Content referencing database[^database] and tool[^tool].

### Subsection

More content with inline links to [NCBI](https://www.ncbi.nlm.nih.gov/).

!!! cite "Important Reference"
    **Key Author.** (Year). Important paper. *Journal*.
    DOI: [link](https://doi.org/xxxxx)

---

## References

[^intro]: Author, A. (Year). Background paper. *Journal*, vol, pages.
          DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)

[^database]: Organization. (Year). *Database Name*.
             https://example.com
             Accessed: Month Day, Year

[^tool]: Developer, A. (Year). *Software* (vX.Y).
         https://github.com/user/repo

---

[← Previous](prev.md) | [↑ Index](index.md) | [Next →](next.md)
```

---

## 🔧 Emacs Integration

In your Emacs org-mode files, you can maintain citations in BibTeX format:

```bibtex
@article{smith1981,
  author = {Smith, T.F. and Waterman, M.S.},
  title = {Identification of common molecular subsequences},
  journal = {Journal of Molecular Biology},
  year = {1981},
  volume = {147},
  number = {1},
  pages = {195--197},
  doi = {10.1016/0022-2836(81)90087-5}
}
```

Then use org-cite to insert citations, export to markdown, and copy to MkDocs.

---

**See Also:**
- [Citation Guide](citation-guide.md) - Full documentation
- [Citation Example](citation-example.md) - BLAST module with citations
- [Emacs Workflow](emacs-workflow-cheatsheet.md) - Integration tips
