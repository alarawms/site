---
draft: true
---

# الفصل الثاني: أساسيات البرمجة
**لغة Python للمعلوماتية الحيوية**

> **المفهوم الأساسي**: تتطلب المعلوماتية الحيوية مهارات حسابية للتعامل مع البيانات البيولوجية وتحليلها واستخراج الرؤى منها. Python توفر الأدوات—وأنت توفر التفكير البيولوجي.

---

## نظرة عامة على الفصل

**المدة**: 3 أسابيع
**المتطلبات الأساسية**: الفصل الأول (فهم أنواع البيانات البيولوجية)
**مستوى البرمجة**: مناسب للمبتدئين (لا يلزم خبرة سابقة)

### أهداف التعلم

بإكمال هذا الفصل، ستكون قادراً على:

1. ✓ كتابة برامج Python للتعامل مع التسلسلات البيولوجية
2. ✓ استخدام مكتبة Biopython للمهام الشائعة في المعلوماتية الحيوية
3. ✓ تحليل ملفات بتنسيقات FASTA و FASTQ و GenBank
4. ✓ الوصول إلى قواعد بيانات NCBI برمجياً باستخدام Entrez
5. ✓ كتابة نصوص برمجية قابلة لإعادة الاستخدام وموثقة جيداً للمعلوماتية الحيوية
6. ✓ تطبيق أفضل الممارسات للبيولوجيا الحسابية القابلة لإعادة الإنتاج

---

## المواضيع

### 1. أساسيات Python

#### 🐍 لماذا Python للمعلوماتية الحيوية؟

**المزايا:**
- ✓ **بناء جملة قابل للقراءة** - الكود يشبه الكود الوهمي
- ✓ **نظام بيئي غني** - Biopython و NumPy و Pandas و Matplotlib
- ✓ **تفاعلي** - اختبر الأفكار بسرعة في Jupyter notebooks
- ✓ **مجتمع** - موارد واسعة للمعلوماتية الحيوية

!!! info "Python مقابل اللغات الأخرى"
    - **R**: أفضل للإحصاء/التصور
    - **Perl**: نصوص المعلوماتية الحيوية القديمة (يتم استبدالها)
    - **Python**: أفضل توازن لسير عمل المعلوماتية الحيوية

---

#### 📦 أنواع البيانات للبيانات البيولوجية

=== "النصوص (التسلسلات)"
    ```python
    # تسلسل DNA كنص
    dna_seq = "ATGCGATCGTAGCTAGCT"

    # عمليات النصوص
    length = len(dna_seq)  # 18
    first_codon = dna_seq[0:3]  # "ATG"
    gc_count = dna_seq.count('G') + dna_seq.count('C')  # 10

    # طرق النصوص
    rna_seq = dna_seq.replace('T', 'U')  # "AUGCGAUCGUAGCUAGCU"
    ```

    **لماذا النصوص؟** التسلسلات البيولوجية هي بيانات نصية بطبيعتها

=== "القوائم (المجموعات)"
    ```python
    # قائمة بأسماء الجينات
    genes = ["BRCA1", "TP53", "EGFR", "MYC"]

    # عمليات القوائم
    genes.append("KRAS")  # إضافة عنصر
    genes.sort()  # ترتيب أبجدياً
    first_gene = genes[0]  # الوصول بالفهرس

    # قائمة بقيم التعبير
    expression = [145.3, 523.8, 189.2, 856.1]
    mean_expr = sum(expression) / len(expression)
    ```

    **لماذا القوائم؟** لتخزين قيم متعددة (أسماء جينات، عدادات، إحداثيات)

=== "القواميس (الربط)"
    ```python
    # ربط أسماء الجينات بقيم التعبير
    gene_expression = {
        "BRCA1": 145.3,
        "TP53": 523.8,
        "EGFR": 189.2,
        "MYC": 856.1
    }

    # عمليات القواميس
    brca1_expr = gene_expression["BRCA1"]  # 145.3
    gene_expression["KRAS"] = 234.5  # إضافة إدخال جديد

    # التحقق من وجود الجين
    if "TP53" in gene_expression:
        print(f"TP53 expression: {gene_expression['TP53']}")
    ```

    **لماذا القواميس؟** طبيعي للعلاقات مفتاح-قيمة (جين→تعبير، كودون→حمض أميني)

---

#### 🔁 التحكم في التدفق

**اتخاذ القرارات مع البيانات البيولوجية:**

```python
def classify_gc_content(sequence):
    """تصنيف التسلسل حسب محتوى GC."""
    gc_count = sequence.count('G') + sequence.count('C')
    gc_percent = (gc_count / len(sequence)) * 100

    if gc_percent < 40:
        return "AT-rich"
    elif gc_percent < 60:
        return "Balanced"
    else:
        return "GC-rich"

# مثال
seq = "ATGCGATCGTAGCTAGCT"
classification = classify_gc_content(seq)
print(f"Sequence is {classification}")  # "GC-rich"
```

**التكرار عبر البيانات البيولوجية:**

```python
# معالجة تسلسلات متعددة
sequences = ["ATGCGT", "GCGCGC", "ATATAT"]

for seq in sequences:
    gc = classify_gc_content(seq)
    print(f"{seq}: {gc}")

# الإخراج:
# ATGCGT: Balanced
# GCGCGC: GC-rich
# ATATAT: AT-rich
```

---

#### 🔧 الدوال: أدوات المعلوماتية الحيوية القابلة لإعادة الاستخدام

```python
def reverse_complement(dna_seq):
    """
    إرجاع المتمم العكسي لتسلسل DNA.

    Args:
        dna_seq (str): تسلسل DNA (A, T, G, C)

    Returns:
        str: تسلسل المتمم العكسي

    Example:
        >>> reverse_complement("ATGC")
        'GCAT'
    """
    # ربط المتممات
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    # بناء المتمم
    comp_seq = ''.join([complement[base] for base in dna_seq])

    # العكس
    return comp_seq[::-1]

# اختبار
original = "ATGCGATCG"
rev_comp = reverse_complement(original)
print(f"Original: {original}")
print(f"RevComp:  {rev_comp}")
```

!!! tip "مبادئ تصميم الدوال"
    1. **غرض واحد** - دالة واحدة، مهمة واحدة
    2. **أسماء وصفية** - `reverse_complement` وليس `rc`
    3. **وثائق** - اشرح ماذا ولماذا وكيف
    4. **تلميحات الأنواع** - وثق أنواع الإدخال/الإخراج المتوقعة

---

### 2. Biopython: مكتبة المعلوماتية الحيوية

#### 📚 مقدمة إلى Biopython

**Biopython** توفر هياكل بيانات وأدوات لـ:
- معالجة التسلسلات
- تحليل تنسيقات الملفات (FASTA، GenBank، PDB)
- الوصول إلى قواعد البيانات (NCBI، UniProt)
- محاذاة التسلسلات
- علم الوراثة العرقية

**التثبيت:**
```bash
pip install biopython
```

---

#### 🧬 كائنات Seq: أفضل من النصوص

=== "استخدام Seq الأساسي"
    ```python
    from Bio.Seq import Seq

    # إنشاء كائن Seq
    dna_seq = Seq("ATGCGATCGTAGCTAGCT")

    # عمليات بيولوجية
    rna_seq = dna_seq.transcribe()
    print(rna_seq)  # AUGCGAUCGUAGCUAGCU

    # المتمم العكسي
    rev_comp = dna_seq.reverse_complement()
    print(rev_comp)  # AGCTAGCTACGATCGCAT

    # الترجمة
    protein = dna_seq.translate()
    print(protein)  # MRSSS*
    ```

=== "لماذا Seq مقابل النص؟"
    ```python
    # قيود النصوص
    dna_string = "ATGCGT"
    # لا توجد طرق بيولوجية
    # dna_string.transcribe()  # ❌ AttributeError

    # مزايا Seq
    dna_seq = Seq("ATGCGT")
    rna_seq = dna_seq.transcribe()  # ✓ يعمل
    protein = dna_seq.translate()    # ✓ يعمل

    # Seq يتحقق من العمليات البيولوجية
    protein_seq = Seq("MKTAYIAK")
    # protein_seq.transcribe()  # ❌ خطأ: لا يمكن نسخ البروتين
    ```

=== "عمليات Seq"
    ```python
    from Bio.Seq import Seq

    seq = Seq("ATGCGATCGTAGCT")

    # عد النيوكليوتيدات
    print(f"A: {seq.count('A')}")  # 3
    print(f"G: {seq.count('G')}")  # 4

    # محتوى GC
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    print(f"GC%: {gc_content * 100:.1f}")  # 57.1%

    # البحث عن أنماط
    position = seq.find("TCG")
    print(f"TCG found at position: {position}")  # 6

    # التقطيع
    first_codon = seq[0:3]  # ATG
    second_codon = seq[3:6]  # CGA
    ```

---

### 3. العمل مع تنسيقات الملفات البيولوجية

#### 📄 تنسيق FASTA

**الهيكل:**
```
>seq_id description
ATGCGATCGTAGCTAGCTGATCGATCG
TCGATCGATCGTACGATCGATCGATCG
>another_seq more info
GCGCGCGCGCGCGCGCGC
```

**تحليل FASTA:**

```python
from Bio import SeqIO

# قراءة تسلسل واحد
for record in SeqIO.parse("sequence.fasta", "fasta"):
    print(f"ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Sequence: {record.seq}")
    print(f"Length: {len(record)}")
```

**مثال عملي: التصفية حسب الطول**

```python
from Bio import SeqIO

def filter_by_length(input_file, output_file, min_length):
    """
    تصفية التسلسلات حسب الحد الأدنى للطول.

    Args:
        input_file: ملف FASTA الإدخال
        output_file: ملف FASTA الإخراج
        min_length: الحد الأدنى لطول التسلسل
    """
    sequences = []

    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq) >= min_length:
            sequences.append(record)

    # كتابة التسلسلات المصفاة
    SeqIO.write(sequences, output_file, "fasta")
    print(f"Kept {len(sequences)} sequences >= {min_length} bp")

# مثال الاستخدام
filter_by_length("all_seqs.fasta", "long_seqs.fasta", min_length=500)
```

---

#### 🧬 تنسيق FASTQ (مع درجات الجودة)

**الهيكل:**
```
@seq_id
ATGCGATCGTAGCT
+
IIIHHGGGFFFEEE
```

**درجات الجودة:**
- أحرف ASCII تشفر الجودة (درجات Phred)
- `I` = جودة عالية (Q=40، 99.99% دقة)
- `E` = جودة أقل (Q=36، 99.97% دقة)

**تحليل FASTQ:**

```python
from Bio import SeqIO

for record in SeqIO.parse("reads.fastq", "fastq"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")

    # درجات الجودة (درجات Phred)
    qualities = record.letter_annotations["phred_quality"]
    mean_quality = sum(qualities) / len(qualities)
    print(f"Mean quality: {mean_quality:.1f}")
```

**تصفية الجودة:**

```python
def filter_by_quality(input_fastq, output_fastq, min_quality=30):
    """الاحتفاظ بالقراءات عالية الجودة فقط."""
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

#### 🧬 تنسيق GenBank (تعليقات توضيحية غنية)

**GenBank يحتوي على:**
- التسلسل
- الميزات (جينات، CDS، محفزات)
- المراجع
- معلومات الكائن الحي

**تحليل GenBank:**

```python
from Bio import SeqIO

# قراءة ملف GenBank
record = SeqIO.read("NC_000913.gb", "genbank")

print(f"ID: {record.id}")
print(f"Description: {record.description}")
print(f"Organism: {record.annotations['organism']}")
print(f"Sequence length: {len(record.seq)}")

# استخراج الميزات
for feature in record.features:
    if feature.type == "CDS":  # تسلسل الترميز
        gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
        location = feature.location
        print(f"Gene {gene_name} at {location}")
```

---

### 4. الوصول البرمجي إلى NCBI

#### 🌐 Entrez: واجهة برمجة NCBI

**ما هو Entrez؟**
- واجهة برمجة موحدة لجميع قواعد بيانات NCBI
- وصول برمجي إلى GenBank و PubMed و SRA وغيرها
- مجاني لكن يتطلب تسجيل البريد الإلكتروني

**الإعداد:**

```python
from Bio import Entrez

# دائماً اضبط بريدك الإلكتروني (مطلوب من NCBI)
Entrez.email = "your.email@example.com"
```

!!! warning "سياسة استخدام NCBI"
    - **قدم بريدك الإلكتروني** - مطلوب من NCBI
    - **حدد الطلبات** - حد أقصى 3 في الثانية (10/ثانية مع مفتاح API)
    - **استخدم Entrez.read()** - تحليل استجابات XML
    - **خزن النتائج** - لا تعيد التنزيل دون داعٍ

---

#### 🔍 البحث في قواعد بيانات NCBI

**مثال: البحث في PubMed**

```python
from Bio import Entrez

Entrez.email = "your.email@example.com"

# البحث في PubMed
handle = Entrez.esearch(db="pubmed",
                        term="CRISPR AND 2023[PDAT]",
                        retmax=10)
record = Entrez.read(handle)
handle.close()

print(f"Found {record['Count']} articles")
print(f"First 10 PMIDs: {record['IdList']}")
```

**مثال: البحث في قاعدة بيانات النيوكليوتيدات**

```python
# البحث عن تسلسلات BRCA1
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

#### 📥 جلب السجلات من NCBI

**جلب سجل GenBank:**

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

# جلب سجل GenBank
handle = Entrez.efetch(db="nucleotide",
                       id="NM_007294",  # BRCA1 mRNA
                       rettype="gb",
                       retmode="text")

# تحليل كـ GenBank
record = SeqIO.read(handle, "genbank")
handle.close()

print(f"ID: {record.id}")
print(f"Description: {record.description}")
print(f"Length: {len(record.seq)} bp")
print(f"Organism: {record.annotations['organism']}")

# استخراج ميزات CDS
for feature in record.features:
    if feature.type == "CDS":
        print(f"Coding sequence: {feature.location}")
```

**التنزيل الدفعي:**

```python
def download_sequences(id_list, output_file):
    """
    تنزيل تسلسلات متعددة من NCBI.

    Args:
        id_list: قائمة معرفات GenBank
        output_file: ملف FASTA الإخراج
    """
    Entrez.email = "your.email@example.com"

    # جلب جميع التسلسلات دفعة واحدة (فعال)
    handle = Entrez.efetch(db="nucleotide",
                           id=id_list,
                           rettype="fasta",
                           retmode="text")

    # التحليل والحفظ
    records = SeqIO.parse(handle, "fasta")
    count = SeqIO.write(records, output_file, "fasta")
    handle.close()

    print(f"Downloaded {count} sequences")

# مثال
ids = ["NM_007294", "NM_000546", "NM_001904"]
download_sequences(ids, "genes.fasta")
```

---

### 5. أفضل الممارسات لكود المعلوماتية الحيوية

#### ✅ تنظيم الكود

**هيكل نصي برمجي جيد:**

```python
"""
analyze_sequences.py - تحليل محتوى GC لتسلسلات FASTA

الاستخدام:
    python analyze_sequences.py input.fasta output.csv
"""

from Bio import SeqIO
import sys

def calculate_gc_content(sequence):
    """حساب نسبة محتوى GC."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def analyze_fasta(input_file, output_file):
    """
    تحليل محتوى GC لجميع التسلسلات في FASTA.

    Args:
        input_file: مسار ملف FASTA الإدخال
        output_file: مسار ملف CSV الإخراج
    """
    results = []

    for record in SeqIO.parse(input_file, "fasta"):
        gc = calculate_gc_content(str(record.seq))
        results.append(f"{record.id},{len(record.seq)},{gc:.2f}")

    # كتابة النتائج
    with open(output_file, 'w') as f:
        f.write("seq_id,length,gc_content\n")
        for line in results:
            f.write(line + "\n")

    print(f"Analyzed {len(results)} sequences")

def main():
    """نقطة الدخول الرئيسية."""
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

#### 📝 التوثيق

!!! tip "مستويات التوثيق"
    1. **وثائق الوحدة** - ماذا يفعل النص البرمجي
    2. **وثائق الدالة** - ماذا تفعل كل دالة
    3. **تعليقات سطرية** - لماذا (وليس ماذا) للمنطق المعقد
    4. **README** - كيفية استخدام النص البرمجي

**مثال التوثيق:**

```python
def align_sequences(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    إجراء محاذاة زوجية للتسلسل باستخدام Needleman-Wunsch.

    Args:
        seq1 (str): التسلسل الأول
        seq2 (str): التسلسل الثاني
        match (int): درجة البقايا المتطابقة (افتراضي: 1)
        mismatch (int): عقوبة عدم التطابق (افتراضي: -1)
        gap (int): عقوبة الفجوات (افتراضي: -2)

    Returns:
        tuple: (aligned_seq1, aligned_seq2, alignment_score)

    Example:
        >>> align_sequences("ATGC", "ATGC")
        ('ATGC', 'ATGC', 4)

    Raises:
        ValueError: إذا كانت التسلسلات فارغة
    """
    # التنفيذ هنا...
```

---

#### 🔬 اختبار الكود الخاص بك

**اختبار بسيط:**

```python
def test_gc_content():
    """اختبار حساب محتوى GC."""
    # حالة الاختبار 1: كل GC
    assert calculate_gc_content("GCGCGC") == 100.0

    # حالة الاختبار 2: لا GC
    assert calculate_gc_content("ATATAT") == 0.0

    # حالة الاختبار 3: مختلط
    result = calculate_gc_content("ATGC")
    assert 49.0 < result < 51.0  # ~50%

    print("All tests passed!")

test_gc_content()
```

---

#### 🚀 اعتبارات الأداء

??? note "كود فعال مقابل غير فعال"
    **غير فعال:**
    ```python
    # قراءة الملف بالكامل في الذاكرة
    sequences = []
    for record in SeqIO.parse("huge_file.fasta", "fasta"):
        sequences.append(record)  # ❌ يحمل كل شيء

    # معالجة التسلسلات
    for seq in sequences:
        print(len(seq))
    ```

    **فعال:**
    ```python
    # المعالجة واحداً تلو الآخر (دفق)
    for record in SeqIO.parse("huge_file.fasta", "fasta"):
        print(len(record))  # ✓ معالجة فورية، لا تخزين
    ```

---

## 📝 التمارين

### التمرين 1: إحصائيات التسلسل

اكتب دالة تحسب إحصائيات شاملة لتسلسل DNA:

```python
def sequence_stats(dna_seq):
    """
    حساب الإحصائيات لتسلسل DNA.

    يجب أن تُرجع قاموساً يحتوي على:
    - length
    - gc_content (%)
    - at_content (%)
    - nucleotide_counts (dict)
    """
    # كودك هنا
    pass

# اختبار
seq = "ATGCGATCGTAGCTAGCT"
stats = sequence_stats(seq)
print(stats)
# الإخراج المتوقع:
# {
#   'length': 18,
#   'gc_content': 55.56,
#   'at_content': 44.44,
#   'nucleotide_counts': {'A': 4, 'T': 4, 'G': 6, 'C': 4}
# }
```

### التمرين 2: محلل FASTA

اكتب نصاً برمجياً:
1. يقرأ ملف FASTA
2. يحسب محتوى GC لكل تسلسل
3. يكتب النتائج إلى ملف CSV بأعمدة: `id,length,gc_content`

### التمرين 3: تنزيل NCBI

اكتب دالة:
1. تبحث في NCBI عن اسم جين
2. تنزل أفضل 5 نتائج
3. تحفظها في ملف FASTA

---

## 📚 القراءات

### مطلوب

1. **Python for Biologists** - Martin Jones (الفصول 1-5)
   *التركيز*: أساسيات Python مع أمثلة بيولوجية

2. **Biopython Tutorial and Cookbook**
   *التركيز*: SeqIO، كائنات Seq، Entrez

### تكميلي

3. **Python Documentation** - الأنواع والوظائف المدمجة
4. **NCBI E-utilities Documentation** - مرجع Entrez API

---

## ✅ التقييم الذاتي

بعد إكمال هذا الفصل، يجب أن تكون قادراً على:

- [ ] كتابة نصوص Python باستخدام المتغيرات والحلقات والدوال
- [ ] استخدام كائن Seq من Biopython للتعامل مع التسلسلات
- [ ] تحليل ملفات FASTA و FASTQ و GenBank
- [ ] الوصول إلى قواعد بيانات NCBI برمجياً باستخدام Entrez
- [ ] حساب إحصائيات التسلسل الأساسية (محتوى GC، الطول، التركيب)
- [ ] تصفية التسلسلات بناءً على معايير (الطول، الجودة، المحتوى)
- [ ] كتابة دوال المعلوماتية الحيوية الموثقة جيداً وقابلة لإعادة الاستخدام
- [ ] اختبار الكود الخاص بك بتأكيدات بسيطة
- [ ] اتباع أفضل الممارسات لتنظيم الكود

!!! tip "الممارسة تصنع الإتقان"
    الطريقة الوحيدة لتعلم البرمجة هي بكتابة الكود. أكمل جميع التمارين وجرب بياناتك البيولوجية الخاصة.

---

## 🔗 الارتباط بالفصول المستقبلية

!!! info "لماذا البرمجة مهمة"
    **الفصل الثالث** (قواعد البيانات) يتطلب:
    - مهارات Entrez للاستعلام عن NCBI
    - تحليل الملفات لملفات قواعد البيانات المحلية
    - هياكل البيانات لتخزين النتائج

    **الفصل الرابع** (محاذاة التسلسلات) يتطلب:
    - معالجة التسلسلات
    - إدخال/إخراج الملفات لقراءة/كتابة المحاذاة
    - الدوال لتنفيذ الخوارزميات

    **الفصل الخامس** (الإحصاء) يتطلب:
    - هياكل البيانات (القوائم، القواميس)
    - تحليل الملفات لمصفوفات العد
    - Python للحوسبة الإحصائية

**الفصل التالي**: [الفصل الثالث: قواعد البيانات البيولوجية](module-3.md) - الآن بعد أن أصبحت قادراً على البرمجة، ستتعلم الوصول إلى البيانات ودمجها من المستودعات العامة.

---

[↑ فهرس المقرر](index.md) | [← الفصل الأول](module-1.md) | [التالي: الفصل الثالث →](module-3.md) | [🌐 عرض بالإنجليزية](/en/courses/foundation-of-bioinformatics/module-2/)
