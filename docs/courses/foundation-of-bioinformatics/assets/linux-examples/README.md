# Linux Training Example Files
**للفصل الثاني: أساسيات Linux**

---

## 📦 الملفات المتوفرة

### 1. `sequences.fasta`
**ملف FASTA بسيط - 3 تسلسلات**

تسلسلات DNA قصيرة لثلاثة جينات:
- gene1 (BRCA1)
- gene2 (TP53)
- gene3 (EGFR)

**الاستخدام:**
```bash
# عد عدد التسلسلات
grep -c ">" sequences.fasta

# عرض أسماء الجينات
grep ">" sequences.fasta

# استخراج الأسماء فقط
grep ">" sequences.fasta | sed 's/>//'
```

---

### 2. `genome_reads.fastq`
**ملف FASTQ - 5 قراءات تسلسل**

قراءات تسلسل DNA مع درجات الجودة (Phred scores).

**الاستخدام:**
```bash
# عرض أول 10 أسطر
head genome_reads.fastq

# عد القراءات (كل قراءة = 4 أسطر)
wc -l genome_reads.fastq  # ثم قسم على 4

# البحث عن قراءة معينة
grep "SRR001666.3" genome_reads.fastq
```

---

### 3. `genome.fasta`
**ملف FASTA كبير - 5 تسلسلات جينومية**

تسلسلات أطول من الكروموسومات البشرية:
- BRCA1 region (chromosome 17)
- TP53 region (chromosome 17)
- EGFR region (chromosome 7)
- Mitochondrial DNA
- Ribosomal RNA gene

**الاستخدام:**
```bash
# عد التسلسلات
grep -c ">" genome.fasta

# حساب الحجم
wc -l genome.fasta
wc -c genome.fasta

# البحث عن كروموسوم معين
grep "chromosome 17" genome.fasta
```

---

### 4. `genes.txt`
**قائمة جينات - 10 جينات سرطانية**

جدول مفصول بـ TAB:
- اسم الجين
- الاسم الكامل
- الوظيفة

**الاستخدام:**
```bash
# البحث عن جين BRCA1
grep "BRCA1" genes.txt

# البحث بتجاهل حالة الأحرف
grep -i "cancer" genes.txt

# عد الجينات
wc -l genes.txt
```

---

### 5. `annotations.txt`
**ملف تعليقات جينية - 10 عناصر**

مواقع العناصر الجينية على الكروموسومات:
- Promoters
- Genes
- Exons
- Introns

**الاستخدام:**
```bash
# البحث عن exons فقط
grep "exon" annotations.txt

# عرض مع أرقام الأسطر
grep -n "gene" annotations.txt

# عد العناصر حسب النوع
grep -c "promoter" annotations.txt
```

---

## 💾 كيفية التحميل

### الطريقة 1: تحميل الكل (مضغوط)
```bash
# تحميل كملف ZIP
wget https://malarawi.sa/courses/foundation-of-bioinformatics/assets/linux-examples.zip

# فك الضغط
unzip linux-examples.zip

# الدخول للمجلد
cd linux-examples
```

### الطريقة 2: تحميل ملف واحد
```bash
# مثال: تحميل sequences.fasta
wget https://malarawi.sa/courses/foundation-of-bioinformatics/assets/linux-examples/sequences.fasta
```

### الطريقة 3: من GitHub مباشرة
```bash
# استنساخ المستودع
git clone https://github.com/alarawms/site.git
cd site/docs/courses/foundation-of-bioinformatics/assets/linux-examples/
```

---

## 🎯 التمارين المقترحة

### تمرين 1: الأساسيات
```bash
# عد كل شيء
wc -l *.fasta
wc -l *.fastq
wc -l *.txt

# عرض محتويات
head sequences.fasta
tail genes.txt
```

### تمرين 2: البحث والتصفية
```bash
# ابحث عن BRCA
grep "BRCA" genes.txt
grep "BRCA" genome.fasta

# عد تكرارات
grep -c "gene" annotations.txt
```

### تمرين 3: الأنابيب (Pipes)
```bash
# دمج أوامر
grep ">" genome.fasta | wc -l
grep "BRCA" genes.txt | cut -f1

# ترتيب وإزالة التكرار
cut -f1 annotations.txt | sort | uniq
```

### تمرين 4: إعادة التوجيه
```bash
# حفظ النتائج
grep ">" genome.fasta > gene_headers.txt
grep "exon" annotations.txt > exons_only.txt

# الإضافة لملف
echo "New analysis" >> results.txt
```

---

## 📚 ملاحظات

- جميع الملفات نصية بتنسيق UTF-8
- التسلسلات البيولوجية حقيقية (مقتطفات من قواعد بيانات NCBI)
- الأحجام صغيرة لسهولة التعلم (<10KB)
- آمنة للتحميل والمشاركة

---

## 🔗 روابط مفيدة

- [الفصل الثاني: أساسيات Linux](../../module-2.md)
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
- [FASTQ Format Specification](https://en.wikipedia.org/wiki/FASTQ_format)
- [FASTA Format Specification](https://en.wikipedia.org/wiki/FASTA_format)

---

**تم إعداده للتدريب على أساسيات Linux في المعلوماتية الحيوية**
مركز العلوم والبحث والتطوير - 2026
