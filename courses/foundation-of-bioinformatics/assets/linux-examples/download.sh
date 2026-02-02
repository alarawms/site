#!/bin/bash
# Download script for Linux training examples
# الفصل الثاني: أساسيات Linux

echo "📥 تحميل ملفات التدريب على Linux..."
echo ""

BASE_URL="https://malarawi.sa/courses/foundation-of-bioinformatics/assets/linux-examples"

# Create directory
mkdir -p linux-examples
cd linux-examples

# Download files
echo "⬇️  تحميل sequences.fasta..."
wget -q "${BASE_URL}/sequences.fasta"

echo "⬇️  تحميل genome_reads.fastq..."
wget -q "${BASE_URL}/genome_reads.fastq"

echo "⬇️  تحميل genome.fasta..."
wget -q "${BASE_URL}/genome.fasta"

echo "⬇️  تحميل genes.txt..."
wget -q "${BASE_URL}/genes.txt"

echo "⬇️  تحميل annotations.txt..."
wget -q "${BASE_URL}/annotations.txt"

echo "⬇️  تحميل README.md..."
wget -q "${BASE_URL}/README.md"

echo ""
echo "✅ تم التحميل بنجاح!"
echo ""
echo "الملفات المتوفرة:"
ls -lh

echo ""
echo "📖 اقرأ README.md للحصول على تعليمات الاستخدام:"
echo "   cat README.md"
echo ""
echo "🎯 جاهز للبدء بالتدريب!"
