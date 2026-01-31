#!/bin/bash
# Test if the site builds correctly

cd /home/alarawms/org/site

echo "🔨 Testing MkDocs build..."
echo ""

# Check if in correct directory
if [ ! -f "mkdocs.yml" ]; then
    echo "❌ Error: mkdocs.yml not found!"
    echo "Run this script from /home/alarawms/org/site/"
    exit 1
fi

echo "✓ Found mkdocs.yml"
echo ""

# Try to build
echo "Building site..."
mkdocs build --verbose

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ BUILD SUCCESSFUL!"
    echo ""
    echo "Site built to: $(pwd)/site/"
    echo "Files generated: $(find site -type f | wc -l)"
    echo ""
    echo "This means your site WILL work on GitHub!"
    echo "The issue is just GitHub Actions configuration."
else
    echo ""
    echo "❌ BUILD FAILED!"
    echo ""
    echo "This is why GitHub deployment isn't working."
    echo "Fix the errors above, then push again."
fi
