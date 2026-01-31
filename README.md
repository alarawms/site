# Bioinformatics and Data Science Website

Arabic-primary bilingual website for bioinformatics education and research documentation, built with MkDocs Material.

## Languages

- **Arabic** (Primary): Root level (`/`)
- **English** (Secondary): `/en/` subdirectory

## Features

- 📝 Bilingual blog system
- 🎓 Educational courses
- 📊 Interactive charts and diagrams
- 🖼️ Image galleries with lightbox
- 🔍 Full-text search (Arabic & English)
- 📱 Responsive RTL/LTR layout

## Local Development

### Prerequisites

- Python 3.x
- pip

### Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Serve locally
mkdocs serve

# Build site
mkdocs build
```

### Project Structure

```
.
├── docs/                    # Site content
│   ├── blog/               # Arabic blog (root)
│   ├── courses/            # Arabic courses
│   ├── en/                 # English content
│   ├── assets/             # Images, icons
│   ├── javascripts/        # Custom JavaScript
│   └── stylesheets/        # Custom CSS
├── mkdocs.yml              # MkDocs configuration
├── requirements.txt        # Python dependencies
└── .github/workflows/      # GitHub Actions for deployment
```

## Deployment

The site is automatically deployed to GitHub Pages on every push to the `main` branch via GitHub Actions.

## License

Content: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Built With

- [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
- [MkDocs Plugins](https://www.mkdocs.org/dev-guide/plugins/)
