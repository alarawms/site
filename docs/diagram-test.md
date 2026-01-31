# Diagram Format Test - Stacked Bilingual Labels

## Current Arabic Version
```mermaid
graph LR
    A[DNA<br/>التخزين] --> B[RNA<br/>النقل]
    B --> C[البروتين<br/>الفعل]
```

---

## Option 1: Simple Two-Line Stacked (Recommended)

```mermaid
graph LR
    A["DNA<br/>الحمض النووي"] --> B["RNA<br/>الحمض النووي الريبوزي"]
    B --> C["Protein<br/>البروتين"]
    C --> D["Function<br/>الوظيفة"]
    D --> E["Phenotype<br/>النمط الظاهري"]
```

**Format:**
- Line 1: English scientific term
- Line 2: Arabic translation
- Clean, minimal, easy to read

---

## Option 2: With Function Labels

```mermaid
graph LR
    A["DNA<br/>الحمض النووي<br/><small>Storage / التخزين</small>"] --> B["RNA<br/>الحمض النووي الريبوزي<br/><small>Transfer / النقل</small>"]
    B --> C["Protein<br/>البروتين<br/><small>Action / الفعل</small>"]
```

**Format:**
- Line 1: English term
- Line 2: Arabic translation  
- Line 3: Role/function (smaller text)

---

## Option 3: Arabic-First Order

```mermaid
graph LR
    A["الحمض النووي<br/>DNA"] --> B["الحمض النووي الريبوزي<br/>RNA"]
    B --> C["البروتين<br/>Protein"]
    C --> D["الوظيفة<br/>Function"]
    D --> E["النمط الظاهري<br/>Phenotype"]
```

**Format:**
- Line 1: Arabic translation (primary)
- Line 2: English scientific term (reference)

---

## Which do you prefer?

**My recommendation: Option 1** (simple, clean, follows scientific convention with English first)
