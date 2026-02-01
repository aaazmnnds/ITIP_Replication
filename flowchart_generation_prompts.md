# Flowchart Generation Prompts
# These prompts were used to generate the ITIP flowcharts using AI image generation

## Figure 1: ITIP Framework Flowchart

**Prompt:**
```
Create a clean, professional flowchart diagram for the ITIP (Information-Theoretic Interaction Pruning) framework. The flowchart should show these steps in vertical sequence with arrows:

1. INPUT box: "Data with Missing Values (X)"
2. Step box: "Construct Missing Indicators Z"
3. Step box: "Impute Missing Values → X_imp"
4. Step box: "Generate Interactions I_j = X_imp_j × Z_j"
5. Decision diamond: "For each interaction j"
6. Step box: "Compute Conditional IG: IG(I_j | Z_j) using BAS"
7. Decision diamond: "IG(I_j | Z_j) ≥ ε?"
   - YES arrow → "Keep I_j"
   - NO arrow → "Prune I_j"
8. Step box: "Train Final Model with Pruned Features"
9. OUTPUT box: "Prediction Model"

Use professional colors: blue for process boxes, orange for decision diamonds, green for kept interactions, red for pruned. Include clear arrows showing flow. Modern, clean design suitable for academic publication.
```

**Output:** `figure1_itip_framework.png`

---

## Figure 2: Simulation Study Flowchart

**Prompt:**
```
Create a clean, professional flowchart for a simulation study showing two distinct phases:

PHASE 1 - DATA GENERATION (left side, light blue background):
1. "Generate X_full ~ N(0, Σ)"
2. "Generate Missingness Z ~ Bernoulli(0.3)"
3. "Generate Outcome Y ~ Logistic(X_full, Z)"
4. "Create X_incomp (mask missing values)"
Arrow labeled "Simulated Dataset" pointing to Phase 2

PHASE 2 - ITIP APPLICATION (right side, light green background):
5. "Impute: X_incomp → X_imp"
6. "Create Interactions: I_j = X_imp_j × Z_j"
7. "Compute IG(I_j | Z_j) for all j"
8. "Apply Threshold ε"
9. "Prune Low-IG Interactions"
10. "Results: Kept vs Pruned"

Use clear separation between phases with a vertical dashed line. Professional academic style with arrows showing flow. Label phases clearly.
```

**Output:** `figure2_simulation_flow.png`

---

## How to Regenerate

If you need to modify or regenerate these flowcharts:

1. **Using AI Image Generation Tools:**
   - Use the prompts above with tools like DALL-E, Midjourney, or similar
   - Adjust colors, layout, or text as needed

2. **Using Diagramming Software:**
   - **draw.io (diagrams.net)**: Free, web-based
   - **Lucidchart**: Professional diagramming
   - **Microsoft Visio**: Enterprise option
   - **TikZ (LaTeX)**: For publication-quality diagrams directly in LaTeX

3. **Using Python (Matplotlib/Graphviz):**
   - See `generate_flowcharts.py` (if you want programmatic generation)

---

## Files Generated

- `figure1_itip_framework.png` - ITIP Framework flowchart
- `figure2_simulation_flow.png` - Simulation Study flowchart

Both are referenced in `itip_jmlr.tex`:
- Figure 1: After Algorithm 1 (Section 4)
- Figure 2: In Section 6 (Illustrative Example)
