## PRISM
[![Typing SVG](https://readme-typing-svg.demolab.com?font=Fira+Code&weight=500&size=24&duration=2000&pause=1000&color=B48EAD&width=550&lines=Screen.+Repair.+Interpret.+Repeat.)](https://git.io/typing-svg)
A CLI-based bioinformatics framework for gene variant pathogenicity screening and computational gene therapy candidate identification.
Integrates an AI interpreter to generated biologically grounded hypotheses based on PRISM result data, and proposes experimental follow-ups.

PRISM is designed as a research & prototyping tool, emphasizing result interpretability and understanding, reproducibility and future scope identification rather than autonomous decision-making.


### Overview
A CLI and result-specialized derivative of GEM (see my GitHub / project portfolio), paired with AI-enhanced result interpretability acheived through algorithm refactoring. 
PRISM addresses the gap in balancing classification accuracy and biological understanding and insight.

Note: All interpretations are explicitly constrained to remain hypothesis-driven and non-authoritative


### Key Features
- Mutation Fingerprint generation pipeline: computationally translates gene variants into numerical data capturing a wide range of feature changes across the following domains (w/ some examples):
  - Biochemical properties (thermodynamic stability, AA product hydrophobicity)
  - Structural effects (stop codons, AA composition changes)
  - Regulatory motif disruption (CREB, INR, Kozak)
  - Genetic domain changes (coding and non-coding % compositions)
  - Repeat Instability (expansions / contractions)
- High-performance Models: ClinicalModel and DiscriminatorModel both possess ROC and PR AUC metrics within the 90% range and 82% accuracy
  - ClinicalModel possesses the lowest false-negative rate, while DiscriminatorModel has the best balanced discriminative power (see GEM)
- Pathogenicity Screening
  - ML-based gene variant classification with probability outputs and mutation fingerprint recording
- Computational Gene Repair (ReGen)
  - Iterative generation of benign candidate variants
- Screening and Repair Mutation Fingerprint Recording, allows for AI-assisted interpretation
- AI-Assisted Data Interpretation
  - Explains predictions and the repair results using provided fingerprint records
  - Identifies dominant feature signals and explicitly states uncertainties and limitations
  - Proposes computational and wet-lab follow-up experiments (future scope)

### Project Structure
```bash
└───src
    └───prism
        ├───database
        │   └───pwm_database
        │       ├───AA_motifs
        │       └───DNA_motifs
        ├───DataSift_configs
        ├───gene_databank
        ├───model_folder
        │   ├───ClinicalModel
        │   └───DiscriminatorModel
        ├───ReGen_candidates
        │   ├───benchmark_repaired_RPs
        ├───Screen_results
        │   ├───benchmark_screened
        |───a02_1_CompositeDNA_Toolkit.py
        |───a02_2_CompositeProt_Toolkit.py
        |───a02_3_DNAMatrix_Toolkit.py
        |───a02_4_ProtMatrix_Toolkit.py
        |───a03_LookingGlass.py
        |───a04_ReGen.py
        |───b00_bio_library.py
        |───b01_utility.py
        |───cli.py
        |───config.yaml
        |───DiagnoSR.py
        |───Interpreter.py
        |───paths.py
        └───__init__.py
```

### Installation
PRISM is a research prototype intended to be run from source

prerequisites: 
- Python 3.10+
- pip and pipenv
- ollama + a local AI model, this version uses gpt-oss
  - if you decide to use another model, you'll have to change the model_biollm variable in Interpreter.py. I've labelled and commented around it for you near the top of the file, you can't miss it.

Go to the download ZIP in this repository and extract the ZIP file to a directory of your choice.

You should see the following files:
```python
Pipfile
Pipfile.lock
pyproject.toml
src/
```

Navigate to the project root
```python
cd path/to/PRISM/mine/is/Users/Elliot/BIOINFORMATICS FOLDER/PRISM
```

Install dependencies with:
```python
pipenv install
```

Install PRISM in editable mode:

```
pip install -e .
```
You can now access this program from anywhere in the terminal and go into the source file yourself and play around with the AI prompt or interpretation logic.

If you want to uninstall this from your computer you can type:
```
pip uninstall prism
```

### Workflow
This is an example run from my computer, note that 'benchmark.fasta' contains a fake gene I made purely for testing.

```commandline
PS C:\Users\Elliot> prism list gns
====[Browsing all available genes]====
benchmark.fasta
test.fasta
video_test.fasta
PS C:\Users\Elliot>
```

Let's say I want to screen benchmark.fasta
```commandline
PS C:\Users\Elliot> prism screen benchmark
[[Screening Initiated]]
===[[Initializing File Diagnostic]]===
selecting model: clinicmod
validating input fasta file: benchmark.fasta
found in gene_databank | now loading C:\Users\Elliot\BIOINFORMATICS FOLDER\PATHFINDER\src\prism\gene_databank\benchmark.fasta...
===[[Diagnostic Completed]]===

Screening variant file: C:\Users\Elliot\BIOINFORMATICS FOLDER\PATHFINDER\src\prism\gene_databank\benchmark.fasta
Loading ClinicalModel Sift: 8442 bytes
[Extracting highest prob. AA sequences]: 100%|██████████████████████████████████████████| 1/1 [00:00<00:00, 665.23it/s]
[Generating AA chain mutation fingerprints]: 100%|███████████████████████████████████████| 1/1 [00:00<00:00, 40.21it/s]
[Generating DNA mutation fingerprints]: 100%|████████████████████████████████████████████| 1/1 [00:00<00:00, 38.15it/s]
[Extracting broad domain changes]: 100%|█████████████████████████████████████████████████| 1/1 [00:00<00:00,  4.84it/s]
[Generating DNA motif fingerprints -- Position Weight Matrix Signals * Gaussian-weighted composite scoring + Cluster Co
[Generating AA profile fingerprints -- Regex & Position Weight Matrix Signals + Cluster Composite Scoring]: 100%|█| 1/1
PS C:\Users\Elliot>
```

Let's take a look at the result file it just made:

```commandline
PS C:\Users\Elliot> prism list scr
====[Browsing all screened genes]====
benchmark_screened
video_test_screened
PS C:\Users\Elliot> prism access scr benchmark_screened
             Name Predicted_Class Prob_Benign Prob_Pathogenic
0  benchmarkgene1               1  0.30137402        0.698626
PS C:\Users\Elliot>
```

I use the list function with scr for screened results and see that my benchmark gene has a 70% chance of being pathogenic. 
I can try and repair it with the following command:
```commandline
PS C:\Users\Elliot> prism repair benchmark --iterations 3
```
I've omitted the output here since it would take up a considerable amount of space.
I'll take a look at the repair log with the following commands:

```commandline
PS C:\Users\Elliot> prism list rpr
====[Browsing all repaired genes]====
benchmark_repaired_RPs
video_test_repaired_RPs
PS C:\Users\Elliot> prism access rpr benchmark_repaired_RPs
================================================================================
ReGen Analysis Results: ClinicalModel | benchmark | benchmarkgene1
================================================================================

ORIGINAL VARIANT STATS:
Ref Sequence: GCTGCTGGACCTGCC
Alt Sequence: AAAAAAAAAAAAAAAAAA
Benign % chance: 30.137402

ANALYSIS SUMMARY:
|- Starting Score: 0.301374
|- Original Length: 15 bp
|- Final Variants: 1
|- Benign Threshold Variants: 0
|- ReGen config: 3 iterations, 1 copies

MAX BENIGN VARIANTS PER ITERATION:
--------------------------------------------------
Score: 50.294387340545654 | Length: 21 bp
Benign % increase: 20.15698552131653
   Sequence:
    AAAAAAAAAAAAAAAAAATCA

Score: 73.50447177886963 | Length: 24 bp
Benign % increase: 43.3670699596405
   Sequence:
    AAAAAAAAAAAAAAAAAATCATCA

Score: 73.50447177886963 | Length: 24 bp
Benign % increase: 43.3670699596405
   Sequence:
    AAAAAAAAAAAAAAAAAATCATCA

FINAL VARIANTS:
--------------------------------------------------
Score: 73.50447177886963 | Length: 24 bp
Benign % increase: 43.3670699596405
   Sequence:
    AAAAAAAAAAAAAAAAAATCATCA
```

So now we have both screened and repaired benchmarkgene1, using all necessary commands! We can now leverage our local AI interpreter on either of the result directories.
I'm going to ask it to interpret the repair results 
```commandline
PS C:\Users\Elliot> prism interpret rpr benchmark_repaired_RPs
Loading Interpreter...
Interacting with Interpreter...
**Interpretation of the Regen repair output**

| Metric | Original | Variant (24 bp) | ΔBenign % |
|--------|----------|-----------------|-----------|
| Score | 0.301374 | 73.504471 | +43.4 % |
| Length | 15 bp | 24 bp | +9 bp |

The most striking change in the predictive model is the jump from a 30 % benign probability to ~73 % benign probability.  This rise is reflected in the fingerprint features that the Regen algorithm maximises (see the “Composite” and “Boundary” fingerprints).  Below are the fingerprint changes that most strongly accompany this shift, the hypotheses for how they could influence the prediction, and the ambiguities that remain.
---
```
I've omitted the rest due to the size of the prompt. 

For a full breakdown of all workflow parameters, type 'prism -h' in the command line.

That concludes the test run!


### Limitations & Notes
- **Research Prototype**: PRISM is a prototype for research and educational purposes, not for clinical decision-making
- **AI Interpretations**: AI-generated hypotheses require expert validation and experimental verification
- **Local AI Performance**: Interpretation generation speed depends on your hardware and chosen model (50B parameters or 20B)
- **Model Flexibility**: While designed for gpt-oss, you can swap in other ollama-pulled models by modifying `Interpreter.py`


Shield: [![CC BY-NC 4.0][cc-by-nc-shield]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg


#### Questions or Collaboration?
- e224chan@uwaterloo.ca
- elliotchan119@gmail.com