# Bioinformatics Algorithms – 02-604 (Spring 2024)

This repository contains code written for the course **02-604: Fundamentals of Bioinformatics** (Spring 2024). The scripts implement core bioinformatics algorithms using Python (with some Go in Week 1), including global alignment, motif finding, hierarchical clustering, and suffix trees.

Each Jupyter Notebook corresponds to a weekly assignment or quiz. You can refer to the course syllabus and weekly quizzes here for context:  
🔗 [Course Syllabus and Quizzes](https://cogniterra.org/course/296/syllabus)

---

## Weekly Highlights

### 🔬 `week3_code` – Protein Fragmentation for Mass Spectrometry  
This script generates all possible contiguous substrings (peptides) from an input protein sequence. It simulates the type of fragmentation that occurs in mass spectrometry, helping to model peptide identification based on mass.

---

### 🧬 `week5_code` – Global Alignment Implementation  
Implements the **Needleman-Wunsch** algorithm for global sequence alignment. The code aligns two sequences by maximizing a scoring function that rewards matches and penalizes mismatches and gaps—critical for comparing DNA, RNA, or protein sequences.

---

### 🌳 `week8_code` – UPGMA for Phylogenetic Tree Construction  
Implements the **UPGMA (Unweighted Pair Group Method with Arithmetic Mean)** algorithm for hierarchical clustering. It builds a phylogenetic tree based on a distance matrix, assuming a constant molecular clock to infer evolutionary relationships.

---

### 🌿 `week10_code` – Suffix Tree Construction and DFS Traversal  
Constructs a **suffix tree** for a given string and performs **Depth-First Search (DFS)** traversal. Suffix trees are efficient data structures for fast pattern matching and are widely used in genome analysis and sequence comparison tasks.

---

Feel free to explore the notebooks to learn more about how each algorithm is implemented.
