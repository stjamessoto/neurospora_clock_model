---
title: "Understanding the Neurospora Clock Network Project"
subtitle: "A Plain-Language Guide to the Data, the Methods, and Why the Numbers Come Out the Way They Do"
author: "Neurospora Clock Network Project"
date: "2026"
---

# Executive Summary

This project takes the fitted regulatory network from Al-Omari et al. (2022) — 11
regulator proteins controlling roughly 3,000 downstream genes in the circadian clock
of the bread mold *Neurospora crassa* — and subjects it to statistical tests the
original paper never ran, while being explicit about where and why this project's own
numbers diverge from the paper's.

**What was reproduced.** The core network construction, the feedforward-loop counting
formula, the degree/hub-size statistics, and the function-level binding-strength table
were all rebuilt directly from the paper's methods and cross-checked against its
reported numbers (see the validation table in Part 5.2 below). Some numbers match
closely or exactly (e.g. the smallest hub, rok-1's 54 targets); others diverge
substantially.

**The central finding.** Every substantial divergence from the paper — fewer
multiply-regulated genes, a lower feedforward-loop count, and the complete absence of
any regulator pair, triple, or quartet sharing *more* targets than chance — traces back
to one fact: the exported connection matrix records only each gene's single *dominant*
regulator, discarding the full continuous ensemble estimate behind it. Because of this,
no gene in this export is ever controlled by more than 2 regulators at once, which makes
it mathematically impossible for any 3 or more regulators to ever share a target gene in
this data. This one structural property of the exported data — not any biological claim
about the real regulatory network — explains essentially every "diverges" row in the
validation table.

**Novel statistical extensions, and what they found:**

- A degree-preserving randomization test (1,000 rewired networks) shows the real
  feedforward-loop count (30) is statistically indistinguishable from pure chance
  (Z = -0.41, two-sided p = 0.775) — a materially different conclusion than the paper's
  own qualitative comparison to an unrelated organism.
- Exact statistical tests across all 550 possible regulator combinations (55 pairs, 165
  triples, 330 quartets) find 60 significant results after correcting for running so
  many tests at once — and every one of them is a *depletion* (regulators sharing
  *fewer* targets than chance), never an enrichment. More strikingly, every one of the
  495 triples and quartets has exactly zero observed shared targets — a structural
  certainty, not merely a statistical pattern.
- A from-scratch investigation of a separate, previously unexplained simulated dataset
  (`MuThetaDataSim_m4.txt`) confirmed what its 10-value "Mu" block represents
  (conditional ensemble averages over 2,000 MCMC runs) and definitively ruled out the
  leading hypothesis for its 6-value "Theta" block (best possible accuracy of 28.7%,
  worse than the trivial 56.4% baseline of always guessing the most common regulator).
- A sequence-based search across the full *Neurospora* proteome found that none of the
  11 clock regulators themselves have a detected evolutionary "sibling" gene (paralog),
  though several of their target genes do.

**What the rest of this document is.** Everything below walks through all of this in
plain language for a reader with no statistics or biology background — defining every
term before using it, explaining every test with a concrete analogy, and being explicit
throughout about what is confirmed, what is inferred, and what remains genuinely open.

---

# About This Document

This document exists for one reason: everything else in this project (the code, the
dashboard, the README) is written for someone who already knows what a p-value is,
what a transcription factor does, and why anyone would care about a fungus's internal
clock. This document is written for everyone else.

If you have never taken a statistics class and have never studied molecular biology,
you should still be able to read this document start to finish and come away
understanding: what this project is, what data it uses, what questions it asks, how it
answers them, and — most importantly — *why the results come out the way they do,
especially the surprising ones*. Every technical term is defined the first time it is
used. Every statistical test is explained with a plain-language analogy before the
formal version. Nothing is assumed.

If you only read one section, read **Part 5.2, "Why is everything depleted?"** — it is
the single most important finding in this project, it explains a pattern that shows up
in three different analyses, and it is explainable without any statistics background
at all.

---

# Part 1: The Biology, in Plain Terms

## 1.1 What is *Neurospora crassa*, and why do scientists study its clock?

*Neurospora crassa* is a species of bread mold — an orange fungus that grows on stale
bread and, more usefully for science, grows extremely well in a laboratory. Scientists
have used it as a "model organism" (a simple, easy-to-study stand-in for more
complicated living things) for over 80 years, because its genetics are well understood
and it does something very convenient for clock researchers: it visibly changes its
growth pattern on a roughly 24-hour cycle, even when kept in constant darkness with no
external time cues at all.

That 24-hour internal rhythm is a **circadian clock** ("circadian" literally means
"about a day," from Latin *circa diem*). Humans have one too — it is why you feel
sleepy at roughly the same time every night and why jet lag is unpleasant. Every
circadian clock, in every organism from fungi to humans, works on the same basic
principle: a small set of genes produce proteins that, after a delay, shut off their
own production, then get cleared out, then the whole cycle restarts. That
delayed-negative-feedback loop is what makes the rhythm self-sustaining instead of a
one-time pulse.

Because *Neurospora*'s clock is genetically simple and easy to manipulate, it has been
one of the primary systems used to work out *how* circadian clocks work at the
molecular level in general — discoveries made in this mold have repeatedly turned out
to apply, with modifications, to the clocks of plants, insects, and mammals.

## 1.2 The core clock, and the much bigger question this project is about

At the center of *Neurospora*'s clock are three genes, usually written as *frq*,
*wc-1*, and *wc-2* (short for "frequency" and "white collar 1/2"). Their proteins form
the delayed feedback loop described above: WC-1 and WC-2 join together to form a
protein pair called the **White Collar Complex**, or **WCC** for short. WCC turns on
the *frq* gene; the FRQ protein then, after building up, turns off its own gene
(indirectly, by interfering with WCC) — and because FRQ eventually gets broken down,
the whole cycle restarts roughly every 24 hours. WCC also directly senses light, which
is how the clock resets itself to match the real day/night cycle instead of drifting.

That three-gene loop is the clock's "engine." But an engine on its own does not do
anything useful — a clock is only useful to an organism if it *controls other things*:
telling the organism when to grow, when to release spores, when to ramp up particular
metabolic pathways, and so on. The genes whose activity rises and falls on the clock's
schedule, *because* WCC and its downstream partners turn them on and off, are called
**clock-controlled genes**, sometimes shortened to **CCGs**. This project — and the
paper it is built on — is entirely about mapping *which* genes are clock-controlled,
and *which* piece of the clock machinery controls each one.

## 1.3 "Regulators" and "target genes"

Two words are used constantly throughout this project, so it is worth being precise
about them from the start:

- A **regulator** is one of a small number of proteins (11 of them, in this project's
  data — WCC plus 10 others) that sit downstream of the core clock and control the
  on/off, up/down activity of other genes.
- A **target gene** is one of the (roughly 3,000) other genes in the genome whose
  activity is controlled by one or more of those 11 regulators.

So the overall picture is a hierarchy: the 3-gene core clock (*frq*/*wc-1*/*wc-2*)
drives WCC, WCC drives the other 10 regulators, and all 11 regulators together drive
about 3,000 target genes. This project is about that last layer — the regulators and
their roughly 3,000 targets — not about the core clock mechanism itself, which is taken
as a known, fixed starting point.

## 1.4 Two different *kinds* of regulator

Not all 11 regulators control their target genes the same way. There are two
categories, and the distinction matters for interpreting several results later in this
document:

- **Transcriptional regulators** (5 of the 11, including WCC) are classic
  "transcription factors" — proteins that physically bind to DNA near a target gene and
  turn its production up or down. This is the textbook picture of gene regulation that
  most people learn in a biology class.
- **Post-transcriptional regulators**, also called **RNA operons** (6 of the 11), work
  completely differently: instead of touching DNA, they bind directly to the *RNA
  message* that has already been copied from a gene, and control how long that message
  survives before being broken down. A gene whose message is protected/stabilized gets
  effectively "louder" (more protein made from it) without its DNA-level control
  changing at all. The name "RNA operon" is borrowed from a related idea in yeast
  biology: one RNA-binding protein can simultaneously control a whole *set* of
  otherwise-unrelated genes, as if they were wired together — a kind of coordination
  that happens after the DNA-reading step, not before it.

The table below lists all 11 regulators used in this project, their type, and (where
known from the original paper) what kind of biological process each one is associated
with.

| Regulator (common name) | Type | What it's associated with (per the original paper) |
|---|---|---|
| **WCC** | Master regulator (transcriptional) | Senses light; drives the core clock; activates all 10 other regulators |
| **pro-1** | Transcriptional | Ribosome/ribosome biogenesis, meiosis, recombination |
| **rpn-4** | Transcriptional | Carbon metabolism, amino acid metabolism, glycerophospholipid metabolism, oxidative phosphorylation; strong positive light response |
| **NCU06108** | Transcriptional | Strong *negative* light response; shares many ribosome-related targets with lhp-1 |
| **NCU07155** | Transcriptional | Nitrogen and sulfur metabolism |
| **lhp-1** | Post-transcriptional (RNA operon) | The single largest hub in the whole network; ribosome, ribosome biogenesis, proteasome, protein export, amino acid biosynthesis, tRNA processing |
| **rrp-3** | Post-transcriptional (RNA operon) | Ribosome, ribosome biogenesis, RNA transport, carbon metabolism, amino acid biosynthesis |
| **NCU03363** | Post-transcriptional (RNA operon) | One of six RNA operons in the paper's model; specific function not detailed |
| **pab-1** | Post-transcriptional (RNA operon) | Poly-A binding protein; uniquely associated with ketone body synthesis/degradation genes |
| **rok-1** | Post-transcriptional (RNA operon) | The smallest hub in the network |
| **has-1** | Post-transcriptional (RNA operon) | A DEAD-box RNA helicase; ribosome biogenesis, meiosis, RNA transport, the spliceosome |

("NCU" numbers, like NCU06108, are standardized gene ID codes for *Neurospora crassa*
— every gene in the organism's genome has one, similar to how every book has an ISBN.
Some regulators also have a shorter, more familiar name like *lhp-1* or *rok-1*; some
don't, and are referred to by their NCU number throughout.)

---

# Part 2: The Original Study This Project Builds On

This entire project is a re-analysis and extension of one published scientific paper:

> Al-Omari, A. M., Griffith, J., Scruse, A., Robinson, R. W., Schüttler, H.-B., &
> Arnold, J. (2022). Ensemble Methods for Identifying RNA Operons and Regulons in the
> Clock Network of *Neurospora Crassa*. *IEEE Access*, 10, 32510-32524.
> https://doi.org/10.1109/ACCESS.2022.3160481

## 2.1 What the original researchers did

The original researchers took real, measured data — 60,759 individual measurements of
messenger RNA levels over time, across the genome — and used it to reconstruct, for
each of the roughly 3,000 target genes, *which* of the 11 regulators control it. This
is a much harder problem than it sounds: you cannot simply look at two genes' activity
rising and falling together and conclude that one controls the other, because there are
thousands of genes moving up and down on similar schedules for many unrelated reasons.
Instead, the researchers built a statistical model with an adjustable "wiring diagram"
(which regulator controls which gene) and adjustable "strengths" for each connection,
and used a technique called **MCMC** (Markov Chain Monte Carlo — explained in the
glossary, Part 4) to search through millions of possible wiring diagrams and find the
ones that best explain the real measured data.

Crucially, they didn't just run this search once and take a single answer. They ran it
as an **ensemble**: thousands of independent runs of the same search procedure, each
one converging on a slightly different candidate wiring diagram, all of them consistent
with the real data to varying degrees. Averaging and summarizing across that whole
ensemble — rather than trusting any single run — is what makes the result a defensible
statistical estimate instead of one arbitrary guess. Their specific method for doing
this ensemble search is called **VTENS** ("Variable Topology Ensemble Method"), and it
is where the phrase "ensemble" throughout this project's documentation comes from.

## 2.2 What came out of that process, and what this project starts from

The end product of the original paper's analysis, exported for this project as
`Connection_Matrix.xlsx`, is a big table: one row per regulator (11 rows), one column
per candidate target gene (originally 3,402 columns, cleaned down to 3,026 — see Part
5.1), and a 1 or a 0 in each cell, meaning "this regulator does / does not control this
gene, according to the fitted model." This 1/0 table — which this project calls the
**connection matrix**, or **T** for short in the code — is the starting point for
almost everything in this project.

It is important to understand what that 1/0 table actually is, because it explains
several surprising results later in this document. The original ensemble search doesn't
naturally produce a clean yes/no answer — it produces, for every (regulator, gene)
pair, a continuous **binding strength** estimate (roughly: "how strongly, on average
across the ensemble, does this regulator seem to control this gene"). To get from that
continuous ensemble output to the clean 1/0 table, the researchers had to make a
decision rule: for each gene, look at all 11 regulators' binding-strength estimates and
call the largest one — the **dominant regulator** — the "real" controller, encoded as a
1, with everyone else getting a 0. `Connection_Matrix.xlsx` is a *thresholded snapshot*
of the full continuous ensemble result, not the full result itself. This one fact — a
continuous estimate collapsed down to "pick the single biggest winner per gene" — turns
out to be the root cause of nearly every place where this project's numbers diverge
from the original paper's own reported numbers (see Part 5.2).

## 2.3 What the original paper concluded

Using this reconstructed network, the original paper reported several headline
findings, most of which this project attempts to reproduce and check in Part 5:

- A largest hub of 768 target genes (belonging to *lhp-1* — bigger than the master
  regulator WCC itself).
- A smallest hub of 54 target genes (belonging to *rok-1*).
- 71 **feedforward loops** (a specific 3-node wiring pattern explained in Part 4) in
  the network.
- A degree distribution (how many regulators control each gene, and how many targets
  each regulator has) roughly following a mathematical pattern called a power law,
  with an average of 2.23 connections per gene.
- The paper's evidence that the 71 feedforward loops were "meaningful" was a single
  comparison to a textbook number from a completely different organism (baker's
  yeast/*E. coli*: "71 is more than the 40 reported there"), not a statistical test run
  on this network's own data.

That last point is the seed of one of this project's two biggest novel contributions
(Part 5.3): the original paper never asked "is 71 more feedforward loops than we'd
expect by pure chance, given how big this specific network's hubs are?" — it only
compared to an unrelated organism. This project runs that missing statistical test.

---

# Part 3: The Data Files This Project Works With

This project uses five real, user-supplied data files, plus the outputs it computes
from them. Understanding what each one is (and, in one case, what it is *not*) matters
for trusting the results downstream.

## 3.1 `Connection_Matrix.xlsx` — the regulator-to-gene wiring diagram

Described above in Part 2.2: 11 regulators x (originally) 3,402 candidate target-gene
columns, thresholded 1/0 values. Before this data could be used, it needed real
cleaning:

- Some row and column labels had stray trailing whitespace (e.g. a regulator name with
  an invisible extra space at the end) that would otherwise silently create duplicate,
  disconnected entries for what should be a single gene or regulator.
- The spreadsheet-reading software automatically appends `.1`, `.2`, etc. to repeated
  column headers rather than merging them, which happens because the same gene
  legitimately appears more than once in the raw export. Un-merging those (recognizing
  `"NCU06108.1"` as really being `"NCU06108"`) recovered 3,029 distinct genes out of
  the original 3,402 raw columns, and combining (via logical OR — "does *any* copy of
  this column say the regulator controls this gene?") the repeated columns for the same
  gene is what produces a clean, non-redundant table.
- Three columns turned out not to be downstream target genes at all, but the three core
  clock genes themselves (*wc-1*, *wc-2*, *frq* — see Part 1.2), included in the raw
  export with odd positional-suffix names. Since the core clock is the *input* to this
  network, not something the regulators control, these three are set aside separately
  rather than counted as targets.
- After all cleaning, the final usable target-gene set is **3,026 genes** (down from
  3,402 raw columns; the original paper reports 3,380).

## 3.2 `binding_strength_by_function.csv` — a real look at the continuous estimates

As explained in Part 2.2, `Connection_Matrix.xlsx` only shows the *winner* per gene —
the single dominant regulator — and throws away the actual continuous binding-strength
numbers behind that decision. This second file, transcribed from a supplementary table
in the original paper, restores a little of that lost detail: for 21 broad functional
categories (ribosome, proteasome, carbon metabolism, and so on — standard categories
from public gene-function databases called KEGG and GO), it gives the continuous
binding-strength estimate for all 11 regulators.

This file needs one important caveat, spelled out clearly in the dashboard itself: a
`0.0000` value in this table does **not** mean "this regulator's effect was measured
and found to be exactly zero, with statistical confidence." It means "this regulator
was not the single largest (dominant) one for this function" — exactly the same
winner-takes-all collapsing described in Part 2.2, just applied at the level of a
functional category instead of an individual gene. There is no way, from this data
alone, to test whether a *non-dominant* regulator's true effect on a given function is
meaningfully different from zero — the ensemble search that produced these numbers was
never designed to answer that question for anything but the single dominant winner.

## 3.3 `MuThetaDataSim_m4.txt` and `ALLNCU.txt` — a per-gene simulated dataset

These two files travel together. `ALLNCU.txt` is a plain list of 2,418 gene ID codes
(2,216 of them distinct — some genes appear more than once). `MuThetaDataSim_m4.txt`
holds 2,418 numeric records in the *exact same order* — confirmed directly against the
data source, not just assumed because the two files happen to be the same length — so
that line 1 of `ALLNCU.txt` and record 1 of `MuThetaDataSim_m4.txt` describe the same
gene, line 2 with record 2, and so on.

Each record has two blocks of numbers:

- **"Mu"** — 10 continuous values per gene, roughly in the range 0-100. These are
  **conditional ensemble averages**, confirmed directly from the data source: recall
  from Part 2.1 that the original binding-strength estimate came from thousands of
  independent MCMC ensemble runs (2,000, specifically). For a given gene and a given
  regulator, most of those 2,000 runs will *not* have that regulator as the gene's
  dominant (largest-binding-strength) regulator — but some subset will. A Mu value is
  the *average* binding strength of one regulator, computed *only* over the ensemble
  runs in which that regulator happened to be dominant for that gene — not a simple
  average over all 2,000 runs. This is a meaningfully different (and more informative)
  number than a plain mean would be, because it answers "when this regulator does win,
  how strongly does it tend to win?" rather than diluting that signal with all the runs
  where it didn't win at all. What is *not* yet confirmed is which of this project's 10
  non-WCC regulators each of the 10 Mu columns corresponds to positionally — the
  general meaning of the numbers is now known, but the specific column-to-regulator
  labeling is still an open question.
- **"Theta"** — 6 non-negative values per gene that always sum to exactly 1.0 (verified
  across every one of the 2,418 records, not assumed). A set of numbers that are all
  non-negative and sum to 1 is, mathematically, a *probability distribution* — Theta is
  some kind of 6-way weighting or likelihood split for each gene. This project tested
  the most obvious hypothesis for what those six categories might represent: the six
  real post-transcriptional ("RNA operon") regulators from Part 1.4. That hypothesis
  was tested rigorously (every one of the 720 possible ways to match up the 6 Theta
  columns to the 6 real regulators was tried) and **rejected** — see Part 5.6 for the
  full result. What Theta actually represents remains an open question.

## 3.4 The UniProt protein-sequence file

A large reference file (10,255 protein sequences) covering essentially every gene in
the *Neurospora crassa* genome, downloaded from UniProt (a major public protein
database). This is the one file in the project that is about the actual physical
sequence of each gene's protein — letters representing amino acids — rather than about
measured activity levels. It is used for exactly one purpose: finding genes that are
evolutionary "siblings" of each other (Part 5.7).

## 3.5 What this project computes and caches

Several analyses in this project (described in Part 5) take real time to compute — the
slowest, a randomization test, takes 8-9 minutes. Rather than re-running them every
time someone opens the dashboard, their results are computed once by a script and saved
to reusable files (in `data/processed/`), which the interactive dashboard then simply
reads. This is a standard practice for keeping an interactive tool responsive without
sacrificing reproducibility — anyone can re-run the original scripts from scratch and
regenerate the exact same cached files.

---

# Part 4: A Glossary of Every Technical Term Used in This Project

This section defines every statistical and technical concept referenced anywhere in
this project, in the order they become useful, each with a plain-language analogy.
Refer back to this section as needed while reading Part 5.

**Ensemble.** A collection of many independent attempts at the same task, whose results
are combined rather than trusted individually — like asking 2,000 different people to
independently estimate the weight of an object and averaging their guesses, instead of
trusting just one guess. Both the original paper's MCMC search (Part 2.1) and this
project's own randomization tests (Part 5.3) are "ensembles" in this sense.

**MCMC (Markov Chain Monte Carlo).** A computational search technique for exploring an
enormous space of possibilities (here: which of the 11 regulators controls each of
~3,000 genes, and how strongly) when it's too large to check every option directly. It
works by starting somewhere and repeatedly making small random adjustments, keeping
adjustments that make the fit to real data better (and sometimes, briefly, accepting
ones that make it slightly worse, to avoid getting permanently stuck on a
so-so answer) — over many iterations this process tends to settle into the regions of
possibility-space that best explain the real data.

**Dominant regulator.** For a given gene, whichever of the 11 regulators has the
largest estimated binding strength for that gene. Only the dominant regulator gets
recorded as a "1" in `Connection_Matrix.xlsx`; every other regulator gets a "0" for
that gene, regardless of how close its own binding-strength estimate was. See Part 2.2.

**Hub / hub size.** A "hub" is a regulator considered together with all of the target
genes it controls, treated as one unit. Its "size" is simply how many target genes it
controls — lhp-1, with 739-768 targets depending on which export you look at, is the
largest hub in this network; rok-1, with 54, is the smallest.

**Degree.** How many connections a single node has. For a *regulator* node, its degree
is the same thing as its hub size. For a *target gene* node, its degree is how many of
the 11 regulators control it — in this project's cleaned data, every gene's degree is
either 1 or 2 (see Part 5.2 for why this matters enormously).

**Feedforward loop (FFL).** A specific 3-node wiring pattern: WCC controls some
regulator "B," and *both* WCC *and* regulator B control the same target gene "C." Named
"feedforward" because the signal reaches gene C by two paths at once (directly from WCC,
and indirectly via B) — a well-studied pattern in biological networks generally,
thought to help filter out brief noise or introduce a deliberate timing delay in how a
target gene responds.

**Null model.** A deliberately "boring," randomized version of your data, built to
share some basic structural properties with the real data (in this project: the same
number of targets per regulator, and the same number of regulators per gene) but with
the specific *wiring* scrambled at random. The point of a null model is to answer the
question "is what I observed in the real data actually special, or is it just what you'd
expect from *any* network with these same basic size properties, wired up randomly?"
This project builds 1,000 such randomized networks (Part 5.3) to test whether the real
FFL count is unusual.

**Z-score.** A single number describing how far an observed value is from what a null
model predicts, measured in units of the null model's own natural variability
(technically, its standard deviation). A Z-score near 0 means "right in line with random
chance." A Z-score whose absolute value is bigger than roughly 2 is generally considered
notably far from chance, in either direction (a positive Z-score means "much higher than
random chance would predict"; a negative one means "much lower").

**p-value.** The probability, if there really were no special effect at all (only random
chance), of seeing a result at least as extreme as the one actually observed. A small
p-value (conventionally, below 0.05) is usually read as "this would be a surprising
coincidence if nothing real were going on, so probably something real is going on." A
p-value is *not* the probability that the finding is "true," and is not, by itself, a
measure of how *large* or important an effect is — only of how *surprising* it would be
under pure chance.

**Hypergeometric test.** The standard statistical test for a very specific, very common
situation: you have a population divided into two groups (say, genes that are a
particular regulator's targets, and genes that aren't); you draw another, separate
subset from that same population (a second regulator's target set); and you want to know
whether the two subsets *overlap* more, or less, than pure chance would predict, given
how big each one is. The classic teaching analogy is drawing colored balls from an urn
without looking: if you know how many red balls and how many total balls are in an urn,
and you draw a handful, the hypergeometric test tells you exactly how surprising any
particular number of red balls in your handful would be. This project's enrichment tests
(Part 5.4) are exactly this test, run on regulators' target-gene sets instead of colored
balls.

**Multiple testing / Benjamini-Hochberg FDR correction, and "q-value."** If you run one
statistical test at the conventional "p < 0.05 counts as significant" threshold, there's
a 5% chance of a false alarm purely by chance. If you run 550 *independent* tests at that
same threshold (as this project's enrichment analysis does, Part 5.4), pure chance alone
would be expected to produce roughly 28 "significant-looking" false alarms — not because
anything real is happening, just from running so many tests at once. The
**Benjamini-Hochberg** procedure is a standard, well-established correction for this
problem: instead of asking "is this one p-value below 0.05," it adjusts every result
into a **q-value**, such that "keep everything with q <= 0.05" guarantees that, of
*everything you keep*, no more than roughly 5% is expected to be a false alarm — no
matter how many tests you ran in total. This project always reports and filters on
q-values, never raw p-values, for exactly this reason.

**Paralog / gene family.** Two genes are "paralogs" of each other if they descend from
a shared ancestor gene that duplicated at some point in evolutionary history, and both
copies are still recognizably similar today. A "family" is a group of such related
genes; its "size" is how many genes belong to it. Finding paralogs requires comparing
actual protein sequences (letter-by-letter, roughly) — it cannot be determined from
activity/expression data alone, which is why this project's paralog analysis (Part 5.7)
needed the separate UniProt sequence file from Part 3.4.

---

# Part 5: What This Project Adds, One Analysis at a Time

Everything up to this point describes background and input data. This part walks
through every original analysis this project performs — what question it asks, how it
answers it, and what the answer turned out to be, explained without assuming any prior
statistics background.

## 5.1 Cleaning the raw data

Covered in Part 3.1. The short version: the raw spreadsheet had whitespace glitches and
mangled duplicate columns that, if used as-is, would have silently under-counted how
many genes each regulator actually controls. After proper cleaning, this project works
with a clean table of 11 regulators by 3,026 target genes.

## 5.2 Validating against the paper — and why is everything depleted?

This is the section to read if you only read one.

The most natural first check for a project like this is: does re-deriving basic numbers
from the exported data match what the original paper itself reported? The table below
shows that check.

| Quantity | Paper's reported value | This project's cleaned export | Verdict |
|---|---|---|---|
| rok-1 target count | 54 | 54 | Exact match |
| lhp-1 target count | 768 | 739 | Close, not exact |
| Singly-regulated genes (targeted by exactly 1 regulator) | 2,686 | 2,797 | Close, not exact |
| Multiply-regulated genes (targeted by 2+ regulators) | 694 | 229 | **Diverges substantially** |
| Total clock-controlled genes | 3,380 | 3,026 | Close, not exact |
| Maximum regulators controlling any one gene | up to 11 (implied by the model) | **2** | **Diverges substantially** |
| Total feedforward loops | 71 | 30 | Diverges |
| Genes with *overlapping* feedforward loops | 5 | 0 | Diverges |
| Average node degree | 2.23 | 2.14 | Close |
| Power-law exponent of the degree distribution | -1.4 | -1.16 | Close |

Some of these are close, some diverge sharply — and it would be easy to treat this as a
list of separate small discrepancies. It is not. **Every one of these divergences traces
back to the single structural fact already introduced in Part 2.2: `Connection_Matrix.xlsx`
only records each gene's single *dominant* (winning) regulator, discarding every other
regulator's estimate for that gene, no matter how close the contest was.**

Here is the plain-language version of why that one fact explains everything else.
Imagine 100 people each vote for their favorite of 11 candidates in an election, and
you are only ever told the name of the single winner in each of 100 separate rooms —
never the full vote tally, never who came second. From that winners-only data alone:

- You would never see two candidates "tied" or "both winning" in the same room — by
  construction, a room has exactly one recorded winner. This directly explains why the
  maximum number of regulators recorded for any one gene, in this cleaned export, is 2
  and not up to 11 — actually, thinking it through further, it explains something even
  stronger: a gene can be recorded with 2 regulators only when WCC (which activates
  every regulator unconditionally, by the model's fixed design — see Part 1.2) and one
  other regulator both "win," since WCC's edge to every other regulator is fixed and not
  part of the same winner-take-all contest as the ~3,000 target genes. For the target
  genes themselves, each one's controlling regulator set is fundamentally a single
  winner, which is why "multiply-regulated" counts crater compared to the paper's
  reported number, which came from the *full* continuous ensemble estimate, not the
  thresholded winners-only snapshot.
- Because fewer genes end up with 2 regulators instead of many, there is less
  "raw material" available for feedforward loops (Part 4's definition requires a gene to
  be controlled by *both* WCC and a second regulator at once) — which is exactly why the
  feedforward-loop count comes out lower (30 vs. the paper's 71), and why *zero* genes
  show the even-more-demanding "overlapping" pattern (controlled by WCC and 2+ other
  regulators at once), where the paper found 5.
- And directly relevant to this project's own novel enrichment tests (Part 5.4): if
  every gene's controlling-regulator-set is fundamentally a single winner, then any two
  regulators' target-gene sets are automatically pushed toward being separate,
  non-overlapping sets — not because the two regulators biologically compete with or
  suppress each other, but purely as an artifact of only ever recording one winner per
  gene. **This is why every statistically significant result in this project's
  enrichment analysis, without a single exception, comes out as "these regulators share
  *fewer* target genes than pure chance would predict" (a "depletion") and never as
  "these regulators share *more* target genes than chance" (an "enrichment").** It is a
  property of how the exported data was constructed, not evidence of any real biological
  antagonism between regulators. Section 5.4 below goes even further: it turns out that
  for *any three or more* regulators at once, the shared-target overlap in this
  particular export isn't just statistically unlikely — it is *exactly zero, always, with
  no exceptions*, for a reason that follows directly and mechanically from this same
  single-winner-per-gene fact.

None of this means the original paper's numbers are wrong, or that this project's
numbers are wrong. It means the exported `Connection_Matrix.xlsx` table is a simplified,
thresholded summary of a richer underlying result, and every downstream analysis run on
that simplified table inherits its simplification. This project treats that as an
important, honestly reported property of the data — not something to quietly work around
or paper over.

## 5.3 Is the feedforward-loop count actually meaningful?

The original paper's only evidence that "71 feedforward loops" was a meaningful,
interesting number was a comparison to a *completely different organism* (a textbook
count of 40 in *E. coli*/yeast literature). It never asked the more relevant question:
**given the sizes of *this specific network's* regulator hubs, how many feedforward
loops would show up even in a totally random rewiring of the same network?** If the
answer is "about the same as what was actually observed," then the loop count isn't
telling you anything special about this network's wiring — it's just a predictable
side effect of how big the hubs happen to be.

This project answers that question directly with a **null model** (defined in Part 4):
it builds 1,000 randomly rewired versions of the network, each one carefully constructed
to keep exactly the same properties fixed as the real network (every regulator keeps the
exact same number of target genes it really has; every gene keeps the exact same number
of controlling regulators it really has) while scrambling *which* regulator connects to
*which* gene. It then counts feedforward loops in every single one of those 1,000
randomized networks, giving a full distribution of "how many loops would you see by pure
chance, given only the sizes involved."

| Quantity | Value |
|---|---|
| Observed feedforward-loop count (this export) | 30 |
| Average count across 1,000 randomized networks | 32.01 |
| Spread (standard deviation) across those 1,000 networks | 4.94 |
| Z-score | -0.407 |
| Chance of seeing this few or fewer loops randomly (one-sided) | 69.0% |
| Chance of seeing a result this extreme in either direction (two-sided) | 77.5% |

**The finding: the real feedforward-loop count (30) is statistically indistinguishable
from pure chance.** A Z-score of -0.407 is far below the "notably unusual" threshold of
roughly 2 in absolute value, and a two-sided p-value of 0.775 means a result this
ordinary would happen about 3 times out of 4 in a purely random network with the same
basic hub sizes. This is a genuinely different conclusion than the paper's own implied
claim (that 71 loops is a meaningfully large, significant number) — though it is fully
consistent with, and arguably a direct consequence of, the structural fact from Part
5.2: with every gene capped at 2 controlling regulators, there simply isn't much room
for the loop count to be anything other than a predictable, mechanical side effect of
how many targets each regulator has, which is exactly the property this null model was
built to hold fixed and control for.

## 5.4 Do regulators avoid or prefer sharing target genes?

This is the project's largest novel statistical analysis. The question: for any group of
2, 3, or 4 of the 11 regulators, do they control significantly *more* genes in common
than you'd expect by chance (suggesting some kind of coordination or redundancy between
them), or significantly *fewer* (suggesting some kind of mutual exclusivity)?

**The test, in plain language.** Picture each regulator's set of target genes as a
handful of colored balls pulled from a big bag of ~3,026 total genes — this is exactly
the "urn" analogy from the hypergeometric test definition in Part 4. If you know exactly
how many targets each regulator really has, the hypergeometric test tells you precisely
how surprising *any* particular amount of overlap between two regulators' target sets
would be, purely from chance, given those sizes. For three or four regulators at once,
this project uses a mathematically exact extension of the same idea (explained in full
technical detail in the `enrichment.py` source file for readers who want the derivation)
rather than an approximation or a computer simulation — every reported p-value in this
analysis is an exact calculation, not an estimate.

**The scale of the test.** With 11 regulators, there are 55 possible pairs, 165
possible groups of three, and 330 possible groups of four — 550 separate tests in total.
Because running that many tests at once creates a real risk of false alarms purely from
volume (explained in Part 4's "multiple testing" entry), every result is corrected with
the Benjamini-Hochberg procedure before being called "significant."

**The result.**

| Group size | Number of possible groups | Number statistically significant (q <= 0.05) |
|---|---|---|
| Pairs (2 regulators) | 55 | 41 |
| Triples (3 regulators) | 165 | 19 |
| Quartets (4 regulators) | 330 | 0 |
| **Total** | **550** | **60** |

**Of all 60 significant results, every single one is a depletion — regulators sharing
*fewer* target genes than chance predicts. Not one is an enrichment.** For example, the
regulators NCU06108 and lhp-1 share only 17 target genes in common, where roughly 141
would be expected purely by chance given how many targets each one has individually —
an outcome so far below chance that it would happen by pure luck roughly 1 time in
10^50 (quantified precisely by the exact hypergeometric calculation, not estimated).

**An even stronger, non-statistical finding for triples and quartets.** Beyond the 60
statistically significant results, a direct check of the raw data reveals something more
absolute: **every single one of the 165 triples and every single one of the 330
quartets — 495 groups in total, with zero exceptions — has an observed overlap of
*exactly* zero target genes.** This is not merely "statistically unlikely" the way the
60 significant pairs/triples above are — it is a mathematical certainty, following
directly from the single-winner-per-gene structural fact explained in Part 5.2. If every
target gene, by construction, has at most 2 controlling regulators recorded against it,
then no gene can *ever* appear in the shared-target list of 3 or more regulators
simultaneously — there is no gene left over to be "shared" by a third or fourth
regulator once the first two have claimed it. The 330 quartet tests exist specifically
to make this structural ceiling explicit and to quantify, with real numbers, exactly how
far short of "chance-level" sharing this zero-overlap outcome falls (the quartet closest
to reaching statistical significance, {pro-1, NCU06108, NCU07155, lhp-1}, still only
reaches a p-value of about 0.12 — nowhere near the conventional 0.05 threshold, let alone
after the Benjamini-Hochberg correction for having run 550 tests at once — because
roughly 2.8 shared target genes would already have been expected by chance for that
particular foursome, and even that modest expectation isn't a large enough number to
make "we observed zero instead" look truly extraordinary).

It's worth being precise about what does and doesn't hold at the pairwise level, since it's
the one place this "always zero" pattern breaks: 45 of the 55 regulator pairs *do* have
at least some real overlap in their target-gene sets (it's specifically once you require
3 or more regulators simultaneously that the export's structure makes overlap
impossible). So the honest summary is: **pairs of regulators sometimes do, and sometimes
don't, share targets — and when they do share fewer than expected, that's a real,
statistically supported pattern in this specific exported data. But three or more
regulators sharing any target at all is not just rare — in this dataset, it cannot
happen, by construction, full stop.**

As with the FFL result in Part 5.3, none of this should be read as a biological claim
about real antagonism or competition between the actual regulator proteins in a living
cell — it is a precise, correctly-calculated statistical description of the *exported,
thresholded* connection matrix specifically, and its interpretation is inseparable from
the winner-take-all construction explained in Part 5.2.

## 5.5 What do the continuous binding strengths actually look like?

Covered in data terms in Part 3.2. In dashboard terms: this project renders the 21
functional-category x 11-regulator binding-strength table as a heatmap and a full data
table, specifically so a reader can see a real, continuous slice of the ensemble's
output — not just the collapsed 1/0 winner-only view everything else in this project
necessarily works from. The accompanying explanation is explicit that a blank/zero cell
here means "not the dominant regulator for this function" and not "confirmed to have
zero effect," for exactly the reasons discussed in Parts 2.2 and 3.2.

## 5.6 The Mu/Theta simulated dataset: one mystery solved, one still open

Covered in data terms in Part 3.3. To summarize the investigative process and its
outcome:

- **Gene identity: confirmed.** Every one of the 2,216 distinct genes referenced in
  `ALLNCU.txt` is a real target gene in this project's own cleaned connection matrix —
  100% overlap, not a coincidental resemblance. This dataset genuinely shares its gene
  universe with the rest of the project.
- **What "Mu" represents: now confirmed**, directly from the data source (see Part
  3.3): a conditional ensemble average — the average binding strength of a regulator,
  computed only over the subset of the 2,000 MCMC ensemble runs in which that regulator
  happened to be each gene's dominant regulator. Still open: exactly which regulator
  each of the 10 Mu columns corresponds to, positionally.
- **What "Theta" represents: tested and refuted, for the one hypothesis tried.** Because
  Theta is a 6-value probability distribution per gene, and this project's data has
  exactly 6 post-transcriptional ("RNA operon") regulators (Part 1.4), the natural first
  hypothesis was that Theta encodes each gene's soft assignment across those 6
  regulators. This was tested rigorously: every one of the 6! = 720 possible ways to
  match Theta's 6 columns to the 6 real regulators was tried, scoring each by how often
  the best-matching regulator (by Theta's largest value) actually agreed with the gene's
  real recorded post-transcriptional regulator, across 887 genes with one unambiguous
  real regulator to check against.

| Quantity | Value |
|---|---|
| Genes evaluated | 887 |
| Best possible accuracy, out of all 720 label matchings tried | 28.7% |
| Accuracy from simply always guessing the most common real regulator (lhp-1) | 56.4% |
| Average accuracy across all 720 possible label matchings (chance level) | 16.7% |

**The best possible matching (28.7% accuracy) is worse than the trivial strategy of
always guessing the single most common real regulator (56.4% accuracy).** In other
words, even given the best possible way of lining up Theta's 6 columns with the 6 real
regulators, Theta's largest value is a *worse* predictor of a gene's real
post-transcriptional regulator than simply ignoring Theta entirely and always guessing
"lhp-1." This decisively rules out the tested hypothesis. What Theta actually represents
remains an open question — flagged honestly as unresolved rather than forced into a
convenient but unsupported story.

## 5.7 Do the clock regulators, or their targets, have evolutionary "siblings"?

This is the one analysis in the project based on actual gene *sequence* data (Part 3.4)
rather than activity/expression data, and it asks a question none of the rest of the
project, and not the original paper either, can answer: **do any of the 11 clock
regulators — or any of their ~3,000 target genes — have close evolutionary relatives
(paralogs, defined in Part 4) elsewhere in the genome?** This matters because a small
family of duplicated, closely related genes can sometimes back each other up
functionally, or represent an evolutionary elaboration of a single ancestral regulatory
role — information invisible to any analysis of activity data alone.

**Method, briefly.** With no specialized bioinformatics search software available in
this environment, this project built a two-stage approach instead: first, a fast,
approximate pre-filter (comparing short 5-letter fragments of each protein sequence to
narrow ~10,255 x 10,255 possible gene pairs down to a manageable few hundred
candidates), and second, real, rigorous sequence alignment (the same general scoring
approach used by industry-standard tools like BLAST) to confirm or reject each
candidate pair. An extra complication surfaced early and was handled explicitly: some
apparent "duplicate pairs" turned out to just be the same physical gene filed under two
different names in the public database (identical sequences, different labels) rather
than true evolutionary duplicates — these were detected and merged into one gene before
the real family-finding step, so reported family sizes reflect genuine duplication, not
database bookkeeping.

**Regulator ID resolution.** Before any sequence comparison could happen, each of the
11 regulators' short paper-given names had to be matched to a real entry in the public
protein database. 8 of 11 were confidently resolved; 3 were not (one likely a naming
collision with an unrelated, well-established gene; one could not be located under any
name variant tried, possibly renamed in a newer genome annotation than the paper used;
WCC itself is a two-protein complex rather than a single gene product, and was excluded
from this specific analysis for that reason). The 3 unresolved regulators were excluded
from the paralog search entirely, rather than guessed at.

**Result.** None of the 8 successfully resolved clock regulators turned out to have any
detected paralog at all — each one is evolutionarily a "singleton" as far as this
analysis could determine. Among the ~3,026 target genes, however, real paralog families
were found, including a family of 3 ADP-ribosylation-factor genes (one of which is
itself a clock-network target gene) and a family of 4 xylanase/xylan-esterase genes
(one of which is also a clock-network target). The full family-size breakdown across
all 9,747 genes considered (after merging duplicate database records) is: 9,538 genes
with no detected family (singletons), 72 families of size 2, 6 of size 3, 5 of size 4,
1 of size 5, 1 of size 6, and 2 of size 8.

**Honest limitations.** This method can miss true but highly divergent (very distantly
related) paralogs entirely, since it has no more sophisticated profile-matching step;
its clustering approach can occasionally over-merge genes that share one narrow
functional feature without being true close relatives; and unlike a real BLAST search,
there is no formal statistical-significance score behind the specific similarity
thresholds used. Every family reported here should be read as a reasonable, biologically
plausible first-pass hypothesis, not a publication-grade evolutionary claim.

## 5.8 Seeing all of it: the interactive dashboard

Every analysis described above is also presented visually in the accompanying Streamlit
dashboard (`app/streamlit_app.py`), across 9 tabs, each with its own in-context
explanation written at the same non-expert-friendly level as this document:

- **Overview** — a short orientation and walkthrough for a first-time visitor.
- **Matrix Explorer** — the raw and cleaned connection matrix itself, browsable.
- **Binding Strengths** — the continuous function-level estimates from Part 5.5.
- **Mu/Theta Sim** — the investigation from Part 5.6, including a per-gene lookup tool.
- **FFLs & Null Model** — the feedforward-loop analysis and randomization test from
  Part 5.3.
- **Enrichment Tests** — the pairwise/triple/quartet co-targeting analysis from Part
  5.4, with filters for group size, direction, and significance.
- **Network Viz** — interactive diagrams: a hub-and-spoke view centered on WCC, a full
  force-directed graph of every regulator and target gene, and a dedicated feedforward-
  loop diagram, in both 2D and 3D.
- **Gene Families** — the paralog analysis from Part 5.7, browsable and searchable.
- **Methodology & Paper Comparison** — the validation table from Part 5.2 and a
  file-by-file map of what's reproduced from the original paper versus original to this
  project.

---

# Part 6: Honest Limitations, Stated Directly

In keeping with this project's guiding principle — report discrepancies and open
questions directly, rather than quietly smoothing them over — here is a consolidated
list of everything in this project that remains uncertain or incomplete:

- The exported connection matrix caps every gene at a single dominant regulator (plus
  WCC), which — as explained at length in Part 5.2 — is very likely an artifact of how
  the data was exported for this project, not a claim that the original paper's own full
  ensemble result only ever supports one controlling regulator per gene. Several
  downstream findings (the FFL null-model result, the enrichment depletion pattern, and
  especially the "triples/quartets are always exactly zero" finding) should be
  understood as properties of *this specific exported snapshot*, not necessarily as
  properties of the underlying biology or of the original paper's full statistical
  result.
- The Mu dataset's column-to-regulator mapping (which of the 10 Mu columns is which
  regulator) remains unconfirmed.
- What the Theta dataset represents remains genuinely unknown — the one hypothesis
  tested for it was clearly rejected, and no alternative hypothesis has yet been tested.
- The paralog/gene-family analysis (Part 5.7) uses an approximate, from-scratch method
  rather than industry-standard search software, and can both miss distant true
  relatives and occasionally over-group unrelated genes that happen to share one
  feature.
- 3 of the 11 regulators could not be confidently matched to a protein sequence record
  at all, and were excluded from the paralog analysis rather than guessed at.

---

# Part 7: Citation

If referencing the original data and model this project is built on, please cite:

> Al-Omari, A. M., Griffith, J., Scruse, A., Robinson, R. W., Schüttler, H.-B., &
> Arnold, J. (2022). Ensemble Methods for Identifying RNA Operons and Regulons in the
> Clock Network of *Neurospora Crassa*. *IEEE Access*, 10, 32510-32524.
> https://doi.org/10.1109/ACCESS.2022.3160481

All statistical extensions described in Parts 5.2-5.7 of this document (the null-model
significance test for feedforward loops, the pairwise/triple/quartet enrichment tests,
and the sequence-based paralog family analysis) are original to this project and were
not performed in the cited paper.
