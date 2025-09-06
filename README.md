***DepMine*** – v.GPEFS_2023Q2_6 – USER MANUAL

13/4/2024

**INSTALLATION**

***Dependencies***

***DepMine*** is a Python 3.x program running on Python versions 3.8 and
above. It requires the following Python libraries to be available :

**pandas, numpy, matplotlib, collections, scipy, timeit,
multiprocessing, re, requests, math, seaborn, random, glob**

Modules not already included in your local Python should be installed
with **pip** or **conda**.

***Parallel Processing***

DepMine makes extensive use of parallel processing on multiple CPUs, and
the speed of time-sensitive operations will be greatly improved by
running the program on a multi-CPU workstation.

***Data Files***

DepMine uses a number of datafiles which are either supplied with the
code or can be downloaded from a range of different sources :

| **File** | **Data** | **Web Page** | **Version** |
|----|----|----|----|
| CRISPRGeneDependency.csv | gene dependency | https://depmap.org/portal/download/all/ | Public 23Q2  |
| Model.csv | cell line information | https://depmap.org/portal/download/all/ | Public 23Q2  |
| OmicsSomaticMutations.csv | mutations | https://depmap.org/portal/download/all/ | Public 23Q2  |
| OmicsCNGene.csv | copy number | https://depmap.org/portal/download/all/ | Public 23Q2  |
| OmicsFusionFiltered.csv | gene fusions | https://depmap.org/portal/download/all/ | Public 23Q2  |
| OmicsExpressionProteinCodingGenesTPMLogp1.csv | expression | https://depmap.org/portal/download/all/ | Public 23Q2  |
| c1.all.v2023.1.Hs.json.txt | chromosomes | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/ | 2023.1.Hs |
| hgnc_complete_set.json | HUGO gene names | https://www.genenames.org/download/archive/ | Current |
| c2.cp.reactome.v2023.2.Hs.symbols.gmt | pathways | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/ | 2023.2.Hs |
| Homo_sapiens.GRCh38.110.refseq.tsv | ENSEMBLE reference sequences | https://www.ensembl.org/Homo_sapiens/Info/Index | GRCh38 |
| ChromosomeLengths.csv | Chromosome lengths | local |  |
| ExpressionThresholds.csv | Expression thresholds | local |  |
| Assignment_Solution_Activities.tsv | Mutation signature assignments | local |  |
| Final_GOF_all_table | Annotated GoFs | local |  |
| AB_Dep_Profiles | cancer profiles | local |  |

On first run, ***DepMine*** will compile two further local files :
AB_GeneDatabaseExpr2023Q2 and AB_CellDatabaseExpr2023Q2 which should be
kept in the same directory as all the other data files and the program
file.

**  
RUNNING**

***DepMine** can be run from a Linux/UNIX shell command line (after
making the distributed program file executable) :*

**chmod a+x DepMine.py**

**./DepMine.py**

***DepMine*** will print the number of CPUs it can detect and then give
the option of the user specifying the individual data files to be used.
Unless you have very good reasons not to, accept the default of ‘**N**’.
for this, and ***DepMine*** will load the default set of data files. If
this is the first run, two local databases will be generated – which
takes a few minutes - and some commentary will be given on the progress
of this.

As a general rule, when ***DepMine*** offers choices to the user, they
will be given in a bracketed list, and the default choice – i.e. what
will be used if you just tap return or type something that is not
offered – will be given in square braces. For example :

**'Aggregate by genes, residues, or not aggregate ? (G,R,\[N\]) : ')**

Again as a general rule, ***DepMine*** only needs the first character of
options to be specified and won’t care about upper or lower case.

Once the data files are loaded – usually within 60 seconds –
***DepMine*** will prompt for a command.

**NOMENCLATURE**

As a general principal, genes that are candidates for therapeutic
intervention on the basis of their dependency in different cell lines,
are referred to as GeneAs, while those that are altered in their
sequence, copy or expression level in cancer cell lines are referred to
as GeneBs. Synthetic sickness/lethality relationships are calculated
between GeneAs and GeneBs. GeneAs are also referred to as targets.**  
COMMANDS**

The available user commands are :

| **CHRO** | Show chromosomal location of a list of genes |
|----|----|
| **CYTO** | Show chromosomal location of the CNV amplifications and deletions in a cell line |
| **DEFI** | Define a *cancer profile* |
| **EXPR** | Analyse expression profile of a gene |
| **FIND** | Find genes whose disruptive mutation and/or deletion and/or low expression switches the dependency of a target gene |
| **GOF** | Find genes whose dependency is switched by the presence of one or more documented gain-of-function mutations |
| **LIST** | Apply **FIND** to a list of genes |
| **LOCU** | Scan genes in a cytogenetic locus to find those where mutational disruption switches the dependency of a target gene |
| **MINE** | Find target genes whose dependency is switched by a defined cancer profile, and optionally refine the *cancer profile* used |
| **PARA** | Set various parameters and thresholds |
| **PROF** | Read/write/organise *cancer profiles* and *cancer profile* collections |
| **QUIT** | Does what it says |
| **SIG** | Analyse the mutational signatures in cell lines matching a *cancer profile* |
| **SWITCH** | Find gene whose dependency is switched in cell lines matching a *cancer profile* |
| **TYPE** | Analyses the tissue type distribution of cell lines matching a *cancer profile* |

The individual commands are described in detail below.

**CHRO** – *Chromosomal Gene Location*

You will be prompted for a number of additional files - just press
return – and then for a filename. This file will typically be generated
by the **FIND** command – see below – but minimally needs to be a CSV
file with a column called ‘GeneB’. This command needs to interrogate
ENSEMBLE – so you need a good internet connection to utilise it. The
command generates a graphic of human chromosomes, with the location of
the GeneB genes indicated.

**CYTO** – *Chromosomal CNV Location*

You will be prompted for a number of additional files - just press
return – and then for a the name of a cell line. This can be a DepMine
identifier e.g. ACH-002968 or a common cell line name e.g. RPE1 or U2OS.
If the name given is not unique a list of the matching cell lines is
given so a unique name can be given.

The command generates a graphic of human chromosomes with amplified
genes indicated in blue and deleted genes in red.

**DEFI** – *Cancer Profile Definition*

You will be prompted for a *genetic profile* definition*,* which can be
just the name of an existing *cancer profile* or an algebraic statement
utilising one or more primitive postulates constructed from the
following operators and combined with set operations : union **\|**,
intersection **&** and exclusion **^**, and **( )**.

<table>
<colgroup>
<col style="width: 10%" />
<col style="width: 27%" />
<col style="width: 42%" />
<col style="width: 19%" />
</colgroup>
<thead>
<tr>
<th><strong>Operator</strong></th>
<th><strong>Operand</strong></th>
<th><strong>Resultant set</strong></th>
<th><strong>Example</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>+</strong></td>
<td>gene</td>
<td>cells where gene is amplified <sup>1</sup></td>
<td><strong>+myc</strong></td>
</tr>
<tr>
<td><strong>+</strong></td>
<td>cytogenetic locus</td>
<td>cells where &gt;=1 gene in locus is amplified <sup>1</sup></td>
<td><strong>+17q12.4</strong></td>
</tr>
<tr>
<td><strong>−</strong></td>
<td>gene</td>
<td>cells where gene is deleted <sup>1</sup></td>
<td><strong>−tp53</strong></td>
</tr>
<tr>
<td><strong>−</strong></td>
<td>cytogenetic locus</td>
<td>cells where &gt;=1 gene in locus is deleted <sup>1</sup></td>
<td><strong>−9p21</strong></td>
</tr>
<tr>
<td><strong>=</strong></td>
<td>gene</td>
<td>cells where gene is diploid <sup>1</sup></td>
<td><strong>=snca</strong></td>
</tr>
<tr>
<td><strong>~</strong></td>
<td>gene</td>
<td>Cells where gene is involved in a gene fusion</td>
<td><strong>~alk</strong></td>
</tr>
<tr>
<td><strong>&gt;</strong></td>
<td>gene</td>
<td>cells where gene is over-expressed <sup>2</sup></td>
<td><strong>&gt;erbb2</strong></td>
</tr>
<tr>
<td><strong>&lt;</strong></td>
<td>gene</td>
<td>cells where gene is under-expressed <sup>2</sup></td>
<td><strong>&lt;tert</strong></td>
</tr>
<tr>
<td><strong>!</strong></td>
<td>gene</td>
<td>cells where gene has a reading-frame-disruptive loss-of-function
mutation</td>
<td><strong>!brca2</strong></td>
</tr>
<tr>
<td><strong>/</strong></td>
<td>gene:mut</td>
<td><p>cells where gene has defined missense mutation</p>
<p>‘?’ matches any change</p></td>
<td><p><strong>/braf:v600e</strong></p>
<p><strong>/kras:g12?</strong></p></td>
</tr>
<tr>
<td><strong>/</strong></td>
<td>gene:mut,mut[,mut]</td>
<td>cells where gene has one of several defined missense mutations</td>
<td><strong>/egfr:l834r,t766m</strong></td>
</tr>
<tr>
<td><strong>/</strong></td>
<td>gene+</td>
<td>Cells where gene has an annotated gain-of-function mutation</td>
<td><strong>/kras+</strong></td>
</tr>
<tr>
<td><strong>/</strong></td>
<td>gene<strong>−</strong></td>
<td>Cells where gene has an annotated loss-of-function mutation</td>
<td><strong>/atm−</strong></td>
</tr>
<tr>
<td><strong>$</strong></td>
<td>gene</td>
<td>cells where gene has no non-synonymous mutations</td>
<td><strong>$rb1</strong></td>
</tr>
<tr>
<td><strong>{</strong></td>
<td>Cosmic signature<sup>3</sup></td>
<td>Cells displaying a defined Cosmic mutational signature</td>
<td><strong>{sbs7a</strong></td>
</tr>
<tr>
<td><strong>{</strong></td>
<td>Cosmic signature aetiology</td>
<td>Cells displaying one or more Cosmic mutational signature attributed
to a defined aetiology</td>
<td><strong>{tobacco</strong></td>
</tr>
<tr>
<td><strong>#</strong></td>
<td></td>
<td>all cells (Universal Set)</td>
<td></td>
</tr>
<tr>
<td><strong>%</strong></td>
<td>tissue</td>
<td>cells derived from tissue</td>
<td><strong>%breast</strong></td>
</tr>
<tr>
<td><strong>*</strong></td>
<td><p>sex</p>
<p>(Male,Female,Unknown)</p></td>
<td>cells with sex of patient</td>
<td><strong>*fem</strong></td>
</tr>
<tr>
<td><strong>;</strong></td>
<td><p>age</p>
<p>(Adult, Paediatric,Unknown)</p></td>
<td>cells with age of patient</td>
<td><strong>;p</strong></td>
</tr>
<tr>
<td>@</td>
<td>Filename</td>
<td>externally sourced CSV file of cell ACH codes</td>
<td><strong>@CHAN-MSI</strong></td>
</tr>
<tr>
<td colspan="4"><ol type="1">
<li><p><em>Definitions of amplification, deletion and diploid status are
based on Mina et al 2020, adapted to the log<sub>2</sub>(CN ratio + 1)
formalism used by DepMap/CCLE.</em></p></li>
<li><p><em>Over- and under-expression are defined per gene, based on
prior analysis of DepMap expression level values across the full
population of cell lines for which data is available (see
below).</em></p></li>
<li><p><em>As defined in</em> <em>(Diaz-Gay et al, 2023).</em></p></li>
</ol></td>
</tr>
</tbody>
</table>

A profile name can be specified directly e.g. **PROFILE_NAME = \<profile
definition\>**, however if no profile name is provided one will be asked
for before the definition can be saved.

In response to the definition, *DepMine* will report the number of cell
lines that are matched by the *cancer profile*, and offer to save/list
the matching cell lines. Finally, you will be prompted to save the
profile by name. If the name is already in use, the choice of
overwriting the previous definition will be given. If a profile is
overwritten, any other profiles that reference it will still work, but
their behaviour will reflect the new definition of the overwritten
profile. If in doubt, do not overwrite existing profiles.

Defined profiles will only be saved to disc, by explicitly saving them
in the **PROF** module.

**EXPR** – *Expression Profiles*

You will be prompted for a gene name and histogram of the distribution
of expression levels for that gene across the entire cell line
collection is generated. Additionally, scatter plots of copy number
versus expression level and gene dependency versus expression level, are
generated.

**FIND** – *Find SSL Partners for a Specified Target Gene*

You will be prompted for a target GeneA name, and then whether mutation
and/or deletion and/or under-expression should be used to specify the
‘defect’ in GeneBs. This can be specified as 1, 2 or 3 of M, D and U, in
any order.

The command will then run a full genome calculation which can take up to
an hour, depending on whether 1, 2 or 3 defect types are being used and
the number of processors in your computer. Once completed, the sorted
list of ‘hits’ that pass the significance threshold is presented and for
each gene in order, you will get the option to plot comparative
dependency distributions for the cells where that gene is defective and
a control set of cells in which it isn’t.

**GOF** – *Find SSL Targets for Gain-of-Function Mutations*

You will be asked whether you want to use the GoF mutations annotated in
the DepMap data, or wish to read in an external GoF dataset. Usually the
DepMap annotations will be sufficient.

You will then be asked if you wish to :

- aggregate GoF mutations at the gene level – i.e. KRAS-G12C, KRAS-G13C,
  KRAS-Q61R will count as equivalent;

- aggregate GoF mutations at the residue level – i.e. KRAS-G12C,
  KRAS-G12V, KRAS-G12A will count as equivalent, but KRAS-G12? and
  KRAS-G13? will be considered separately.

- not aggregate – so all mutations are counted individually.

The command will then run a full genome calculation for each mutation or
aggregate looking for target GeneAs that show high dependency in cells
harbouring those mutations relative to cells in which these mutations do
not occur. For statistical power, only those mutations or aggregates
which occur in \> 5 cell lines are considered.

Depending on whether the command uses aggregated or non-aggregated
mutations and on the number of processors in your computer this can take
a couple of hours to complete. At the end of the run, a p-value sorted
list of ‘hits’ is written to disk for subsequent analysis.

A subordinate command **PGOF** reads back the **GOF** generated hit-list
and plots comparative dependency distributions for each target ‘hit’ in
cells harbouring the GoF mutations, and a control set of cells that
don’t.

**LIST** – *Batch SSL Finding*

You will be prompted for a CSV file that minimally contains pairs of
gene names (GeneA, GeneB) for which SSL will be calculated. An
additional ancillary column can be specified to be carried over into the
output, and the calculation can optionally be restricted to a specified
tissue type.

The command will generate an output CSV file with p-value and Cohen d
value for the comparison of dependency distributions for GeneA in cells
in which GeneB is disruptively mutated, CNV deleted and under-expressed.
The calculation is repeated with GeneA and GeneB swapped.

**LOCU** – *Identifying Direct Targets in Locus Deletions*

This command tries to identify individual GeneBs within a potentially
deleted cytogenetic locus that show a direct SSL with a target GeneA.

The command will prompt for a target gene name and a cytogenetic locus,
and will respond with the names of the genes at that locus. It will then
calculate SSL between the target gene and each gene in the locus but
only in cells where that gene is CNV diploid but has frame-disruptive or
missense LoF mutations, and report the results as a waterfall plot based
on Cohen d-value, marked up to indicate those genes where the SSL with
the target also has a p-value better than the threshold parameter.

**MINE** – *Optimise Cancer Profiles*

The command will prompt for a target gene name and a defined *cancer
profile*, and will then plot comparative dependency distributions for
the target gene in those cells that match the *cancer profile* and a
control set of cells that don’t.

You will then be offered the option to mine the cells matching the
*cancer profile* for factors additional to the original profile that
display significant over/under-representation in the high/low dependency
ends of that distribution. These are in order : sex (*female ,male,
unknown*), age (*paediatric, adult, unknown*), tissue, mutational
signatures, frame-disruptive / missense LoF mutations, GoF mutations,
copy number deletion/amplification (singleton gene and locus),
high/low-expression level (aggregated to pathways).

At each stage the distribution of identified factors can be plotted
overlayed on the parental matching distribution. Once all secondary
factors have been identified these can be reviewed and then amended
cancer profiles embodying all possible combinations of secondary factors
added to the original profile, are evaluated in terms of effect size and
the number of cell lines matched. The Pareto optimal surface for these
two parameters is plotted, and the refined profiles providing optimal
solutions are optionally written to a profile collection, which can be
loaded for subsequent use using **PROF**.

**PARA** – *Specify Parameters*

Allows you to override default values for various parameters :
dependency ‘high’ threshold,

minimal mutational frequency, maximum p-value for significance, minimum
Cohen d-value for a strong effect.

**PROF** – *Cancer Profile Management*

In this and all modules that utilise *cancer profiles*, if no *cancer
profiles* have been loaded previously you will be prompted to load a
cancer profile collection from a list. *DepMine* will build this list
from files in the launch directory of the type **filename.profiles** if
any are present. The default profile collection provided with *DepMine*
is called **AB_Dep_Profiles** and comes initialised with a single
profile – **ALL** – which encompasses all the cell lines in the DepMap
dataset for which any data is present.

Once a profile collection is loaded, the following subcommands are
available :

| **Show** | Lists the top-level definition of a profile including references to other profiles |
|----|----|
| **Delete** | Removes a profile from the loaded collection – if this is referenced by other profiles this will be require confirmation that all affected profiles will be removed. NB this only affects the loaded copy of the profiles |
| **Name** | Allows the name of a profile to be changed. If this is referenced by other profiles the references will be automatically updated to the new profile name |
| **Flatten** | Lists the definition of a profile with all references to other profiles flattened to primitive set postulate functions |
| **Read** | Reads in profiles from a file – optionally these are added to the profiles already loaded, or replaces them. NB any changes to the previously loaded profiles will not be automatically saved prior to replacement by the new set |
| **Write** | Writes the current loaded profiles to a file. |
| **Quit** | Quits this module |

**SIG** – *Mutational Signature Distribution*

You will be prompted for a named *genetic profile* and a histogram of
the distribution of mutational signatures within the cell lines matching
that profile will be generated.

**SWITCH** – *Identify Target Genes for Cancer Profiles*

You will be prompted for a named *cancer profile* and the command will
then perform a full genome calculation for looking for target GeneAs
that show high dependency in cells matching the *cancer profile*
relative to cells that do not match the profile. Depending on the cancer
profile and the number of processors in your computer the analysis will
take 10-20 minutes. On completion a waterfall plot of GeneA ‘hits’ based
on p-value is generated.

You will then be offered the opportunity to plot comparative dependency
distributions for the target genes in those cells that match the *cancer
profile* and a control set of cells that don’t.

The list of GeneAs that pass the requirement for a strong effect with
statistical significance can be saved as a CSV file.

**TYPE** – *Tissue Type Distribution*

You will be prompted for a named *genetic profile* and histograms of the
absolute and relative tissue distribution of the matching cell lines
will be generated

**APPENDIX – Generating Local Files**

ExpressionThresholds.csv can be regenerated from
OmicsExpressionProteinCodingGenesTPMLogp1.csv with the local program
***expr_scorer2.py***

Assignment_Solution_Activities.tsv can be regenerated from cell line vcf
files by https://cancer.sanger.ac.uk/signatures/assignment/

The directory of vcf files required can be generated from
OmicsSomaticMutations.csv with the local program ***MakeVCF.py***
