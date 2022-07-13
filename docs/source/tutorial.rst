Tutorial
========

RNAmediator_plfold
################

Using Genes
**********************

ConstraintPLFold.py can be used for folding a list of sequences with constraints.
Therefore it is recommended to use one fasta file containing sequences and a corresponding bed file containing the constraints. Matching between sequences and constraints is done via the gene identifier (here **ENSG00000270742**) which has to be the same in the bed file (Column 4)and the fasta header and should look as follows:

**FASTA** ::

    >ENSG00000270742:chr1:61124732-61125202(+)
    TTTTTTCTTTATAATTATTCCCCTATTTGAAAAATCAACTTGTATATGAGGCAGCAAACACCTTGCAGAGC...

**BED** ::

    chr1	26	35	ENSG00000270742 .	+


A FASTA File with header following this format can easily be generated using `getfasta` from the `BedTools`_ suite on a gene coordinates BED file with options `-name` and `-full-header` enabled.

.. _BedTools: https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html

The command line call for this simple scenario is supposed to look like this, where `window`, `span` and `ulim` are options for RNAplfold_ and will be used to collect and summarize results in the the `RNAmediator_collect_plfold` step:

.. _RNAplfold: https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html

.. code-block:: shell

    RNAmediator_plfold -s FASTA -x BedFile --window SizeOfWindowToFold --span BasepairSpan -u ulim


By default all results are printed to stdout. The first output is the folded sequence without constraints. This output can be redirected to a file using ``--unconstraint nameforunconstraint``. The following outputs are the matrices for an unpaired constraint and for an paired constraint. These matrices can also be redirected to a file using ``--paired nameforpaired`` or ``--unpaired nameforunpaired``. The folder to which these files should be saved is determined via the ``--outdir`` flag.

Constraints can also be passed using the ``--ono`` (one on one) flag like this:

.. code-block:: shell

    python ConstraintPLFold.py -s FASTA -x ono,BedFile

Matching between constraints and sequences consequently ignores the gene identifiers and instead the first sequence will use the first line in the bed file as a constraint.

If your constraints are specified using genomic coordinates, it is essential to provide another bed file with genomic coordinates.
This file can be passed using the ``--genes GeneBedFile`` and should look like:

**GeneBED** ::

    chr1  61124732  61125202  ENSG00000270742 .	+


Folding without constraint
**************************

In case you simply want to fold sequences without applying any constraints you can use the `scanning` mode of `RNAmediator_plfold` which will only create ``--unconstraint nameforunconstraint`` ouput as follows:

.. code-block:: shell

    RNAmediator_plfold -s FASTA -x scanning --window SizeOfWindowToFold --span BasepairSpan -u ulim


Folding with sliding window constraint
**************************************

In case you want to fold sequences without specific constraint positions you can use the `sliding` mode of `RNAmediator_plfold` (DEFAULT) which will generate a list of constraints covering every nucleotide of the input sequence and fold each independently. This will take long and produce large amounts of output, so make sure this is really what you want. It does, however, provide one with the opportunity to scan for regions with largest impact on structure along a given sequence, e.g. to identify potential miRNA or antagonist binding sites:

.. code-block:: shell

    RNAmediator_plfold -s FASTA -x sliding --window SizeOfWindowToFold --span BasepairSpan -u ulim


Using Transcripts or single gene features
******************************************

If you plan to use transcript isoforms instead of genes it is highly recommended to change the header of your fasta sequences as well as providing local instead of genomic coordinates for constraints in the bed file.

If using single exon or intron sequences, make sure they have unique identifiers (e.g. ENST000000001.5_intron1, ENST000000001.5_intron2 etc.) in both files, using just the transcript ID is not enough to map constraints to the respective sequences in that case.

 For transcripts this can look as follows:

**FASTA** ::

    >ENST00000240304.5::ENST00000240304.5:0-5482(.)
    GUCUUGUCGGCUCCUGUGUGUAGGAGGGAUUUCGGCCUGAGAGCGGG...

**BED** ::

    ENST00000240304.5	1178	1183	ENST00000240304.5	.	.

In this case you should *ALWAYS* provide a 'local' GeneBed file (``--genes``) which looks quite similar to the Constraints File but genes start at position 0 and end at the length of the sequence. It is not important to use this file in the Rissmed_plfold call.

However, it is essential for the `RNAmediator_collect_plfold` step to generate valid genomic coordinate BED files as output.

**GeneBED** ::

    ENST00000240304.5	0	5482	ENST00000240304.5	.	.


RNAmediator_collect_plfold
######################

The methods mentioned in the `RNAmediator_plfold` example will
produce output that can be processed by `CollectConsResults.py`. This will generate BED files storing the probability of being unpaired for spans of nucleotides around (not overlapping) the constraint. Make sure to provide a comma separated pattern for window and span as used in the `RNAmediator_plfold` call (here SizeOfWindowToFold,BasepairSpan) and ulim either similar or lower than in the `RNAmediator_plfold` call. Therefore simply call:

.. code-block::

    Rissmed_collect_plfold -d path/to/ConstraintPLFold/output -g GeneBed --outdir path/to/outdir --pattern SizeOfWindowToFold,BasepairSpan -u ulim

.. note::

    The collection step will only work if you provide coordinates of the gene or transcript via a **GeneBED** file.

The ``-u`` parameter defines the span sizes that are used for the output BED files which might
for example look like this for ``-u 5``:

::

    chr1	110135760	110135765	ENSG00000065135|44542-44544|110135774-110135777	0.02620831	+	9	0.05007846	0.07628677	2.2444807817789587	0.02620831000000001	0.04615814895706037	0.15227569
    chr1	110135761	110135766	ENSG00000065135|44542-44544|110135774-110135777	0.019514269999999993	+	8	0.049938	0.06945227	2.426255722126518	0.019514269999999993	-0.03666320494171017	0.15227569
    chr1	110135762	110135767	ENSG00000065135|44542-44544|110135774-110135777	0.09578122	+	7	0.13804191	0.23382313	1.4457214540425338	0.09578122	0.9069421599395611	0.15227569
    chr1	110135763	110135768	ENSG00000065135|44542-44544|110135774-110135777	0.07929997	+	6	0.06154513	0.1408451	1.5621026141257632	0.07929996999999998	0.7030295094408919	0.15227569
    chr1	110135778	110135783	ENSG00000065135|44542-44544|110135774-110135777	0.31023007	+	-7	0.14063798	0.45086805	0.7213795423617754	0.31023007	3.560189541047878	0.15227569
    chr1	110135779	110135784	ENSG00000065135|44542-44544|110135774-110135777	0.2991335	+	-8	0.1657577	0.4648912	0.7438289320089759	0.2991335	3.4228983161623856	0.15227569
    chr1	110135780	110135785	ENSG00000065135|44542-44544|110135774-110135777	0.29541885	+	-9	0.15833754	0.45375639	0.7515304743397738	0.29541885	3.376939173064934	0.15227569

.. note::

    If you used different plfold parameters (-w, -l) in the RNAmediator_plfold call, you have to adapt the pattern accordingly


Output
#######

`RNAmediator_plfold`
****************

Per default, `RNAmediator_plfold` prints to STDOUT, or dumps `numpy`_ arrays to disk following the naming provided by the user with the ``--unconstraint nameforunconstraint``, ``--paired nameforpaired`` and ``--unpaired nameforunpaired`` options. These `numpy`_ arrays are then used as input for `RNAmediator_collect_plfold`.
If the user prefers to also generate human readable files similar to `RNAplfold` output, the option `--save 1` has to be set. This will generate '.gz' files providing the same output as a commandline call to `RNAplfold`.

.. _numpy: https://numpy.org/doc/stable/reference/generated/numpy.array.html

The standard user will not want to work on this output directly but generate genomic coordinate BED files, summarizing the effect of ligand binding as will be explained next.

`RNAmediator_collect_plfold`
************************

The genomic coordinate BED file generated by `RNAmediator_collect_plfold` contains the following columns:

::

    Chr    Start   End     Constraint      Accessibility_difference        Strand  Distance_to_constraint  Accessibility_no_constraint     Accessibility_constraint        Energy_Difference       Kd_change       Zscore  Accessibility_constraint_pos

.. note::

    The first 6 fields follow standard BED format, with the Difference in Accessibility in Column 5 and the constraint as well as its local and genomic coordinates as Identifier in Column4.
    `Distance_to_constraint` shows the distance between the region that changed and the applied constraint
    `Accessibility_no_constraint` is the accessibility of this region before the constraint was applied
    `Accessibility_constraint` is the accessibility of this region after the constraint was applied
    `Energy_Difference` is the change in 'Free Energy' of the structure after applying the constraint
    `Kd_change` is the influence on the Kd of a potential binding partner after the constraint has been applied
    `Zscore` is the Zscore of accessibility changes in this position after the constraint has been applied in comparison to all changes in all positions
    `Accessibility_constraint_pos` shows the accessibility at the position that was constraint before the constraint was applied


Further steps
##############

The BED file created by CollectConsResults can be used to intersect with other known binding sites on the same gene/transcript. Thus, it is possible to see whether the changes in RNA structure upon binding of one ligand might affect the structure of binding site of other ligands.


RNAtweaks
##########


RNAplfold
**********

For RNAplfold usage two different wrappers exist. One uses the command line version of RNAplfold and the other
uses the ViennaRNA API

.. code-block:: python

    from RNAmediator.RNAtweaks import RNAtweaks
    sequence = "AAATTTTGGGGGGCCCC"
    window = 3  # winsize option of RNAplfold
    span = 3   # span option of RNAplfold
    region = 3  # ulength option of RNAplfold
    constraint = ("paired", 3, 5)
    api_result = RNAtweaks.api_rnaplfold(sequence, window, span, region=region, temperature=37.0, constraint=[constraint])
    cmd_result = RNAtweaks.cmd_rnaplfold(sequence, window, span, region=region, temperature=37.0, constraint=[constraint])

For now only paired and unpaired constraints are supported. The constraints must be a list of Tuples in the format ``("paired"/"unpaired", start, end)`` in contrast to the ViennaRNA constraints these are zero based.
Both calls will produce an identical PLFoldOutput object that reflects the ViennaRNA `_lunp` file.

PLFoldOutput
**************
Object that reflects the ViennaRNA `_lunp` file. The objects supports various functions to get different representations of the file. The recommended usage is to produce an numpy array as follows:

.. code-block:: python

    array = api_result.numpy_array


However, it is also possible to get the text representation of the file, which is usually produced by RNAplfold via:

.. code-block:: python

    array = api_result.get_text(nan="NA", truncated=False)

Hereby nan replaced the non float values with ``"NA"`` and the truncated flag is used to either keep or drop the header lines beginning with "#".





