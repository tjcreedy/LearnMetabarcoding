===================================================
Filtering Amplicon Data
===================================================

------------------------------------------
Introduction
------------------------------------------

In this section we will cover a range of ways to filter sequence reads. It is important to note that for instructional purposes, we are presenting each filtering type separately, with separate commands. Experienced amplicon bioinformaticians will soon notice that this is not the most efficient pipeline, but our primary objective here is clarity in the communication of these concepts in an order that is easy and logical to understand.

------------------------------------------
Quality Filtering
------------------------------------------

Quality filtering aims to remove sequences that contain sequencing errors.

We have seen that the strings of letters in fastq files stand for quality values, but what do these values actually mean? Phred scores are not just some arbitrary values. They are actually logarithmically related to the probability of a base calling error. Given a phred score S, the probability of the base being incorrect P is:

.. math:: 

	P=10^{^{-S}/_{10}}

So a phred score of S=10 gives P=0.1, meaning there is a 10% chance of error. A score of S=20 means a 1% chance of error (P=0.01), and S=30 means a 0.1% chance of error (P=0.001).

Given this knowledge, we can use the Phred scores to calculate the number of expected errors over a sequence, by converting each base’s quality score into a probability and then summing the probabilities.

Given the below fastq, convert the phred scores to probabilities (look up the S values for each character on the fastq wikipedia page), then calculate the number of expected errors.

.. code-block:: rst
	
	@sequence1 
	AGTACTGCGC 
	+
	@EB>9:7/.(&

As you can imagine, there are many ways that Phred scores can be used to filter sequences. The three common ones are:

* Excluding reads with any scores below a threshold

* Trimming ends from reads based on a threshold

* Excluding reads based on number or rate of expected errors 

**footnote in PDF**

All three are commonly used in amplicon filtering.

With paired amplicons, we know that the end of the read is likely to be lower quality; but we have already merged our pairs. The quality scores in the middle of our reads have been adjusted to take account of the duplicate data, so in general it’s unlikely we have any specific poor regions of our reads. Excluding reads entirely based on quality is a more conservative approach and is generally suggested here.

For basic fastq filtering based on minimum score, we will use fastq_quality_filter from the fastx_toolkit package. This is a great little package of handy tools available here `http://hannonlab.cshl.edu/fastx_toolkit/index.html <http://hannonlab.cshl.edu/fastx_toolkit/index.html>`_ 

We will use the fastx_filter function from the VSEARCH software for filtering by expected error and expected error rate. VSEARCH is a software package specifically designed for metabarcoding, based on the USEARCH package but completely free and open source. We’ll see other tools from this useful package later; you can read the documentation here: https://github.com/torognes/vsearch

Let’s try doing some different sorts of quality filtering. Here are some commands.

Exclude any reads with quality score lower than 13 (probability of error ≈ 0.05):

.. code-block:: bash 

	$ fastq_quality_filter -q 13 -p 100 -i ​in.fastq​ -o ​out.fastq

Exclude any reads with fewer than 60% of bases with a quality score equal to or greater than 30 (p =
0.001):

.. code-block:: bash 
	
	$ fastq_quality_filter -q 30 -p 60 -i ​in.fastq​ -o ​out.fastq

Exclude any reads with more than 1 expected errors over the entire read:


.. code-block:: bash 

	$ vsearch --fastx_filter ​in.fastq​ --fastq_maxee 1 --fastaout ​out.fasta

Exclude any reads with more than 0.1 expected errors per base: 

.. code-block:: bash 

	$ vsearch --fastx_filter ​in.fastq​ --fastq_maxee_rate 0.1 --fastaout ​out.fasta

* Compare the number of reads returned with different filters using ``​grep​``. Note that ``fastq_quality_filter`` returns fastq files, but vsearch returns fastas. The regex for fastas should be "^>".

* Can you adjust the filters to get roughly the same number of reads filtered out using the different methods?

* Do you think these are the same reads that are being filtered each time?

Which quality filtering parameter to pick? Well, it depends partly on the nature of the data, partly on the aim of your filtering, and partly on what feels right to you.

*In my opinion, filtering based on the number of expected errors makes sense: there is a logical basis for the selection of the threshold, the removal of reads based on their overall likelihood of error, not some relatively arbitrary threshold of minimum quality score. While obviously this value increases with the length of the read, so could be argued isn’t a comparable value between different fragment lengths, my argument would be that it’s a reflection of the reality of sequencing, and that no matter how long my fragment is, I don’t want any errors. So I generally filter using* ``​--fastq_maxee 1`` ​. *If I suspect from later examination that I still have a lot of sequencing errors, I’ll reduce this to* ``--fastq_maxee 0.5`` . *If in the very rare case I’m simply not getting enough sequences returned, I might loosen this to* ``-​-fastq_maxee 1.5`` *or even* ``2`` *,​ but generally this data isn’t really trustworthy.* 

* Keep whichever one of your filtered fastas you like best. Delete the rest. This file will be the file used for the next step.

If you pick the output of ``fastq_quality_filter`` , you will need to convert this output to fasta using:

.. code-block:: bash 

	$ fastq_to_fasta -i ​in.fastq​ -o ​out.fasta

------------------------------------------
Not Filtering: Dereplication
------------------------------------------

Now that we’re down to just sequences, we can compress our dataset somewhat. This is because we should have many duplicate reads in this concatenated file due to PCR. It is important to retain our current file, because it keeps track of which reads came from which sequence, but now we are gradually preparing for finding the OTUs, which is a whole-dataset problem.

The subsequent steps are much faster if we remove these duplicates; we will record in the file headers how many copies there are of each sequence for later use.

VSEARCH has a dereplication function, which we will use here:

.. code-block:: bash 

	$ vsearch --derep_fulllength ​in.fasta​ --output ​out.fasta​ --sizeout --relabel uniq

* Use ``head`` to see how the sequence headers are formatted.

* Use ``grep`` to count up how many unique sequences we have.

We renamed these because the original names became meaningless once we dereplicated, so we might as well save space. Very importantly, we’ve included a “size” annotation that specifies how many copies each sequence had in the original dataset.

Make sure to keep the input file with all unique sequences! The dereplicated file has lost all of the information linking sequences to samples, but this remains in the input file, and we will use this much later. **For now we will use the dereplicated file for the next step​.** 

------------------------------------------
Indel Filtering
------------------------------------------

In this case, I use the term indel to mean insertions or deletions from PCR or sequencing errors rather than natural mutations. The fundamental assumption here is that the sequenced region is sufficiently conserved that there will not be any naturally-occuring indels sequenced, because these would have been deleterious and the organism would not have survived to have been sampled.

Thus, you should think carefully about how to apply this type of filtering to your own data depending on the barcode region used. Insertions or deletions are easy to spot because they will change the length of the sequence from what is expected based on the primers. While filtering based on length primarily removes indels, it can also be used to remove other reads that are clearly erroneous for other reasons.


Before we start, let’s double-check the length distribution of our reads. We can do this using a command we used before, having adapted the command for fastas (where the sequences are every other line, rather than every 4 lines):

.. code-block:: bash 

	$ sed -n '2~2p' ​file​ | while read l; do echo ${#l} ; done | sort | uniq -c

Oh dear, what’s happened to our reads? Check the first 10 lines of the fasta:

.. code-block:: bash 

	$ head -n 10 ​file

VSEARCH, although it’s great in many respects, outputs files in wrapped format, which means it starts a new line after 80 sequence characters. While this is nicer to look at, this is a pain for using quick-and-easy tools to summarise data on the linux command line. So we must run a quick command first to unwrap this data:

.. code-block:: bash 

	$ perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ​infile​ > ​outfile

Use the output from this in the sed command above.

If we have a very variable region, we might not want to do any filtering at all, or we may know a reasonable range of lengths within which we expect our reads to fall. There are lots of tools for length filtering; we’ll use VSEARCH again - in fact, the ``​--fastx_filter`` command again. Let’s try filtering with quite a wide range:

.. code-block:: bash 

	$ vsearch --fastx_filter ​in.fasta --fastq_minlen 400 --fastq_maxlen 440 --fastaout out.fasta​

In this case, the region of CO1 we use is sufficiently conserved that, on balance of probabilities, any insertions or deletions are due to PCR and/or sequencing errors, and/or maybe errors with the pair merging we did, rather than natural mutations. So, if you have a strict length expectation for your reads, you can exclude any reads longer or shorter than this value.

* Run the filtering again, this time allowing no variants from our target length of 418bp

What if we have a coding region like CO1 where we expect no single-base indels, but overall the region is more variable (or we have a wider range of taxa) and we might see real established whole-codon insertions or deletions?

* Briefly think about how we might specify a filter for this - assuming the same target length of 418bp, what lengths might we allow?

As it happens, it doesn’t seem that there are any programs out there that do this filtering already. One way to do it is filter by each length separately, and then concatenate the results - so, allowing one codon of variation, we would do:

.. code-block:: bash 

	$ for l in 415 418 421; do \
	> vsearch --fastx_filter mbc_concat.fastq --fastq_minlen $l \
	> --fastq_maxlen $l --fastaout mbc_concat_${l}_ctrim.fasta; \ 
	> done && cat *ctrim.fasta

We would want to make sure that worked properly by checking the number of sequences in the relevant files:

.. code-block:: 

	$ grep -c "^>" *ctrim.fasta

------------------------------------------
Frequency filtering / denoising
------------------------------------------

Due to the use of PCR to generate metabarcoding datasets, we should generate plenty of copies of the amplicon of interest, and thus many reads - depending on sequencing depth, of course. We can therefore be relatively confident that sequences that occur at very low frequencies in the dataset are more likely to be errors.

These can be filtered out with a simple threshold, although the selection of this threshold is likely to be dataset-dependent.

We can actually do this as part of the dereplication command we just did. Try running the command from the last section again, with the same input and adding the argument ``minuniquesize 1`` ​to get rid of all singleton sequences. **​Use a different output filename​**, we’re not going to keep this output.

A much more sophisticated approach to filtering errors is denoising. Denoising algorithms use read frequency and sequence composition to infer likely sequencing errors. Instead of doing the size filtering as part of dereplication, we will instead do it as part of a denoising command. We will use the unoise3 algorithm, implemented again in VSEARCH. Your input file here should be the output from dereplicating in the last section.

.. code-block:: bash 

	$ vsearch --cluster_unoise ​in.fasta​ --minsize 4 --unoise_alpha 2 --centroids out.fasta

**footnote** 

The key parameter here is the alpha parameter, which determines the threshold level of dissimilarity between frequent and infrequent reads for exclusion of infrequent reads. Note that we’re using a less conservative minsize threshold than the default of 8 because of the smaller size of our dataset.

* What is the effect on the number of sequences and size distribution of those sequences of varying the alpha parameter? You can get the size distribution by running

.. code-block:: bash 

	$ grep "^>" ​in.fasta​ | sed -e "s/size=\([^;]\)/\1/" | sort | uniq -c

Some researchers argue that denoising should be run at the level of the individual sample, not the dataset as a whole, because the frequency of reads is only meaningful relative to individual pools of amplicons. What do you think? If you’ve got the time, come back to this section later and do the following:

* Run ``vsearch --fastx_filter`` on each separate sample fastq file using a loop

* Run dereplication on each sample fasta separately using a loop

* Run denoising on each sample fasta separately using a loop

* Concatenate the results using sed and re-run dereplication using the following command: 

.. code-block:: bash 
	
	$ vsearch --derep_fulllength ​in.fasta --sizein --sizeout --relabel uniq --output out.fasta

* Compare the total unique read numbers and size distribution to the version produced earlier

------------------------------------------
Point error filtering
------------------------------------------

Filtering by length will remove sequences that have one or more PCR/sequencer-caused insertions or one or more deletions, however in some cases these errors may cancel one another out; or alternatively, PCR or sequencing may induce the equivalent of point mutations, where a single base is misread. Similarly, noncoding gene variants such as numts or pseudogenes may actually have point mutations in comparison to the ‘true’ region.

We can identify some point errors because they will alter the translation of the genetic code in such a way that it becomes meaningless - if the barcode region is a coding region, of course. The obvious error is the introduction of stop codons into the translation. By translating all of our sequences and checking for stop codons, we can easily reject these errors or variants. We use the script filtertranslate.py for this - check the helpfile by running:

.. code-block:: bash 

	$ filtertranslate.py -help

**footnote**

* Figure out what the command is to run it using automatic reading frame detection. Hint: check the usage line to figure out where some of the arguments go. Don’t forget, our samples are insects.

You may want to rename the automatically-named output.

* Have a look at the failed file. Go to an online amino acid translator and paste in a sequence. See what the translation looks like. Can you see the stop(s)?

Other ‘point errors’ are harder to spot. Some will not affect coding at all, which is impossible to distinguish from natural variation. The majority will affect coding, but again distinguishing these natural variation is very hard. This is actually something that some denoising algorithms attempt to do, broadly, but further work on this is ongoing.

------------------------------------------
Chimera filtering
------------------------------------------

You might have noticed that earlier, while running denoising, we noticed some chimeras. We left these in because we used the ``--ampout`` argument, so time to remove them.

We’re still using trusty VSEARCH:

.. code-block:: bash 

	$ vsearch --uchime3_denovo ​in.fasta​ --nonchimeras ​out.fasta

It’s that simple. And with that, we have a file that ideally contains only true biological sequences.

**footnote**















