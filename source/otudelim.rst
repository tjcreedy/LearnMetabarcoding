=======================================================
OTU Delimitation
=======================================================

------------------------------------------
Introduction
------------------------------------------

Now that we have our filtered reads, the next step is to find the Operational Taxonomic Units. The principle here is that the reads we have found are all the true biological sequences (aka haplotypes) present in our samples, however these encompass both intra- and inter-specific variation and we want data that is comparable to species-level ecological data. Therefore, we must try to collapse the intra-specific variation while maintaining inter-specific variation. We will try several different methods for doing this, and compare them.

------------------------------------------
3% clustering
------------------------------------------

The conceptually most simple method is OTU clustering, whereby clusters of similar species are designated as an OTUs, with the centroid of the cluster, i.e. the OTU sequence, being determined based on read abundance. This sort of OTU delimitation is implemented in VSEARCH

.. code-block:: bash 

	$ vsearch --cluster_size ​in.fasta​ --id 0.97 --centroids ​otus_3pc.fasta​ --sizein --relabel otu

The key parameter here is the ``--id``, which is the minimum distance (pairwise identity) between a sequence and the cluster centroid for that sequence to be part of a cluster. At ``-id 0.97``, this is called i.e. sequences in a cluster are no more than 3% different from the centroid. This is the theoretical threshold of intra-specific diversity being applied. 

* What happens if you modify this threshold? Try increasing it and decreasing it.

If you’re familiar with CD-HIT-EST, this is very similar. CD-HIT-EST uses a slightly different calculation of
dissimilarity, which can be used by including the parameter ``iddef 0`` 

------------------------------------------
Linkage-based clustering
------------------------------------------

SWARM is a slightly different clustering algorithm. Rather than using a threshold that applies to all clusters like vsearch’s clustering method, it uses a local linking threshold that is based on number of differences, rather than overall dissimilarity. Run swarm as follows:

.. code-block::bash 

	$ swarm -w ​out.fasta​ -d 1 -z ​in.fasta

* What happens if you run with higher ``-d`` values (they must be integers).

------------------------------------------
Bayesian clustering
------------------------------------------

A final method for clustering is CROP, which was originally designed for 16S and has not been updated in a while. Nonetheless, it uses an interesting algorithm based on assuming sequences are a gaussian mixture and models the process of clustering using Markov-chain Monte Carlo simulations.

CROP does not work on dereplicated reads, it needs the originals, so we must replicate them out again.

.. code-block:: bash 

	$ vsearch --rereplicate ​in.fasta​ --relabel repl --output ​out.fasta

We can then run this in crop. Note it creates several outputs using ``outname`` as a name base.

.. code-block:: bash 

	$ crop -b 40 -z 400 -s -r 0 -i ​in.fasta​ -o ​outname

The options ``-b``, ``-z`` and ``-r`` are for optimising the MCMC process. You can read more about these in their documentation (`​https://github.com/tingchenlab/CROP/wiki/THE-CROP-WIKI <​https://github.com/tingchenlab/CROP/wiki/THE-CROP-WIKI​>`_). You can leave these alone. Option ``​-s`` specifies that we want the equivalent of 3% clustering. We could change this to ``-g`` which is equivalent to 5%. You can alternatively supply values to ``-l`` and ``-u`` , for example ``​-l 1 -u 3`` , which are the lower and upper bounds of similarity levels.


CROP creates a bunch of extra files, you may want to ​mv​ them to their own folder. You could create a folder for each CROP run you do and make sure you’re in that directory when you run it.

* How do the number of OTUs differ between methods?

* For these sequences, 3% dissimilarity between two sequences is 0.03 * 418 = 12.54 differences. How does swarm perform if you use this value?

* Even if methods recover similar numbers of OTUs, are they necessarily recovering the same sequences?

For a quick look, we can use the following command to count the number of exactly identical sequences in all of the input files:

.. code-block:: bash 

	$ cat​ in1.fasta in2.fasta [in3.fasta]​ | \
	> perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' | \ 
	> sed -n '2~2p' | sort | uniq -c | grep -c "2"

.











