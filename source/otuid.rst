==============================================
Identifying OTUs and Generating Outputs
==============================================

------------------------------------------
Introduction
------------------------------------------

Now we have some OTUs in hand, we have two final tasks. We need to get some idea of what OTUs are present in which sample, and we want to try and identify our OTUs to some degree. Firstly, we’ll map our individual reads against the OTUs. Secondly, we’ll try and identify them against a global database using a very broad approach, then against a local database.

------------------------------------------
Mapping reads to OTUs
------------------------------------------

This step takes our trimmed, merged and quality filtered reads and tries to find which OTU each read belongs to. We take our reads prior to the more recent read filtering because this filtering may have removed many true but rare variants of our OTUs, or slight PCR or sequencing-induced error variants. These reads get matched to the closest OTU, but we set a similarity cut-off so that any reads that are major errors will not match to any OTU and be ignored.

We are going to use our 3% clustering OTUs for this at first, using the ``​VSEARCH --usearch_global`` command to search the reads against the OTUs, and output an OTU table.

.. code-block:: bash

	$ vsearch --usearch_global mbc_concat.fasta -db ​otus.fasta ​-id 0.97 -uc ​mapresults.uc

We use the ``-id 0.97`` parameter to set a 3% similarity cutoff, which is obviously appropriate for the 3% OTUs. The algorithm finds the best match in the database for each query sequence and output this information in a ``.uc`` file, which is a tabulated output. Have a look at this file using head or cat. Each line shows the hit statistics for each input read against its best match, if it has one. You should be able to see that each read has a ``sample=`` part to the name that we added earlier to identify its source. So all we need to do now is count up the number of hits for each OTU for each sample.

Thankfully, vsearch can do this for us. Run the command above again, but this time swap out the ``-uc out.uc`` o​ption for ``-otutabout ​out.tsv.`` Have a look at the table output using head or cat. Perfect!

The reads can be mapped to the CROP OTUs in a similar fashion, altering the -id value to fit the clustering threshold used in CROP. Swarm OTUs are not generated in quite as consistent a manner, but an appropriate threshold using vsearch would be ``​-id 0.99``.

For more accurate assignment of reads to swarm OTUs, you ​*could* try and do the following if you have time and the skills.

1. Use swarm’s ``-o`` option while generating OTUS to output a list of unique sequences making up each OTU

2. Use​ vsearch ``--search_exact`` to map the reads against the unique sequences

3. Link this mapping to the list of unique sequences somehow, then to the OTU sequences, and finally count up the number of sequences from each sample for each OTU.

------------------------------------------
Identification using GenBank and MEGAN
------------------------------------------

The most straightforward automated way of assigning taxonomy to OTUs is to BLAST the OTUs against the GenBank nucleotide database. The server has a subset of this database, and the blast executables.

Select one of your OTU files and blast it against nt using the following command:

.. code-block:: bash 
	
	$ blastn -db /AMM/blastdb/nt -query ​in.fasta​ -outfmt 5 -out ​out.xml​ -evalue 0.001

This command generates a very large XML file containing the full record of all the alignments BLAST has found for our OTUs. We need to transfer this to our computer for the next step, but to save bandwidth, let’s first compress the xml using zip:

.. code-block:: bash

	$ zip ​out.xml.zip in.xml

Using your FTP client, or whatever file transfer method you like, transfer the zipped XML file to your computer and extract it.

We will use the software MEGAN to parse the BLAST results and assign taxonomy to our OTUs. This is one of many algorithms out there for doing this - we introduce it here because it has a handy GUI interface and you’re probably bored of looking at the terminal.

Open up MEGAN. It usually takes a little while to open because it has to load the entire NCBI taxonomy into memory and display it. I bet you miss the efficiency of the terminal already.

MEGAN works through each OTU and find the location on the NCBI taxonomy tree for each GenBank sequence the OTU had a hit against. It then uses a Lowest Common Ancestor algorithm to estimate the most appropriate location on the tree for an OTU, and assigns the OTU the appropriate taxonomy.

Once MEGAN has opened and loaded the tree, you should see a very high-level cladogram of living organisms. This is the entire NCBI taxonomy. To map the OTU BLAST data onto this, you need to load the XML. Go to ​File > Import From BLAST​. In the Files tab of the window that appears, use the button to the right of the first box to browse to and select your XML file. MEGAN will automatically fill the third box. It should look something like this:

.. image:: megan_screenshot.png
	:align: center

Go to the ​LCA Params tab at the top. Here you will see the parameters that MEGAN uses when assigning taxonomy using its lowest common ancestor algorithm. For now we’ll leave these as default and just press Apply​.

Once MEGAN has finished you should see a reduced version of the taxonomy tree. It may not be very detailed: at the top bar, select Rank and choose Species. Have a look at the tree.

* Are all the OTUs Coleoptera?

Each circle on the tree is one or more OTUs that have been assigned to a node. The larger the circle, the more OTUs have been assigned to that node. If you click on a node, you’ll see two values. ​Assigned is the number of OTUs assigned to that node, ​Summed is the number of OTUs assigned to that node and all child nodes. If you ​right click on a node and click ​Inspect​, you can see more details about that node and the OTU(s) assigned to it, as well as all the BLAST information. The greyed out BLAST hits are those that aren’t taken into account in the LCA analysis.

You’ll notice that many OTUs have been assigned to internal nodes. Inspect some of these.

* Why do you think the algorithm has assigned them to internal nodes?

* Do you think algorithm is always correct?

To output the taxonomic assignment for all of the OTUs for use in analysis, we need to select all of the nodes. You can do this by going to ​Select > All Nodes​. Then go to ​File > Export > Text (CSV) Format​. For the Choose data to export: field, select ​readName_to_taxonPath​, click OK and select your output location. This generates a comma-separated table with the OTU name and full NCBI taxon path of the assigned node.

Use your FTP client to send this file to your directory on the server. We’ll come back to it later.Go to ​Options > Change LCA Parameters​. Let’s adjust this to only take into account close relatives by modifying the ​Min Percent Identity parameter to ignore any hits below, say, 80% similarity (i.e. ​Min Percent Identity​ = 0.80).

* How does this change things?

* Is this an appropriate selection for this dataset?

Return the ​Min Percent Identity to 0.0 and this time change the ​Top Percent value to 60 and the ​LCA Algorithm to weighted. This allows taking into account many more hits for the LCA algorithm, but weights them according to their score, which is more appropriate for shorter reads. This also may work better for an understudied community that doesn’t get many close hits.

* How does this change the taxonomy assignments?

------------------------------------------
Identification using curated databases
------------------------------------------

BLAST is not fundamentally designed as a taxonomic assignment tool, and MEGAN is forced to work with BLAST’s alignment and matching summary outputs to assign taxonomy. There are other tools out there that are specifically designed for the assignment of taxonomy to anonymous sequences. Rather than comparing a sequence against a database one-by-one, these methods use k-mer approaches to place a sequence within the whole set of references. Then, they use an approach broadly analogous to MEGAN’s LCA method to probabilistically assign taxonomy, generally assigning confidence scores to different taxonomic levels. Such methods include PROTAX, SINTAX, SPINGO and RDPclassifier. It should be noted that like much of the software available for metabarcoding, these tools are often written for use on 16S. It should always be considered whether other loci may not be treated as expected by such tools.

The other downside to the prior method is that GenBank is not necessarily authoritative - it is well known that many sequences available on GenBank are misidentified. This would not be an issue if we were working with a well-known taxon, but when our survey lineages are likely to be poorly covered in GenBank, yet we require a relatively detailed identification, searching against all of GenBank is likely to be less successful. Some researchers have taken sections of GenBank and curated the sequences, removing sequences that are poor quality or unlikely to be correct.

One issue with using these more advanced classification tools is that they often require quite specific reference database structures. Thankfully, many curated database authors have released their databases in these formats.

For CO1, the MIDORI database (Machida, 2017, doi: 10.1038/sdata.2017.27) is a curated version of GenBank’s CO1 sequences. The authors have taken the useful step of creating a server for searching sequences against MIDORI using three of the above classifiers. You can access the server here: `http://reference-midori.info/server.php <http://reference-midori.info/server.php>`_ 

You can see that you can select a program, paste or upload your sequences, and select a database and searching parameters. We don’t want to overload their server with redundant searches, so we’ve already done this step for you. We ran the 3% OTUs against the MIDORI CO1 database using RDP and SINTAX, and you can find the resulting files in ``/AMM/resources/metabarcoding/taxassign/`` under the different program names. We ran SPINGO too, but its outputs require more processing to be comparable, so we’ll just consider RDP and SINTAX. Copy the files to your directory. RDP outputs two files, the “hier_outfile” is a summary and the “usga_classified” is the individual OTU taxonomies.

To quickly get an idea of how many Coleoptera OTUs we have, run the following command on the SINTAX file, the RDP classified file, and the MEGAN output you uploaded:

.. code-block:: bash 

	$ grep -c “Coleoptera”

* Do the different assignment programs agree?

Download these files to your computer using your FTP client and open them up in a text editor or spreadsheet software. The exact format varies, but all they output broadly similar information: the name of the OTU, some taxonomy and a confidence for each taxonomic level. They are fairly intuitive. Compare the MEGAN, RDP and SINTAX classifications for some different OTUs.

* Which programs provide lower-level identifications?

* Are species level identifications likely to be accurate?

* What levels of confidence are given to the order level identifications? Might this be very conservative? Why?

* What other taxa do we apparently have? You will see that we have some obvious non-Coleoptera OTUs, but also some OTUs that have been assigned to other Insect orders. How consistent are these identifications between methods? Are we confident that these really are not Coleoptera?

Note that it’s perfectly feasible that there could have been non-Coleoptera Insect DNA in these samples.

------------------------------------------
Identification using a local database
------------------------------------------

It’s clear that global and curated databases are very useful for assigning broad taxonomy, but don’t do great identifying our OTUs to lower taxonomic levels for these sorts of samples. Usefully, we picked out some morphospecies and sequenced them separately using Sanger sequencing.

These sequences are in ``/AMM/references/canopy_Coleop_COX1_sI.fa``.

Let’s use BLAST to search our OTUs against this fasta of references. There’s no need to copy it to our directory. Run the following BLAST command:

.. code-block:: bash 

	$blastn-query​otus.fasta-​subject/AMM/references/canopy_Coleop_COX1_sI.fa-outfmt 6 -out ​out.txt -​ num_threads 1 -evalue 0.001 -perc_identity 97

Because we know that both our OTUs and our reference set are likely to all be closely related, we’re setting a strict ``-evalue`` and a threshold percentage identity so that we don’t simply get every OTU matching against every reference.

Use ``cat`` to view the output file. We’re using the standard blast tabulated output (``-​ outfmt 6``); you can find out what the columns refer to at `https://www.ncbi.nlm.nih.gov/books/NBK279684/ <https://www.ncbi.nlm.nih.gov/books/NBK279684/>`_. Rather obviously, the first refers to the query sequence, the second the subject, and the third the percent identity. You can see we’ve got some clear hits for some of our OTUs!

* Do any of our OTUs hit multiple different references? Why might this be?

* Do the OTUs matching Coleoptera references correspond to those assigned to Coleoptera using MEGAN?

* Does this give us more information about any of our OTUs compared with the global database search?

------------------------------------------
Mapping reads to references
------------------------------------------

If we were only interested in the species for which we have reference sequences, for example if we were monitoring for a single or set of known species, there isn’t really much need to generate anonymous OTUs. We could simply map our trimmed, merged and quality filtered reads against the references directly, to get counts of our reference species in each sample.

Try using both ``​blastn`` and ``vsearch --usearch_global`` to map the raw reads against the reference sequences, drawing on what you learned from using the previous commands.

Experiment with different thresholds and outputs.

* Experiment with different thresholds and outputs.

* Is one method clearly superior?

**Identifying contigs** 

You should by now have assembled some decent-length contigs - perhaps not complete, but near. However, we have no idea what morphospecies these genomes actually are!

The process of identifying these contigs is called baiting. We use short known sequences, in this case COX1 barcodes of our morphospecies, to identify the much larger complete mitogenomes. We can do this quite simply using BLAST.

The COX1 barcodes are in the /AMM/references/ directory. This is a relatively small dataset so there’s no need to bother copying it over and making an indexed BLAST database. Instead, we just BLAST against the file directly. Pick one of your assembly outputs and run BLAST:

.. code-block:: bash 

	$ blastn -query ​contigs.fa​ -subject /AMM/references/canopy_Coleop_COX1_sI.fa -outfmt 6 -perc_identity 95

You should very rapidly get a BLAST output table, which we can interrogate to see which contigs matched. We are looking for very high-identity matches here - these should be the exact same individuals. So >99%, over the entire length of the barcode.


Run the same BLAST command on all of your contigs files from the different assemblies. See if you can identify a good length contig for each of the five barcodes.

This is as far as we’ll take the assembly steps. The next stage after identifying these contigs is to find and annotate the mitochondrial genes, but this is beyond the scope of this workshop. If this is something you’re interested in learning about, some easy starting points for automated annotation are the MITOS web server `http://mitos.bioinf.uni-leipzig.de/index.py <http://mitos.bioinf.uni-leipzig.de/index.py>`_ and the MitoZ annotate script `https://github.com/linzhi2013/MitoZ <https://github.com/linzhi2013/MitoZ>`_. The latter works well but is still in development and can be tricky to set up and get to work properly.

We have generated complete, annotated versions of these novel mitogenomes, and these will be used in the next session to build phylogenetic trees.




	
