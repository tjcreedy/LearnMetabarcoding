==========================================
Read processing
==========================================

------------------------------------------
Introduction
------------------------------------------

The data for this section is already on the server, in the directory ​/AMM/resources/metabarcoding/​. This directory contains the results of every core step already, for your reference. The starting data is all in the ​0_rawsequences/ directory. Copy this over to your directory now

.. code-block:: bash
	
	$ cp -r /AMM/resources/metabarcoding/0_rawsequences ~/

------------------------------------------
FASTQ files
------------------------------------------

Let’s review a few basic points about sequence file structure and exploring these files using the command line. Change into the directory and list its contents, showing sizes:

.. code-block:: bash

	$ cd 0_rawsequences/

	$ ls -lh

We can see how many lines in each file using the word count function, specifying we want the number of lines:

.. code-block:: bash

	$ wc -l *.fastq

The ​:samp:`*.fastq` here means we want all of the fasta files in the directory. We could replace this with a single
file name.
We would like to review the structure of these files. We could print the whole file to the terminal using the cat command, but these files are over 35000 lines. We expect them to be the same structure all through the file, so the first 10 lines should be fine. Run the following command, replacing ​file with one of the ``file`` names.


.. code-block:: bash

	$ head -n 10 ​file

You will see the fastq format comprising header, sequence and quality scores. A useful point to note is the structure of the file header, specifically that it starts with “@D00”. If the structure of this file is completely new to you, take a few minutes to read the first section on the wikipedia page on fastq format: `https://en.wikipedia.org/wiki/FASTQ_format <https://en.wikipedia.org/wiki/FASTQ_format>`_

To get specific lines from a file, use the sed function:

.. code-block:: bash 

	$ sed -n ‘4,8p’ ​file # prints lines 4-8

* Use this to have a look at some files. If it doesn’t work - well you might have copied and pasted the command, which has formatted quotation marks. You need to type this out.

Note that the R1 and R2 files from the same library have the same read headers, apart from a 1 or 2 in the second part of the name. Reads with the same header were read from the same location on the sequencer, so they are assumed to be the forward and reverse read of the same fragment - these are called **​mate pairs​.** It’s important to ensure that both the forward and reverse read for each fragment are always kept present (“in sync”) for some future processes.

While we can divide the results of our wordcount function above by four to get the number of sequences, it is quicker to use ``​grep`` ​. This is a powerful tool that searches for text in a file using regular expressions, and we can use it to search for all the headers and return a count of them, as follows:

.. code-block:: bash

	$ grep -c "^@D00" ​file

The ​``^@D00`` is a regular expression - the ​^ character is a special symbol in regex meaning “the beginning of a line”. So this regex searches for lines beginning with ``​@D00`` We use more than just ``​@`` because it is also the symbol for a quality score of 31. Any sequences where the first base has a score of 31 will be counted twice - once for the header, once for the quality string, as they are two different lines starting with ``​@`` ​. Hence using multiple characters.

Like the ​wc -l function above, we can run grep on all of our files at once to get the total read numbers for each of our libraries:

.. code-block:: bash 

	$ grep -c "^@D00" *.fastq

We can see that we’re dealing with about 9000 reads per library.

* Do all the libraries have the same R1 and R2 read numbers?

------------------------------------------
Demultiplexing
------------------------------------------

There are many ways in which metabarcoding libraries may be sequenced. We are going to work here with a library preparation pipeline that involved indexing amplicons during initial PCR, such that each sample had a different 6-base index within a library. Once sequenced, we need to use these index sequences to separate out different samples. This process is called ​**demultiplexing​.** In this case, the amplicons were sequenced using paired-end sequencing, meaning the two ends of each fragment were read, working inwards. We have eight sequence files received from our sequencing facility, one for each read direction for each of four libraries.

Each of these libraries is actually three different metabarcoding samples. Each sample within each library was amplified with a different 6-base index. You can see these indices in the indices.txt file:

.. code-block:: 
	
	$ cat indices.txt

The ``​cat`` command simply outputs the entire contents of the files it is given, one after the other. Here we’re using it to just output the content of one file, but look out for this command later!

We can use ``grep`` to illustrate the location of these indices on the sequences. We take the output of the ``head`` command (from the previous section) to show only the first three sequences in a file, send this to the next function using the pipe (​|​) command. The next function is ``grep`` ​, which we use to search for an index sequence using regex. The ``|$`` at the end of the regex means “or the end of the line” - we use this to print every line, but only highlight the parts of interest.

.. code-block:: bash 

	$ head -n 12 ​file​ | grep -E "​INDEX​|$"

E.g.

.. code-block:: bash 

	$ head -n 12 Lib4_R1.fastq | grep -E "AACACC|$"

We’re going to be seeing a lot of pipes, so make sure you understand what they mean. Google “linux pipe tutorial” if you don’t understand.

When we did PCR, the index was part of the primers used, so that the reaction added this sequence to our amplicons when copying them. The primers used for this metabarcoding reaction were:

Forward (R1): **​CCNGAYATRGCNTTYCCNCG**

Reverse (R2): **​TANACYTCNGGRTGNCCRAARAAYCA**

Note the presence of non ATCG bases - these are ambiguities added to the primers to allow them to be less specific and bind to a greater variety of taxa.

* Use the grep command above to search for the primer sequence in a file.

	* Note that ambiguities (any base apart from ATCG) should be replaced with a ​. (full stop) which is a special regex symbol meaning “any character”.

	* Try writing the command yourself before looking at the answer in the footnote below.

	* Look at a few different libraries, both forward and reverse, both index and primer.

You’ll probably see that there are occasions where no index or primer is highlighted on a sequence. This means there was a sequencing error. Look closely and you’ll see that a base is missing or inserted, or just wrong.

We want to split each of these libraries up by index into a separate pair of files for each of the 12 samples, and remove the index sequences. We will do this using the tool cutadapt (`​https://cutadapt.readthedocs.io <​https://cutadapt.readthedocs.io​>`_). This versatile tool allows separating files by index and removal of these indices and primers. It can allow for some error in the index sequence, and can keep read pairs in sync (more on this later).

**Footnote was here in PDF** 

As ever with a new tool, first cast your eye over the help, either online or buy running: 

.. code-block:: bash 
	
	$ cutadapt --help

It’s quite long, but at least read the first section. It’s helpful to think about exactly what our data is and what we want to do:

* We have paired data

* Our indices are at the beginning of the reads

* We have multiple indices (cutadapt calls them ‘adapters’)

* We want to output a different file for each index


Cutadapt has settings for all of these situations. It will allow reading two files as input, and will ensure that pairs of reads in these files are kept in sync. Indices at the beginning of reads are specified using ``-g`` ( ``​-G`` for the second file of reads), and we can specify these multiple times, and give different adapters names. We also specify that the adapters are right at the beginning, with no gaps, using a ``^`` symbol. We can specify that we want output files depending on the combination of adapters found using the ``​-o`` and ``​-p`` options for the first and second files respectively.

Referring to the indices.txt file, we can now construct a command that demultiplexes our Lib1. To avoid a mess of files, I strongly suggest returning up to the parent directory and creating a new directory. Call this something appropriate.

This command assumes that you are in the directory containing the “0_rawsequences” directory and an empty directory called “1_demux”. Reminder: using the ​\ allows us to split the command over multiple lines. You can either type this and press enter afterwards, or you can just ignore it and continue typing the command at the beginning of the next line.

Run the following, then review what is printed to the terminal carefully. Make sure you understand what it’s saying:

.. code-block:: bash 

	$ cutadapt -g T4=^AAGAGG -g T9=^AATCGC -g T11=^AGCTAC \
	> -G T4=^AAGAGG -G T9=^AATCGC -G T11=^AGCTAC \
	> -o 1_demux/{name1}-{name2}_R1.fastq -p 1_demux/{name1}-{name2}_R2.fastq \ 
	> 0_rawsequences/Lib1_R1.fastq 0_rawsequences/Lib1_R2.fastq

* How many files do you expect to get out of this?

List the files in the demux directory, and run the grep command from the previous section to see the read numbers per file. More than you expected?

**Footnote again in PDF**

This is because the command has looked for all adapter combinations. This is nice security against errors. The sequencer has mistakenly associated some reads as the same fragment when they aren’t - they actually come from two different samples, hence some files with two different sample names. And in some cases, no index can be found on one or both of a paired read, probably due to a sequencing error. These are marked as unknown. Happily, all of these errors are in a distinct minority, and the majority of reads have been allocated to files for our samples.

If you add everything up, you’ll notice we’re missing some reads from our original files: these had no mate pair and were completely discarded.

* Construct three more cutadapt commands to demultiplex the other three libraries, placing the outputs into the same demux directory.

You should now have lots of files in that demux directory. It’s good practice to keep track of how demultiplexing performed: you could put the output of a grep command into a file to keep a record.

Let’s get rid of the files we don’t need. You’ve doubled the amount of storage you’re using - here the files aren’t very large but if you were doing this with a standard dataset, directories would fill up quickly. Navigate to the demux folder, very carefully copy the following command and run it. It works through the files, extracting the first and second sample name, then deletes the file if they don’t match. You do not need to type any ​#comments​, or add the extra spaces - this is just to make it clearer.

.. code-block:: bash 

	$ for f in *; do \                          # loop through files
	> s1=${f%-*} ; s2=${f%_*} ; s2=${s2#*-}; \  # extract sample names
	> if [ $s1 != $s2 ]; \ 			    # check if different
	> then rm $f; \ 			    # delete if different
	> else mv $f ${f#*-}; \			    # keep if same
	> fi; \ 				    # end if clause
	> done 					    # end loop 

An explanation for this code is below. This isn’t a crucial bioinformatics step, it’s just to tidy things up. You don’t need to understand everything about this command, although the loop syntax is going to crop up frequently.

Then also run:

.. code-block:: bash 

	$ rm unknown*

This will delete the files beginning with unknown. These contain the sequences for which no index was identified - we’re not interested in them.

**footnote again**

In the above command, the ``for`` part sets up a loop that works through each file name in the directory, using ``f`` to store the current name. It then uses something called parameter substitution to strip the names down to the first and second sample names, storing these in two variables. It then asks if the two sample names are different - if so, the current file is deleted, otherwise (i.e. they’re the same), the file is renamed (moved from one name to another), again using parameter substitution to strip out the unneeded parts of the name. You do not need to understand this.

Cutadapt has a lot more parameters for searching. It can be error-tolerant, allowing you to permit some mismatches between your indices and the sequences. Our adapters are all at least 3 bases different, so in theory we could allow a one-base difference to try and get more reads for our samples. The relevant option is ``-e`` ​, which is the maximum error rate (0-1), i.e. the total proportion of errors allowed in our indices. The default, 0.1, would allow 10% errors, but since our indices are only 6 bases this rounds to 0 errors allowed.

* Create two more new directories and try ``-e`` values that allow one error or two errors, putting the outputs into those two directories. Explore the read numbers using ``grep``.  

* How do the read numbers vary? Is being more error-tolerant sensible?

Feel free to explore more of the parameters, for example the minimum overlap parameter ``​-0`` 

* In this case, I believe the default settings are appropriate. What do you think?

------------------------------------------
Primer Removal
------------------------------------------

As well as demultiplexing, cutadapt removed the indices from our samples. You can check this using the grep command from the previous section. Alternatively, you can use the following command to explore the length distribution of the sequences. The ``sed`` command ​p​rints only every 4th line, starting at the 2nd. This is sent to a ``while`` loop, which reads each line and stores it in ``l`` ​. The loop outputs the number of characters in the ``l`` variable, one per line. These are then sorted into alphanumeric order, and then each unique number is counted to get the distribution. Run this on a pre- and post- demultiplexed file.
 
.. code-block:: bash

	$ sed -n '2~4p'​ file​ | while read l; do echo ${#l} ; done | sort | uniq -c

You should see that the average sequence length has reduced by 6.

We can also use cutadapt to remove primers. We cannot be certain that the amplified region of the primer sequence is exactly identical to that region on our source DNA, because primers do not always bind perfectly. So this region must be discarded.

This process is very similar to demultiplexing, except we only have one sequence to remove, rather than three, and we only want one output file for each input file. The power of cutadapt’s paired-file-aware approach is that again, we can filter out any mate pairs that don’t have both primers - this is definitely a mark of a sequencing error! First, create another new directory in our parent directory for the trimming output.

Making sure you’re in the parent directory, try and adapt our demultiplexing command to trim the primers given earlier from one of the demultiplexed file pairs. You will want to add the parameter ``--discard-untrimmed`` ​. We could have added this to demultiplexing to remove all “unknown” files as well.

Cutadapt is aware of ambiguous bases so it’s fine to use them as-is. The primers should have been consecutive with the indices, so now must be at the start of the reads: thus you can use ​^ to anchor the sequence as before. You don’t need to name the primer sequences (``XX=``), and you don’t need to use ``{name}``  in the output - the file name will do. Try running it, if it doesn’t work, check the answer.

Make sure to look over the output from cutadapt because this is very informative. You’ll notice now that some errors are being allowed, since these sequences are longer and so the default 10% allows 2 errors in these primers. Additionally you’ll notice that unlike with the indices, the length of sequence removed has varied slightly. We’ll come back to this.

We want to run this on all of our files, ideally without writing the command over and over. We can put this in a loop using bash. Since we need all of the unique samples, we first need to design a command for listing these:

.. code-block:: bash 
	
	$ ls 1_demux/* | cut -d_ -f1 | sort | uniq

This extracts the part of each name before the first ​_ and finds the unique ones. We can then use this as the basis of a loop:

.. code-block:: bash 

	$ while read s; do \
	> cutadapt -g ^CCNGAYATRGCNTTYCCNCG -G ^TANACYTCNGGRTGNCCRAARAAYCA \
	> -o 2_trimmed/${s}_R1.fastq -p 2_trimmed/${s}_R2.fastq --discard-untrimmed \ > 1_demux/${s}_R1.fastq 1_demux/${s}_R2.fastq; \
	> done < <(ls 1_demux/ | cut -d_ -f1 | sort | uniq)

The read command reads from the command piped in at the end, and the while command works through this bit-by-bit. The ``​$s`` refers to the sample name - note that ``​${s}`` is used where we want to add a ​_ immediately after, otherwise bash will look for a variable called ``$s_R1``.

Check your trimmed directory to make sure you have all of your files, and check back through the terminal output to make sure that you didn’t miss any errors. As always, review your read numbers.

**footnote** 

------------------------------------------
Quality settings
------------------------------------------

Because primers are a region on our sequence that we have some ​*a priori*​ knowledge about, this is a good opportunity for filtering sequences with errors.

* Try running primer trimming with the strictest settings, a 100% match in length and bases - is this sensible?

Many metabarcoding pipelines trim primers by just trimming a number of bases equal to the primer length off the beginning of each read. You can do this in cutadapt using the ​-u ​option, which you would need to do for each direction separately:

.. code-block:: bash 

	$ cutadapt -u ​n​ -o ​out.fastq​ -i ​in.fastq

* Try running this for the forward and reverse reads. If you look at the number of reads, this clearly retains more. What might be the downsides of doing it this way?

------------------------------------------
Pair Merging
------------------------------------------

These reads were sequenced such that 300bp of our fragment was read from each end. There should be 418bp of sequence between the primer pair.

* You know from previous sections how much index and primer sequence was trimmed from each end. How much sequence should we have left in each direction?

* How much should these sequences overlap?

This overlap between these is definitely enough to be able to confidently assemble each pair together, and each pair of the primer-trimmed files can now become one, by merging the mate pairs together by finding where the sequence overlaps.

An important consideration when pair merging is that, in general, sequence reads decline in quality along their length. This quality information is stored in the FASTQ file - it would be a lot of work to check all of our sequences manually, but there are some ways to summarise the quality scores. It’s good practice to do this with your reads at this point

One very popular option is to use the program FastQC. You can run it on a single file like this:

.. code-block:: bash 

	$ fastqc ​in.fastq

You can view the output of FastQC by using your FTP software (e.g. FileZilla) to access the server, navigating to the directory in which you ran fastqc and opening the html file with the same name as your file. There are lots of useful graphs, but the first one is key here.

If you don’t want to have to use your FTP software, one fun way of checking a file’s quality with output on the command line is fastqe:

.. code-block:: bash 

	$ fastqe --bin ​in.fastq

The ends away from the primer are the parts that the merging algorithm is attempting to match, and you’ll notice these tend to be lower quality. This leads to two considerations:

1. It is important to be aware of quality scores when merging, since erroneous sequence could cause the merging of reads that are not the same fragment, forming a chimera

2. However, merging can actually provide some validation of quality and improve quality scores for the overlap region if it matches.

As always, check out the helpfile for PEAR:

.. code-block:: bash 
	
	$ pear --help

Not too many options. Try running PEAR on a single pair of your primer-trimmed files from the previous section, after creating a new directory to hold merged files. The output should be the first part of the name - PEAR will add to this for the output files

.. code-block: bash 

	$ pear -f ​in_R1.fastq​ -r ​in_R2.fastq​ -o ​merged/out

PEAR gives us some really informative information on the terminal, make sure to review it. As you’ll have seen from the helpfile, there are many different thresholds that can be applied to the merging.

* Run your command again, applying a sensible minimum overlap based on what we calculated earlier.

* Run your command again, trying out some different p values and quality thresholds. The former is the threshold probability of an overlap being incorrect, the latter is the threshold score for trimming off low quality read ends.

I generally set the quality threshold for trimming low quality ends to 26, and the minimum overlap to a value around 20-30 bases less than the presumed overlap length. These sequences are pretty good quality so these settings won’t reject much more than the default in this case. Clean up your experimentation.

Like with cutadapt, we can run PEAR in a loop to apply it to all of our samples: again, we list our files, find our unique samples, then loop on these and use each sample name to run an iteration of PEAR. Here my command assumes your files are in a directory called 2_trimmed and you’re outputting to a directory called 3_merged - you should modify these as necessary:

.. code-block:: bash 

	$ while read l; do pear -f 2_trimmed/${l}_R1.fastq -r 2_trimmed/${l}_R2.fastq -o 3_merged/$l -q 26 -v 100; done < <(ls 2_trimmed/ | cut -d_ -f1 | sort | uniq)

Make sure to review those terminal outputs! Then list the contents of your merged directory, you’ll see files for the assembled, unassembled and discarded reads. We don’t need the latter two, so let’s clean up:

.. code-block:: bash 

	$ cd 3_merged
	$ rm *discarded* *unassembled* && rename -e "s/assembled\.//" *
	$ cd ../

The ``&&`` here runs both commands on this line one after the other. The ``rename`` command uses regular expressions to remove the “assembled followed by a full stop” part from any files with that in the name.

------------------------------------------
Pair concatenation
------------------------------------------

What happens if our fragment length is so much longer than our read length that it doesn’t overlap? For example, if we had sequenced these fragments using a 2x150bp metric instead of 2x300bp? Let’s simulate it.


Using the output from primer trimming, we can apply cutadapt to remove 150bp from the 3’ end of each read for a pair of files. First create a folder for this experiment (mine is called trimmed_150), then run these commands:

.. code-block:: bash 

	$ cutadapt -u -150 -o trimmed_150/T11_R1.fastq trimmed/T11_R1.fastq
	$ cutadapt -u -150 -o trimmed_150/T11_R2.fastq trimmed/T11_R2.fastq

* Try running PEAR on these two files.

We obviously can’t merge them, because there’s nothing left that overlaps. So instead we have to perform concatenation - stitching the forward and reverse reads together. This achieves the same aim as PEAR, converting into a single sequence rather than two independent reads.

The sequencer has introduced some length variation - some reads have recovered a little more data than others. PEAR took care of this automatically when it was forming a single read using the overlap, but concatenation is more basic. With non-overlapping data we need to concatenate like with like so that we generate a consistent (pseudo-)region of DNA. Otherwise we would generate spurious insertions/deletions when comparing our sequences.

In most cases, we would trim or discard reads such that the forward and reverse reads were each a fixed length, and then stitch each mate pair together to form a pseudo-locus. While there’s really a missing section of DNA in the middle, we have summarised our molecular information into a single fragment that is perfectly usable for many metabarcoding applications.

Let’s review the length distribution of our sequences to select a fixed length for each direction:

.. code-block:: bash 

	$ sed -n '2~4p' ​in.fastq​ | while read l; do echo ${#l} ; done | sort | uniq -c

Generally, we would select something around the central tendency, to retain as much data as possible.

* Select a length for each read. It does not need to be the same value.


We know the end with the primer is the accurate end, so we trim bases from the other end. We use cutadapt for this. The -l argument trims reads down to a value, and the -m argument specifies a minimum length. Set them as the same value:

.. code-block:: bash 

	$ cutadapt -l ​n​ -m ​n​ -o ​out.fastq​ ​in.fastq

Run this on your forward and reverse read file.

PEAR can stitch our mate pairs, and it reverse-complements the reverse reads for us.

.. code-block:: bash 

	$ pear -i -f ​in_R1.fastq​ -r ​in_R2.fastq -​ o ​outname

Oh dear. The problem is that we removed short reads without removing their mates. This gives us the opportunity to test using a tool for mate-pairing - i.e. making sure two files are in sync. We will use pairfq for this (`https://github/sestaton/pairfq <https://github/sestaton/pairfq>`_) :

.. code-block:: bash 

	$ pairfq makepairs -f ​in_R1.fastq -r ​in_R2.fastq -fp ​out_R1.fastq -rp ​out_R2.fastq -fs ​out_R1_unpaired.fastq​ -rs ​out_R2_unpaired.fastq

As always, use grep to check out file read numbers. Then run pear again to create the mate pairs.

* Are these concatenated sequences as reliable as our merged sequences? Why not?

You could at this point create a parallel dataset of concatenated short reads. Later, you can come back to these and work through the rest of the pipeline. Compare how these sequences behave in the future steps, particularly chimera filtering and OTU delimitation.

As an aside, if we have reads that are just too short to overlap, or too short to overlap well (e.g. < 10-20bp overlap), one option is to edit the reads such that the small missing region is padded with Ns. Reads that do overlap are merged if possible, or trimmed to be consecutive. Reads that are too short have N added to go up to the right length, and then the reads are stitched. This only applies to regions where we can be reasonably confident of a consistent, predictable length between primers. One issue would be the selection of an OTU delimitation method that took account of the ambiguous regions, otherwise the sequences would in general be more similar to one another than expected. For this reason this would only usually be done for projects with a small missing region.

------------------------------------------
Data concatenation
------------------------------------------

Now have complete set of sequence reads with extraneous sequence removed and read pairs brought together.

Metabarcoding works on all sequence reads from across the dataset to find OTUs, and it’s more efficient if we do the following steps with all reads compiled into one file. In theory this can be done at any point in the pipeline, but this seems like a convenient point.

We could simply use the ``cat`` command to concatenate all the files, but then we’d lose information about what reads come from what sample and this is important later. So we add the sample name to the header of each read prior to concatenation using the following command. Before you run this, make sure you’re in the directory containing your merged sequences.

.. code-block:: bash 

	$ for f in *; do \   # ​↓​ the space here should be included!
	> sed -e "s/\(^@D00.*\) .*$/\1;sample=${f%.*};/" $f \
	> >> ../mbc_concat.fastq; \
	> done

This command loops through all the files in the directory. For each file, the ``sed`` command using regular expression substitution to add a string to the end of each header line containing the sample name, and adds the entire file to the end of a single file in the parent directory called ``​mbc_concat.fastq​``. Regular expressions are a very powerful part of coding that I suggest you learn at some point if you don’t know about them already. There are plenty of tutorials on the internet.

Return to the parent directory, use ``grep`` to count the number of sequences in ``​mbc_concat.fastq`` ​and view the ``​head​`` of the file. We will use this file in the next practical.

















