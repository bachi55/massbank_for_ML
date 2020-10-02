===========================
Data for Method Development
===========================

To develop our method we require *annotated* tandem mass spectra (:math:`\text{MS}^2`) associated with retention times
(RTs). Our method assumes that in practice multiple (:math:`\text{MS}^2`, RT)-tuples are measured simultaneously within
one experiment. Therefore, our training dataset needs to contains sub-sets of compounds for which the :math:`\text{MS}^2`
and RTs have been measured under the same experimental conditions, such as the liquid chromatographic (LC) system,
mass spectra ionization, etc.

Massbank
========

`Massbank <https://github.com/MassBank/MassBank-data/releases/tag/2020.09>`_ is a community driven public repository
containing annotated (tandem) mass spectra. For most Massbank records a rich set of meta-data is provided. These include
RTs, experimental conditions and structure information on measured molecule. The database (DB) is divided into
contributors, e.g. a particular laboratory, and furthermore into acquisitions, e.g. different experimental setups within
one laboratory.

Extracting (:math:`\text{MS}^2`, RT)-tuples from Massbank
---------------------------------------------------------

Using the `massbank2db <https://github.com/bachi55/massbank2db>`_ (v0.4.2) we build an SQLite database (DB) containing
Massbank records of release 2020.09. Details of resulting DB structure, can be found in the massbank2db package
description. We will review the main points here again.

Grouping Records based on their Experimental Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only *within* each contributor and acquisition, e.g Eawag with EA and EQ, the records are grouped, as we assume that
(especially) RTs are harder to compare between laboratories, i.e. contributors, and that each contributor has a reason
to introduce multiple acquisition (prefixes). A summary of all contributors and acquisition (Prefix of ID) `is provided
in the Massbank repository <https://github.com/MassBank/MassBank-data/blob/main/List_of_Contributors_Prefixes_and_Projects.md>`_.

Furthermore, records are grouped by their mass spectrometry (MS) and liquid chromatographic (LC) configuration. A side
note here: We are only concerned with LC systems, but not with gas chromatography (GC). Massbank contains records using
the latter analysis technique, which is ignore in our experiments.

The result of the first level of grouping as described above is illustrated in the following example:

    (AU, Athens_Univ)
        -> AU_000  ~ positive, LC configuration 1, MS instrument 1

        -> AU_001  ~ positive, LC configuration 2, MS instrument 1

        -> AU_002  ~ negative, LC configuration 3, MS instrument 1

        -> AU_003  ~ negative, LC configuration 4, MS instrument 2

        -> ...

    Here, for example, AU_001 is a new dataset identifier allowing to access all records uploaded by Athens_Univ with
    acquisition prefix AU where LC configuration 2 was used together with MS instrument 1 operating in positive ionization
    mode.

Filtering out "incomplete" Records
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For our method development we consider only Massbank records providing sufficient meta data. That means, we require:

- Records are required to have a RT reported
- Records are required to have the LC configuration specified
- The measured compound (molecule) needs to be in PubChem (copy from June 2019) either linked by CID or InChIKey

  - allows usage of PubChem's standardized SMILES (structure) representation
  - ensures same SMILES standardization as (later) for molecular candidates used

- Records are required to have an :math:`\text{MS}^2`
- Remove datasets, e.g. AU_001, with less than 50 (:math:`\text{MS}^2`, RT)-tuples (after the previous filtering steps)

Ensuring a set of complete, accurate and sufficiently sized datasets motivates the proposed data filtering.

Basic Dataset Information
-------------------------

(Check `basic_dataset_statistics.csv <msms_rt_ssvm/datasets/basic_dataset_statistics.csv>`_)

.. csv-table:: Datasets after grouping with basic information.
    :file: basic_dataset_statistics.csv
    :header-rows: 1
