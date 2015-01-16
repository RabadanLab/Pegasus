Pegasus
=======

Annotation and Prediction of Oncogenic Gene Fusions in RNAseq

REQUIREMENTS

Pegasus requires the following software packages to be installed on your system:
- JAVA
- PERL
- PYTHON 2.7 with the following libraries: 
      - numpy_1.8
      - pandas_0.12
      - scikit-learn_0.14
      
All the python libraries are included into the ANACONDA package.
( https://store.continuum.io/cshop/anaconda )

- Report from fusion detection tools. Version supported:
      - chimerascan 0.4.5
      - defuse 0.4.3
      - bellerophontes 0.4.0* 


INSTALL 
In order to install Pegasus, download, untar and unzip the Pegasus_dist.0.3.tgz tarbal into the install folder:

tar xvzf Pegasus_dist.0.3.tgz

Please note that the following files have to be downloaded separately and added to the install folder:

/path/to/Pegasus/hg19.fa
/path/to/Pegasus/hg19.fa.fai
/path/to/Pegasus/Homo_sapiens.GRCh37.60.chr.gtf

The dimension of these files is considerable and authors decided to provide them separately. You can find this files in Pegasus web site (See MANUAL.txt for details).


RUNNING PEGASUS ON TESTDATA

You can run Pegasus from the Pegasus folder with the following command:

./scripts/pegasus.pl [options]

the options are:

-s data_input_file: path to data input configuration file (mandatory). 
-d config_file: path to Pegasus configuration file (mandatory).
-k keepDB: 1 the database is cleaned, 0 otherwise (optional, default 1).
-p sge: 1 run on a SGE system, 0 otherwise (optional, default 0)
-o log_folder: path to a log folder (mandatory).

Example of running Pegasus:
> ./scripts/pegasus.pl -s data.txt -d pegasus_config.txt -o logs -p 0

In Pegasus tarball you also find an example of data configuration file (data.txt) and Pegasus configuration file (pegasus_config.txt).

data.txt and pegasus_config.txt are respectively an example data configuration file and an example Pegasus configuration file that you find in the Pegasus tar ball.

note: Pegasus is designed to be interrupted in any stage. Running Pegasus a second time will imply that the already implemented steps will be skipped. 
To clean Pegasus run it is highly recommended to remove the content of the logs and results folder.

Es. :
rm -r /path/to/Pegasus/logs/*
rm -r /path/to/Pegasus/results/*


For further details on manual configuration refer the MANUAL.txt file in the Pegasus tar ball.

