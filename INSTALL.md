# How to install CNVScan

## Common Process

   1. Download the compressed library.

   2. Uncompress the file into a directory which you want.

      $ tar xvfz CNVScan-2.0.1.tgz


   3. Download and install BioPython (http://biopython.org/wiki/Download)
      If you have pip3 installed it is easy: 
      $ pip3 install biopython

      If you don't have pip3 installed you can download it:
      * Ubuntu : sudo apt-get install python3-pip
      * Centos:
        # First command requires you to have enabled EPEL for CentOS7
        sudo yum install python34-setuptools
        sudo easy_install pip

   4. Download GFF parser for Biopython using pip3:
      $ sudo pip3 install bcbio-gff



Install CNVScan

   1. Go to the folder where you just unzipped your file

      $ cd CNVScan-2.0.1


   2. Launch CNVScan

$ python3 pCNVScan2.py <options>

## Common Questions

    1. I got the error: ImportError: No module named 'BCBio'

       You are missing the GFF parser for Biopython (http://biopython.org/wiki/GFF_Parsing). Check the step 4 in installation process.

    2. I got the error: ImportError: No module named 'Bio'

       You are missing BioPython for python3. Check if you are using the good pip for installation. See Installation step 3.


