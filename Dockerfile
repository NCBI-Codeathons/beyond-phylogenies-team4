######Pipeline dependencies & apps:#########

FROM amazonlinux:latest

# metadata
LABEL base.image="amazonlinux:latest"
LABEL dockerfile.version="1"
LABEL software="FlaCo pipeline"
LABEL software.version="v1"
LABEL description="Florida Covid transmission cluster determination and dissection"
LABEL website="https://github.com/NCBI-Codeathons/beyond-phylogenies-team4/tree/main/bin"
LABEL license="None"
LABEL maintainer="Team4"
LABEL maintainer.email="foobar@gmail.com"

RUN yum -y install openssl-devel libxml2-devel libcurl-devel harfbuzz-devel fribidi-devel freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel wget tar gcc-7.3.1 gcc-gfortran g++ libstdc++-devel fontconfig-devel bzip2-devel pcre2-devel make which perl git java libffi libffi-devel

#Get Python3
RUN wget https://www.python.org/ftp/python/3.9.6/Python-3.9.6.tgz && tar -xvf Python-3.9.6.tgz && cd Python-3.9.6/ && ./configure --enable-optimizations && make altinstall && cd ..

#conda, pangolin, usher and nextflow (python dependencies)
RUN wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh && \
        bash Anaconda3-2022.05-Linux-x86_64.sh -b && \
        anaconda="/root/anaconda3/bin/conda" && \
        eval "$(${anaconda} shell.bash hook)" && \
	/usr/local/bin/python3.9 -m pip install --upgrade pip && \
        /usr/local/bin/python3.9 -m pip install biopython numpy setuptools && \
        conda deactivate && \
        git clone https://github.com/yatisht/usher.git && \
        cd usher/ && \
        eval "$(${anaconda} shell.bash hook)" && \
        conda create -n usher-env && \
        conda activate usher-env && \
        conda config --add channels defaults && \
        conda config --add channels bioconda && \
        conda config --add channels conda-forge && \
        conda install usher && conda deactivate && cd .. && \
        eval "$(${anaconda} shell.bash hook)" && \
        git clone https://github.com/cov-lineages/pangolin.git && \
        cd pangolin/ && conda env create -f environment.yml && \
        conda activate pangolin && \
        /usr/local/bin/python3.9 -m pip install . && conda deactivate && cd .. && \
        eval "$(${anaconda} shell.bash hook)" && \
        wget -qO- https://get.nextflow.io | bash && \
        conda deactivate


#Get latest R from CRAN (don't trust default download, it may be <4.1.0)
RUN wget https://cran.r-project.org/src/base/R-4/R-4.2.1.tar.gz && \
	tar xfz R-4.2.1.tar.gz && cd R-4.2.1/ && \
	./configure --with-readline=no --with-x=no && \
	make && \
	./bin/R --slave -e 'install.packages(c("optparse","phytools","devtools","rlist", "tidyverse","BiocManager","ggplot2","ape","phangorn","shiny","dplyr","reshape2"),repos="https://cloud.r-project.org",dependecies=TRUE)' && cd .. 

ENV PATH=${PATH}:"/R-4.2.1/bin"

#iqtree2:
RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.1.2/iqtree-2.1.2-Linux.tar.gz && tar xfz iqtree-2.1.2-Linux.tar.gz

ENV PATH=${PATH}:"/iqtree-2.1.2-Linux/bin"

#MAFFT:
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.490-gcc_fc6.x86_64.rpm && rpm -Uvh mafft-7.490-gcc_fc6.x86_64.rpm 

#BLAST:
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz && tar xfz ncbi-blast-2.13.0+-x64-linux.tar.gz && chmod 777 /ncbi-blast-2.13.0+/bin/*

ENV PATH=${PATH}:"/ncbi-blast-2.13.0+/bin"

#FastTree
RUN wget http://www.microbesonline.org/fasttree/FastTree && chmod a+x FastTree 

#Minimap2
RUN git clone https://github.com/lh3/minimap2 && cd minimap2 && make && cd ..

ENV PATH=${PATH}:"/minimap2"

#ViralMSA
RUN wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py" && sed -i '1{s@.*@#! /usr/local/bin/python3.9@}' ViralMSA.py && chmod a+x ViralMSA.py

#Flaco
RUN git clone https://github.com/salemilab/flaco.git && chmod 777 /flaco/bin/*

ENV PATH=${PATH}:"/flaco/bin"

#RaxML
RUN git clone https://github.com/stamatak/standard-RAxML.git && cd standard-RAxML/ && make -f Makefile.SSE3.gcc && chmod 777 raxmlHPC-SSE3 && cd ..

ENV PATH=${PATH}:"/standard-RAxML"

#Extra script
RUN wget https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/src/mask_alignment_using_vcf.py && mv mask_alignment_using_vcf.py mask_aln_using_vcf.py && \
	sed -i '1{s@.*@#! /usr/local/bin/python3.9@}' mask_aln_using_vcf.py && chmod a+x mask_aln_using_vcf.py

# And the GitHub scripts!
RUN git clone https://github.com/NCBI-Codeathons/beyond-phylogenies-team4.git && chmod 777 /beyond-phylogenies-team4/bin/*

ENV PATH=${PATH}:"/beyond-phylogenies-team4/bin"

ENV PATH=${PATH}:"/"
WORKDIR /

COPY . .
