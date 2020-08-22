FROM biocorecrg/centos-perlbrew-pyenv3-java:centos7

# File Author / Maintainer
MAINTAINER Toni Hermoso Pulido <toni.hermoso@crg.eu>
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 
ARG FASTQC_VERSION=0.11.9
ARG STAR_VERSION=2.7.3a
ARG QUALIMAP_VERSION=2.2.1
ARG MULTIQC_VERSION=1.8
ARG SAMTOOLS_VERSION=1.10
ARG TOOL_MULTIQC_VERSION=1.1
ARG HTSEQ_VERSION=0.11.1
ARG BEDTOOLS_VERSION=2.29.2
ARG HTSLIB_VERSION=1.10.2

#INSTALLING BEDTOOLS
RUN bash -c 'curl -k -L https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static.binary > /usr/local/bin/bedtools'
RUN chmod +x /usr/local/bin/bedtools

#INSTALLING FASTQC
RUN bash -c 'curl -k -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip > fastqc.zip'
RUN unzip fastqc.zip; chmod 775 FastQC/fastqc; ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc

# Installing STAR
RUN bash -c 'curl -k -L https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz > STAR.tar.gz'
RUN tar -zvxf STAR.tar.gz
RUN cp STAR-${STAR_VERSION}/bin/Linux_x86_64/* /usr/local/bin/
RUN rm STAR.tar.gz

# Installing QUALIMAP
RUN bash -c 'curl -k -L https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v${QUALIMAP_VERSION}.zip > qualimap.zip'
RUN unzip qualimap.zip
RUN rm qualimap.zip

RUN ln -s $PWD/qualimap_v${QUALIMAP_VERSION}/qualimap /usr/local/bin/

# Installing MULTIQC // Latest dev version is much better. 
# RUN pip install -Iv https://github.com/ewels/MultiQC/archive/v${MULTIQC_VERSION}.tar.gz 
RUN pip install --upgrade --force-reinstall git+https://github.com/ewels/MultiQC.git

# Installing samtools
RUN yum install -y xz-devel.x86_64
RUN bash -c 'curl -k -L https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 > samtools.tar.bz2'
RUN tar -jvxf samtools.tar.bz2
RUN cd samtools-${SAMTOOLS_VERSION}; ./configure; make; make install; cd ../
RUN rm samtools.tar.bz2

#INSTALLING htslib
RUN bash -c 'curl -k -L https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 > htslib.tar.bz2'
RUN tar -jvxf htslib.tar.bz2
RUN cd htslib-${HTSLIB_VERSION}; ./configure; make; make install; cd ../
RUN rm htslib.tar.bz2

#Installing kenttools
RUN yum install -y libpng12
RUN bash -c 'curl -k -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig > /usr/local/bin/bedGraphToBigWig'
RUN bash -c 'curl -k -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedSort > /usr/local/bin/bedSort'
RUN bash -c 'chmod +x /usr/local/bin/bed*'

#Adding perl script for improving multiQC report
RUN bash -c 'curl -k -L  https://github.com/CRG-CNAG/make_tool_desc_for_multiqc/archive/v${TOOL_MULTIQC_VERSION}.tar.gz > tool_ver.tar.gz'
RUN tar -zvxf tool_ver.tar.gz
RUN mv make_tool_desc_for_multiqc-${TOOL_MULTIQC_VERSION}/make_tool_desc_for_multiqc.pl /usr/local/bin/ 
RUN chmod +x /usr/local/bin/make_tool_desc_for_multiqc.pl
RUN rm -fr make_tool_desc_for_multiqc-* v${TOOL_MULTIQC_VERSION}.tar.gz 
RUN pip install --upgrade pip
RUN pip install cython 
RUN pip install pysam 

#Installing HTSEQ
RUN pip install numpy matplotlib
RUN pip install HTSEQ==${HTSEQ_VERSION}

#cleaning
RUN yum clean all
RUN rm -fr *.zip *.gz *.bz2 
RUN rm -fr STAR-* bedtools2  samtools-*
RUN rm -rf /var/cache/yum
