custom_logo: '/software/bi/biocore_tools/logo/biocore-logo.png'
custom_logo_url: 'https://biocore.crg.eu'
custom_logo_title: 'Bioinformatics Core Facility @ CRG'

extra_fn_clean_trim: 
   - '-trimmed'
   - '_log'
   - '_read1'
   - '_read2'

extra_fn_clean_trim: 
   - '-trimmed'


table_columns_visible:
        FastQC: 
                percent_gc: False    
                percent_duplicates: False
        HTseq Count:
                percent_assigned: False
module_order:
 - fastqc:
        name: 'FastQC (raw)'
        path_filters:
            - '*raw_fastqc.zip'
 - fastqc:
        name: 'FastQC (trimmed)'
        info: 'This section of the report shows FastQC results after adapter trimming.'
        target: ''
        path_filters:
            - '*-trimmed*_fastqc.zip'
 - star
 - htseq
 - read_counts 
