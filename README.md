    ProgramName:       /home/hxzk/script/pipeline_v3.0/bin/TS_to_genome.py
         Version:       V1.0
         Contact:       Liaoyw yiliu4234@163.com
    Program Date:       2024.4.25
          Modify:       -
     Description:       This is a software transforming enst position to genome position
           Usage:       python3 /home/hxzk/script/pipeline_v3.0/bin/TS_to_genome.py -g <gtf file> -b <bed file>  -o <outputfile>
                        ~/script/pipeline_v3.0/bin/TS_to_genome.py -g /yjdata/database/star/gtf/human/Homo_sapiens.GRCh38.99.gtf -b RNAIP_RNA_TestVsInput_peaks.xls -o peak_genome.txt
         Options:
            -h --help              Show this help message.
            -g --gtf   gtf file
            -b --bed  bed file
            -o --output outputfile
