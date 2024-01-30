import os
import subprocess

print("[*] Starting quantification of reads with kallisto")

def kallisto_quant(directory):
    logfile = open(directory + '/kallisto_quant.txt','w')

    print("[*] START with dir : {}".format(directory))
    kallisto_index = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/AEP_Genome_Refseq/curvibacter_kallisto_index"
    gtf_file = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/AEP_Genome_Refseq/curvibacter_kallisto_index"

    if os.path.isfile(kallisto_index) == False:
        print("[-] ERROR can't find kallisto index file : {}".format(kallisto_index))
        logfile.write("[-] ERROR can't find kallisto index file : {}\n".format(kallisto_index))
        logfile.close()
        return 1

    files = []
    for file in os.listdir(directory):
        if file.endswith(".fastq.gz") == True:
            files.append(file)

    samples = []
    for file in files:
        if "dedup_forward" in str(file) or "dedup_reverse" in str(file):
            if file.split("_read_")[-1].split(".fastq.gz")[0] not in samples:
                samples.append(file.split("_read_")[-1].split(".fastq.gz")[0])

    print("[*] length samples: {} --> total files {}".format(len(samples),len(files)))
    for sample in samples:
        fw_file = directory + '/dedup_forward_read_'+ sample + '.fastq.gz'
        r_file = directory + '/dedup_reverse_read_' + sample + '.fastq.gz'

        if os.path.isfile(fw_file) == False or os.path.isfile(r_file) == False:
            print("[-] ERROR forward or reverse file not found for sample: {}".format(sample))
            logfile.write("[-] ERROR forward or reverse file not found for sample: {}\n".format(sample))
            logfile.close()
            raise FileNotFoundError("[-] ERROR either {} or {} not found ...".format(fw_file,r_file))
        elif os.path.isdir(sample) == True:
            print("[+] output already created for sample : {}".format(sample))
            logfile.write("[+] output already created for sample : {}\n".format(sample))
            print("[*] Skipping sample : {}".format(sample))
        else:
            print("[*] Starting kallisto procedure for sample : {}".format(sample))
            logfile.write("[*] Starting kallisto procedure for sample : {}\n".format(sample))
            cmd = "kallisto quant -i {} -o {} -t 8 -b 50 {} {}".format(kallisto_index,
                                                                directory + '/' + sample,
                                                                fw_file,
                                                                r_file)

            print("[*] {}".format(cmd))
            proc = subprocess.Popen(cmd, shell=True)
            returncode = proc.wait(timeout=7000)  # 66 Minutes
            if(returncode != 0):
                logfile.write("[-] error during kallisto - returncode for Popen != 0 ...\n")
                logfile.close()
                raise Exception("[-] error during kallisto ..")
            else:
                logfile.write("[+] DONE kallisto procedure of sample {}\n".format(sample))
                print("[+] DONE kallisto procedure of sample {}".format(sample))
    logfile.close()
    return 0

dirCurvi = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/AEP_TranscriptomicData_RAW/stringent_trimmed_files_final/deduplicated_reads"

kallisto_quant(dirCurvi)
