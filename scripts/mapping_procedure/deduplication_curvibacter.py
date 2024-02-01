import subprocess
import os

def deduplication(directory:str):

    if os.path.isdir(directory) == False:
        print("[-] Directory does not exist!")
        return 1

    if os.path.isdir(directory+"/deduplicated_reads") == False:
        print("[-] Output directory does not exist!")
        return 1

    logfileName = directory + "/deduplicated_reads/dedup_logfile.txt"
    logfile = open(logfileName,'w')
    print("[+] Current Working Directory : {}".format(directory))
    files = []
    for file in os.listdir(directory):
        if file.endswith(".fastq.gz") == True:
            if file.startswith("fw_paired_out") or file.startswith("r_paired_out"):
                files.append(file)

    samples = []
    for file in files:
        sample = file.split("_R")[0].split("_out_")[-1]
        print("[*] New Sample : {}".format(sample))
        if sample not in samples:
            samples.append(sample)

    for sample in samples:
        forward_read = directory + "/fw_paired_out_" + sample + "_R1_001.fastq.gz"
        reverse_read = directory + "/r_paired_out_" + sample + '_R2_001.fastq.gz'

        if os.path.isfile(forward_read) == False or os.path.isfile(reverse_read) == False:
            logfile.write("error reverse or forward read does not exist! aborting mission! \n")
            print("[-] Forward {} or reverse file does not exist! Abort Deduplication ...".format(forward_read,reverse_read))
            return 1

        output_filename = "deduplicated_" + sample + '_both_reads.fastq.gz'
        print("[*] Working on forward {} and reverse {} files - trying to create deduplication file {}".format(forward_read, reverse_read, output_filename))
        output = directory + "/deduplicated_reads/" + output_filename

        if os.path.isfile(output):
            print("[*] Output file: {} already exist, skipping file: {}".format(output,file))
            logfile.write("output file {} already exists ...  skipping ...\n".format(output))
        else:
            print("[*] Starting Deduplication Command")

            cmd = "~/bbmap/dedupe.sh threads=8 in1={} in2={} out={} ac=f".format(forward_read,reverse_read,output)
            proc = subprocess.Popen(cmd, shell=True)
            returncode = proc.wait(timeout=2400) # 40 Minutes
            if(returncode != 0):
                raise Exception("[-] error during deduplication ..")
            else:
                print("[+] DONE deduplication of sample {}".format(sample))
                logfile.write("done deduplication of sample {}\n".format(sample))
            print("[*] Starting Reformatting Command")

            dedup_fw_read = directory + "/deduplicated_reads/dedup_forward_read_" + sample + '.fastq.gz'
            dedup_rv_read = directory + "/deduplicated_reads/dedup_reverse_read_" + sample + '.fastq.gz'


            cmd = "~/bbmap/reformat.sh threads=8 in={} out1={} out2={}".format(output,dedup_fw_read,dedup_rv_read)
            proc = subprocess.Popen(cmd, shell=True)
            returncode = proc.wait(timeout=2400) # 40 Minutes
            if(returncode != 0):
                raise Exception("[-] error during reformat..")
            else:
                print("[+] DONE reformat of sample {}".format(sample))
                logfile.write("done reformat of sample {}\n".format(sample))
    logfile.close()
    return 0

dirCurvi = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/AEP_TranscriptomicData_RAW/stringent_trimmed_files_final"
deduplication(dirCurvi)
