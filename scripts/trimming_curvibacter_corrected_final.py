import subprocess
import os

def trimming(directory):
    print("[+] Current Working Directory : {}".format(os.getcwd()))
    files = []
    for file in os.listdir(directory):
        if file.endswith(".fastq.gz") == True:
            files.append(file)

    samples = []
    for file in files:
        if file.split("_R")[0] not in samples:
            samples.append(file.split("_R")[0])

    if (len(samples)+len(samples))%2 != 0:
        raise Exception("[-] number of files is odd and does not correspond to paired-end reads")

    for sample in samples:
        fw_file = sample + '_R1_001.fastq.gz'
        r_file = sample + '_R2_001.fastq.gz'
        if (str(fw_file.split("_")[0])+str(fw_file.split("_")[1])) != (str(r_file.split("_")[0])+str(r_file.split("_")[1])):
            raise Exception("[-] forward file sample {} does not match reverse file sample {}".format(str(fw_file.split("_")[0]+fw_file.split("_")[1]),r_file.split("-")))
        else:
            print("[*] starting trimming process of sample number {}".format(r_file.split("-")[0]))
            fw_paired_output = directory + "stringent_trimmed_files_final/fw_paired_out_" + fw_file
            fw_unpaired_output = directory + "stringent_trimmed_files_final/fw_unpaired_out_" + fw_file
            r_paired_output =  directory + "stringent_trimmed_files_final/r_paired_out_" + r_file
            r_unpaired_output = directory + "stringent_trimmed_files_final/r_unpaired_out_" + r_file
            clip_setting = "ILLUMINACLIP:TruSeq3_polyG.fa:2:30:10:2:keepBothReads"
            invoc = "java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar"
            cmd = "{} PE -threads 10 -phred33 {} {} {} {} {} {} {} CROP:60 HEADCROP:11 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(
                invoc,directory + fw_file, directory + r_file,fw_paired_output,fw_unpaired_output,r_paired_output,r_unpaired_output, clip_setting)
            print("[*] starting trimming ...")
            proc = subprocess.Popen(cmd, shell=True)
            returncode = proc.wait(timeout=2400)  # 40 Minutes
            if(returncode != 0):
                raise Exception("[-] error during trimming ..")
                return 1
            else:
                print("[+] DONE trimming of sample {}".format(sample))
    return 0

dirCurvi = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/AEP_TranscriptomicData_RAW/"
trimming(dirCurvi)
