import os
import subprocess

def get_files(input_directory:str):
	try:
		files = os.listdir(input_directory)
		samples = []
		for file in files:
			if file.startswith("dedup_") and file != "dedup_logfile.txt":
				# example file name:
				# dedup_forward_read_G1_S11_L001.fastq.gz
				sample = file.split("read_")[1].split(".")[0]
				if sample not in samples:
					samples.append(sample)
		if len(samples) <= 0:
			raise Exception("not enough samples in directory")
		else:
			return samples
	except Exception as e:
		print("[-] ERROR fetching samples, with exception: {}".format(e))
		
		
def star_aligning(star_output_dir:str, input_directory:str, fastq_files:list):
	try:
		logfile = star_output_dir + "star_rsem.log"
		with open(logfile,"w") as logfile:

			star_genome_dir = "/gpfs/project/lubec100/CuTrOmics/CurvibacterMonocolonization/genome_files/rsem_reference/"
			if os.path.isdir(star_genome_dir) == False:
				print("[-] STAR-RSEM genome directory does not exist!")
				logfile.write("ERROR:STAR-RSEM genome directory does not exist!\n")
				raise Exception("STAR-RSEM genome directory does not exist!")
			else:
				print("[*] Genome directory exist!")
				logfile.write("INFO:genome directory exist\n")

			star_genome_dir = star_genome_dir + "Curvibacter_reference"

			output_directory = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/star_rsem_output/"
			if os.path.isdir(output_directory) == False:
				print("[-] STAR-RSEM output directory does not exist!")
				logfile.write("ERROR:STAR-RSEM output directory does not exist!\n")
				raise Exception("STAR-RSEM output directory does not exist!")
			else:
				print("[*] Output directory exists!")
				logfile.write("INFO:output directory exist\n")

			logfile.write("INFO:starting alignment procedure\n")
			print("[*] Starting STAR alignment procedure")
			for sample in fastq_files:
				output_directory_sample = output_directory + sample
				fw_read = input_directory + "/dedup_forward_read_" + sample + ".fastq.gz"
				rev_read = input_directory + "/dedup_reverse_read_" + sample + ".fastq.gz"

				if os.path.isfile(fw_read) == False or os.path.isfile(rev_read) == False:
					print("[-] ERROR one of the read files for sample {} does not exist!".format(sample))
					logfile.write("ERROR:one of the read files for {} does not exist!\n".format(sample))
					raise Exception("one of the read files does not exist!")
				else:
					print("[+] working with sample: {}".format(sample))
					logfile.write("INFO:processing sample {}\n".format(sample))
					cmd = "rsem-calculate-expression -p {} --paired-end --star --star-gzipped-read-file {} {} {} {}".format(8,fw_read, rev_read, star_genome_dir, output_directory_sample)
					print("[*] Starting Popen process for rsem")
					logfile.write("INFO:executing {}\n".format(cmd))
					proc = subprocess.Popen(cmd, shell=True)
					returncode = proc.wait()
					if returncode != 0:
						print("[-] ERROR during subprocess call, returncode != 0")
						logfile.write("ERROR:returncode for subprocess call != 0\n")
						raise Exception("returncode != 0")
					else:
						print("[+] DONE with sample {}".format(sample))
						logfile.write("INFO:done mapping and count estimation for sample: {}\n".format(sample))
				print("DONE processing sample {}".format(sample))
		return 0
	except Exception as e:
		print("[-] ERROR during star rsem procedure with exception: {}".format(e))

dirCurvi = "/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/AEP_TranscriptomicData_RAW/stringent_trimmed_files_final/deduplicated_reads"

samples = get_files(dirCurvi)
star_aligning("/gpfs/project/lubec100/CuTrOmics/CurvibacterTranscriptome/star_rsem_output/", dirCurvi, samples)
