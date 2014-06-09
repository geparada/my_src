bash pre_processing.sh &&
bash process_smallRNAS.sh > >(tee -a log.txt) 2>&1
