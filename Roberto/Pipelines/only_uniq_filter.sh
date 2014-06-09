# Before you proceed, make sure your alignment file is sorted by the read name column (or any column with a unique identifier for each read)
# Separate the read name column with "cut". The parameter of -f should be the column of the read name. This is 1 for .sam, while it is 4 for .bed
cut -f1 example.sam > reads.txt

# Count the number of occurances (i.e. mapped locations) of each read name 
uniq -c reads.txt > count.txt

# Add the count to the last column of each line. "join" merges the two files based on column 1 of the first file and column 2 of the second file. If the original file is in .bed format, use "-1 4" instead.
join  -1 1 -2 2 example.sam count.txt > example_with_count.txt

# Keep only lines with a count of 1
awk '$NF == 1' example_with_count.txt > example_unique.txt

# Remove the last column
awk '{delim = ""; for (i=1;i<=NF-1;i++) {printf delim "%s", $i; delim = OFS}; printf "\n"}' example_unique.txt > example_unique.sam
