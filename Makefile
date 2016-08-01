all: brat_bw build_bw trim acgt-count remove-dupl
brat_bw: BRAT_bw_nova.cpp
	g++ -O3 -w -o brat_bw BRAT_bw_nova.cpp
build_bw: Build_bw_nova.cpp
	g++ -O3 -w -o build_bw Build_bw_nova.cpp
trim: Filter_and_trim.cpp
	g++ -O3 -w -o trim Filter_and_trim.cpp
acgt-count: Coverage_ACGT_nova.cpp
	g++ -O3 -w -o acgt-count Coverage_ACGT_nova.cpp
remove-dupl: Remove_copy_duplicates.cpp
	g++ -O3 -w -o remove-dupl Remove_copy_duplicates.cpp

clean: brat_bw build_bw trim acgt-count remove-dupl
	rm brat_bw build_bw trim acgt-count remove-dupl

