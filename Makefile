all: closure gauss fft prime lusplit

lusplit: src/utils.h src/lusplit.cc
	mpic++ -g -Wall -std=c++11 src/lusplit.cc -o sim_lusplit_parallel -fopenmp
prime: src/utils.h src/prime.cc
	mpic++ -g -Wall -std=c++11 src/prime.cc -o sim_prime_parallel -fopenmp
closure: src/utils.h src/closure.cc
	mpic++ -g -Wall -std=c++11 src/closure.cc -o sim_closure_parallel -fopenmp
gauss: src/utils.h src/gauss.cc
	mpic++ -g -Wall -std=c++11 src/gauss.cc -o sim_gauss_parallel -fopenmp
fft: src/utils.h src/fft.cc
	mpic++ -g -Wall -std=c++11 src/fft.cc -o sim_fft_parallel -fopenmp

