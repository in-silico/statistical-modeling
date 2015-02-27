
FLAGS := -lgplib -larmadillo -std=c++11

all: mvgauss.test gpreg.test

mvgauss.test: test_mvgauss.cpp
	g++ test_mvgauss.cpp -o mvgauss.test $(FLAGS)

gpreg.test: test_gpreg.cpp
	g++ test_gpreg.cpp -o gpreg.test $(FLAGS)
