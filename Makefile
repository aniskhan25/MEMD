CXX ?= g++
CXXFLAGS ?= -std=c++11 -O2 -Wall -Wextra
CPPFLAGS ?= -Iinc -Isrc
SOURCES := $(wildcard src/MEMD_*.cc) $(wildcard src/UTY_*.cc)
TARGET := memd

.PHONY: all run demo-data benchmark clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ -o $@

demo-data:
	python3 scripts/generate_demo_data.py --n 1001 --ndim 16 --out data/demo_signal.csv

benchmark:
	python3 scripts/benchmark_cpu.py

run: $(TARGET) demo-data
	mkdir -p output
	./$(TARGET) data/demo_signal.csv output/memd

clean:
	rm -f $(TARGET)
	rm -f output/memd_*.txt
