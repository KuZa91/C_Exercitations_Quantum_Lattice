OPTGCC   = -std=c++11 -O3 -w -mtune=native -fprefetch-loop-arrays -funroll-all-loops -ffast-math
OPT0     = -std=c++11 -O0 -g


all: INFO

INFO:	
	@echo
	@echo "Compilation options are:"
	@echo
	@echo "     out     : optimized output example code"
	@echo "     int     : optimized integration example code"
	@echo "     cls     : optimized classes example code"
	@echo "     mat     : optimized matrix example code"
	@echo
	@echo "     out0   	: Non-Optimized output example code"
	@echo "     int0  	: Non-Optimized integration example code"
	@echo "     cls0  	: Non-Optimized classes example code"
	@echo "     mat0 	: Non-Optimized matrix example code"
	@echo

out:	ScalFieldHMC.c
	g++ $(OPTGCC) ScalFieldHMC.c -o ScalFieldHMC.run

clean:
	@rm -f -r ./DataSets/*.dat
