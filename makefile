MATLAB_ROOT=/opt/MATLAB/R2021a
MEX=$(MATLAB_ROOT)/bin/mex

all: qcsp_detector.mexa64

qcsp_detector.mexa64: qcsp_detector.cpp lib/libqcsp_passed.a
	$(MEX) qcsp_detector.cpp -Llib -lqcsp_passed

clean:
	rm qcsp_detector.mexa64