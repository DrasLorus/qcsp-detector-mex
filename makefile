# Makefile to produce the mex file on Linux. May work on MacOS, will not work on Windows.
# Use 'make MATLAB_ROOT=path/to/matlab/root' where 'path/to/matlab/root' is the folder
# which contains 'bin/matlab'.
MATLAB_ROOT=/opt/MATLAB/R2021a
MEX=$(MATLAB_ROOT)/bin/mex

HEADERS=\
	includes/CCorrAbsMaxGeneric.hpp \
	includes/CCorrelationEngineGeneric.hpp \
	includes/CDetectionStateGeneric.hpp \
	includes/CDetectorGenericInterface.hpp \
	includes/CDetectorSerialGeneric.hpp \
	includes/CIterativeAdderGeneric.hpp \
	includes/CNormGeneric.hpp \
	includes/CScoreAccumulatorGeneric.hpp \
	includes/CScoreProcessorGeneric.hpp


all: qcsp_detector.mexa64

qcsp_detector.mexa64: qcsp_detector.cpp lib/libqcsp_passed.a includes $(HEADERS)
	$(MEX) qcsp_detector.cpp -Llib -lqcsp_passed

lib/libqcsp_passed.a includes:
	@echo "  To acquire ./lib/libqcsp_passed.a and headers in ./includes/:\n  -- use QCSP_PASSED targets qcsp_passed and detector_export_public_headers.\n"
	@exit 1

clean:
	rm qcsp_detector.mexa64
