ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

default: ExecuteMakeKStarSignalOnlyHistograms ExecuteFitKStarSignalOnlyShapes ExecuteMakeKStarWrongAssignmentHistograms ExecuteCombineKStarAssignments ExecuteMakeKStarSBHistograms ExecuteFitKStarSB

ExecuteMakeKStarSignalOnlyHistograms: MakeKStarSignalOnlyHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKStarSignalOnlyHistograms.cpp \
		-o ExecuteMakeKStarSignalOnlyHistograms \
		$(ROOTLIBS)

ExecuteFitKStarSignalOnlyShapes: FitKStarSignalOnlyShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSignalOnlyShapes.cpp \
		-o ExecuteFitKStarSignalOnlyShapes \
		$(ROOTLIBS)

ExecuteMakeKStarWrongAssignmentHistograms: MakeKStarWrongAssignmentHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKStarWrongAssignmentHistograms.cpp \
		-o ExecuteMakeKStarWrongAssignmentHistograms \
		$(ROOTLIBS)

ExecuteCombineKStarAssignments: CombineKStarAssignments.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		CombineKStarAssignments.cpp \
		-o ExecuteCombineKStarAssignments \
		$(ROOTLIBS)

ExecuteMakeKStarSBHistograms: MakeKStarSBHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKStarSBHistograms.cpp \
		-o ExecuteMakeKStarSBHistograms \
		$(ROOTLIBS)

ExecuteFitKStarSB: FitKStarSB.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSB.cpp \
		-o ExecuteFitKStarSB \
		$(ROOTLIBS)

ExecuteFitKStarSBFloatWidth: FitKStarSBFloatWidth.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSBFloatWidth.cpp \
		-o ExecuteFitKStarSBFloatWidth \
		$(ROOTLIBS)
