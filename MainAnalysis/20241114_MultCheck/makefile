Execute: plotManyZ.cpp
	g++ plotManyZ.cpp -o Execute \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o


ExecuteHistFillerGen: HistogramFillerGen.cpp
	g++ HistogramFillerGen.cpp -o ExecuteHistFillerGen \
		`root-config --libs --cflags` \
		-I$(ProjectBase)/CommonCode/include \
		$(ProjectBase)/CommonCode/library/*.o