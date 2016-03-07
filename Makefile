CXX = icc
MPICXX = mpicxx
SOURCE = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCE))
EXECUTABLE = 2D_Unsteady_Diffusion
#CXXFLAGS =-O0 -g -Wno-deprecated -Wall $(VTK_CPPFLAGS)
CXXFLAGS =-O3 -openmp -Wno-deprecated -Wall $(VTK_CPPFLAGS) -I/usr/include/vtk-5.8
LDFLAGS = $(VTK_LDFLAGS)
LIBS = -lvtkCommon -lvtkFiltering -lvtkGraphics -lvtkIO -lvtkRendering -lvtkWidgets

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(MPICXX) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE) $(LIBS)
	@echo DONE!

-include $(OBJECTS:.o=.d)

%.o: %.cpp
	$(MPICXX) -c $(CXXFLAGS) $*.cpp -o $*.o
	@$(MPICXX) -MM $(CXXFLAGS) $*.cpp > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -rf *.o *.d *.vtk $(EXECUTABLE)
	@echo ALL CLEANED UP!

rebuild:
	make clean
	make
