CHARM_HOME ?= ~/charm/
STRUCTURE_PATH = ../utility/structures
OPTS = -g -I$(STRUCTURE_PATH) -DGROUPCACHE=0 -DDELAYLOCAL=0 -DCOUNT_INTRNS=0 -DDEBUG=0
CHARMC = $(CHARM_HOME)/bin/charmc $(OPTS)
LD_LIBS = -L$(STRUCTURE_PATH) -lTipsy

BINARY = simple
OBJS = Main.o Reader.o Particle.o BoundingBox.o

all: $(BINARY)

$(BINARY): $(OBJS)
	$(CHARMC) -language charm++ -o $(BINARY) $(OBJS) $(LD_LIBS)

proj: $(OBJS)
	$(CHARMC) -language charm++ -tracemode projections -o $(BINARY) $(OBJS) $(LD_LIBS)

$(BINARY).decl.h: $(BINARY).ci
	$(CHARMC) $(BINARY).ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o $(BINARY) charmrun

common.h: $(STRUCTURE_PATH)/Vector3D.h $(STRUCTURE_PATH)/SFC.h Utility.h

Main.o: Main.C $(BINARY).decl.h common.h Reader.h TreePiece.h BoundingBox.h BufferedVec.h TreeElement.h CacheManager.h Node.h Resumer.h Traverser.h Driver.h UserNode.h GravityVisitor.h DensityVisitor.h PressureVisitor.h CountVisitor.h
	$(CHARMC) -c $<

CacheManager.h: $(BINARY).decl.h

Reader.o: Reader.C Reader.h
	$(CHARMC) -c $<

Particle.o: Particle.C Particle.h
	$(CHARMC) -c $<

BoundingBox.o: BoundingBox.C BoundingBox.h
	$(CHARMC) -c $<


test: all
	./charmrun ./simple -f ../inputgen/100k.tipsy +p3 ++ppn 3 +pemap 1-3 +commap 0 ++local

