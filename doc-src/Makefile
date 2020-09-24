OBJECTS:=${subst ../plots, imgs, $(patsubst %.pdf,%.svg,$(wildcard ../plots/stability/*.pdf))}

$(info $(OBJECTS))

all: $(OBJECTS)


imgs/stability/%.svg: ../plots/stability/%.pdf
	inkscape -l $@ $<
