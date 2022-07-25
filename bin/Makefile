all : stats results zip

stats :
	RUN=$(shell head -1 RUN); outdir=$${RUN}-stats; mkdir -p $$outdir; cd AppResults; collect-stats.sh $$RUN ../SAMPLES ../$$outdir

results :
	RUN=$(shell head -1 RUN); cd AppResults; collect-results.sh $$RUN ../SAMPLES; cp -r consensus ..

zip : 
	RUN=$(shell head -1 RUN); ZIP=$${RUN}-stats.zip; \
	cp consensus/$${RUN}-results.* $${RUN}-stats; \
	zip -r $${ZIP} $${RUN}-stats
