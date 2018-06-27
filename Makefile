all:
	mkdir -p build; cd build; cmake ../ItrAs; make -j8; cd ..

test:
	./build/tools/ItrAs -in ItrAs/tests/Staphylococcus_aureus_AllPathsCor.fa -mink 25 -step 10 -maxk 95 -out TestResults

help:
	./build/tools/ItrAs --help

clean:
	rm -rf build/ TestResults/