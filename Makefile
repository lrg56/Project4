all: raytrace.c
	gcc raytrace.c -o raytrace

clean:
	rm -rf raytrace *~

test:
	./raytrace 400 400 input.json output.ppm
