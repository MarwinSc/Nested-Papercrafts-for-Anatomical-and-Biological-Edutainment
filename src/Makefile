all: lib.so

clean:
	rm -f *.o *.so *.html

lib.so:
	g++ -shared -Wl,-install_name,lib.so -o lib.so -fPIC C_demoCode.cpp
