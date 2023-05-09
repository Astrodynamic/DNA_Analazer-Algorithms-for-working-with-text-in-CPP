all: build

build: install

rebuild: clean install

install:
	cmake -S . -B ./build
	cmake --build ./build

unistall:
	find ./ -name "build" -type d -exec rm -rf {} +

clean: unistall
	rm -rf ./lib/*

tests:
	cmake -S ./test/ -B ./test/build
	cmake --build ./test/build
	./test/build/TEST

cppcheck: install
	@cd build/; make cppcheck

clang-format: install
	@cd build/; make clang-format

.PHONY: all build rebuild unistall clean cppcheck clang-format tests