all: clean app

clean:
	find ./ -name "build" -type d -exec rm -rf {} +

test:
	cmake -S ./tests -B ./tests/build
	cmake --build ./tests/build

app:
	cmake -S . -B build
	cmake --build build

run:
	./build/DNA_Analyzer